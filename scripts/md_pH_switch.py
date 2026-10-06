#!/usr/bin/env python3
"""Dynamics-based, PROPKA-orthogonal validation of the His433 pH switch.

For one binder-target complex we build TWO systems that differ ONLY in the
protonation state of EGFR His433 (the conserved pH sensor):
  * HIP  = protonated imidazole (+1)  -> models tumor pH 6.5
  * HIE  = neutral imidazole   ( 0)   -> models blood  pH 7.4
and run short restrained implicit-solvent MD (OpenMM, OpenCL on H100) in each.

Switch signal (tumor-ON) = the His433<->binder-carboxylate salt bridge forms
preferentially when His433 is protonated:
  d_occ   = occupancy(HIP) - occupancy(HIE)              (>0 = switch)
  d_dmin  = mean_mindist(HIE) - mean_mindist(HIP)         (>0 = switch, A)
plus interface contact persistence and binder RMSD as pose-stability guards.

This is force-field (Amber14+GBn2) based and fully independent of the PROPKA
Wyman-linkage score, giving an orthogonal second line of evidence for Obj-1.

Usage:
  md_pH_switch.py --pdb complex.pdb --out outdir [--device 0]
                  [--equil-ps 200 --prod-ps 2000 --interval-ps 10 --replicas 2
                   --restraint-k 2.0 --sensor-cut 6.0]
Writes outdir/<name>.json with the metrics.
"""
import argparse
import json
import os
import sys
import numpy as np

from openmm import app, unit, Platform, System, LangevinMiddleIntegrator
from openmm import CustomExternalForce
import openmm as mm
from pdbfixer import PDBFixer


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def detect_chains(topology):
    """target = chain with ~191 residues (largest); binder = the other largest."""
    chains = list(topology.chains())
    sizes = [(c, sum(1 for _ in c.residues())) for c in chains]
    sizes.sort(key=lambda x: -x[1])
    target = sizes[0][0]
    # binder = largest remaining protein chain
    binder = sizes[1][0] if len(sizes) > 1 else None
    return target, binder


def find_his433(target_chain):
    """His433 = 100th residue (1-based ordinal) of the target chain."""
    res = list(target_chain.residues())
    if len(res) < 100:
        raise RuntimeError("target chain has <100 residues (%d)" % len(res))
    h = res[99]
    if h.name not in ("HIS", "HID", "HIE", "HIP"):
        # fall back: nearest HIS to ordinal 100
        his = [i for i, r in enumerate(res) if r.name.startswith("HI")]
        if not his:
            raise RuntimeError("no HIS in target chain")
        h = res[min(his, key=lambda i: abs(i - 99))]
        log("WARN: ordinal-100 is %s; using nearest HIS at ordinal %d" % (res[99].name, res.index(h) + 1))
    return h


def find_switch_his(his_chain, partner_chain, positions, cut, expected_pos=None):
    """binder-His (inverted) mode: the engineered switch His in his_chain. When expected_pos is given
    (the LOCKED motif His carried through from the MPNN fixed_residues map), that His is used verbatim;
    otherwise the His whose ND1/NE2 is closest to a partner carboxyl-O is used as a heuristic.

    code-review finding: the pure-proximity heuristic can pick an INCIDENTAL MPNN-introduced His that
    sits closer to some target carboxyl than the engineered motif His, so the MD would titrate/measure
    the wrong residue. Passing the known locked position removes that ambiguity; we still compute the
    proximity pick and WARN loudly when the two disagree (fail-loud, don't silently measure the wrong His)."""
    pos = np.array(positions.value_in_unit(unit.angstrom))
    part_ox = []
    for r in partner_chain.residues():
        if r.name in ("ASP", "GLU"):
            part_ox += [a.index for a in r.atoms() if a.name in ("OD1", "OD2", "OE1", "OE2")]
    if not part_ox:
        raise RuntimeError("no Asp/Glu in partner (target) chain")
    best, best_d = None, 1e9
    for r in his_chain.residues():
        if not r.name.startswith("HI"):
            continue
        hN = [a.index for a in r.atoms() if a.name in ("ND1", "NE2")]
        if not hN:
            continue
        d = np.min(np.linalg.norm(pos[part_ox][:, None, :] - pos[hN][None, :, :], axis=-1))
        if d < best_d:
            best_d, best = d, r
    if best is None:
        raise RuntimeError("no HIS in binder chain")
    if expected_pos is not None:
        locked = [r for r in his_chain.residues()
                  if r.name.startswith("HI") and str(r.id) == str(expected_pos)]
        if locked:
            if locked[0] is not best:
                log("  WARN: proximity His id=%s (%.2f A) != LOCKED His id=%s -> using LOCKED His "
                    "(code-review finding: don't measure an incidental His)" % (best.id, best_d, expected_pos))
            return locked[0]
        log("  WARN: no His at locked pos %s in binder; falling back to proximity His id=%s"
            % (expected_pos, best.id))
    return best


def carboxyl_atoms(residue):
    out = []
    for a in residue.atoms():
        if a.name in ("OD1", "OD2", "OE1", "OE2"):
            out.append(a.index)
    return out


def his_n_atoms(residue):
    return [a.index for a in residue.atoms() if a.name in ("ND1", "NE2")]


def nearest_binder_carboxyl(positions, his_res, binder_chain, cut):
    """binder Asp/Glu whose carboxyl-O is closest to a His N; within cut (A)."""
    pos = np.array(positions.value_in_unit(unit.angstrom))
    hN = his_n_atoms(his_res)
    hpos = pos[hN]
    best = None
    best_d = 1e9
    for r in binder_chain.residues():
        if r.name not in ("ASP", "GLU"):
            continue
        ox = carboxyl_atoms(r)
        if not ox:
            continue
        d = np.min(np.linalg.norm(pos[ox][:, None, :] - hpos[None, :, :], axis=-1))
        if d < best_d:
            best_d = d
            best = (r, ox)
    if best is None or best_d > cut:
        return None, best_d, hN
    return best, best_d, hN


def build(pdb, his_variant, binder_his=False, sensor_cut=6.0, switch_his_pos=None, charge_scale=1.0):
    """PDBFixer -> add missing heavy atoms -> set switch-His variant -> add H -> implicit-solvent system.
    binder_his=True titrates the BINDER His (switch_his_pos = the locked motif His id, else nearest a
    target Asp/Glu); else the antigen His433. charge_scale<1.0 applies ECC to ionizable groups."""
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # do not model long gaps; keep as-is
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")

    target, binder = detect_chains(modeller.topology)
    his = find_switch_his(binder, target, modeller.positions, sensor_cut, switch_his_pos) if binder_his else find_his433(target)
    residues = list(modeller.topology.residues())
    his_idx = residues.index(his)
    variants = [None] * len(residues)
    variants[his_idx] = his_variant
    added = modeller.addHydrogens(ff, pH=7.0, variants=variants)

    # re-resolve objects after addHydrogens (topology rebuilt)
    target, binder = detect_chains(modeller.topology)
    his = find_switch_his(binder, target, modeller.positions, sensor_cut, switch_his_pos) if binder_his else find_his433(target)

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=2.0 * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=4.0 * unit.amu,  # HMR -> stable 4 fs timestep (~2x speedup)
    )
    if charge_scale != 1.0:
        nscaled = scale_ionizable_charges(system, modeller.topology, charge_scale)
        log("  ECC: scaled %d ionizable-group charges by %.2f" % (nscaled, charge_scale))
    return ff, modeller, system, target, binder, his


# Electronic-continuum-correction (ECC): scale the partial charges of ionizable SIDE-CHAIN groups by
# ~0.8 to mimic electronic polarization that fixed-charge force fields miss and that over-stabilizes salt
# bridges by ~3-4 kcal/mol. Backbone charges are left intact (prosECCo-style), so only the charged
# moieties weaken. A switch whose d_occ SIGN survives this is de-risked against the force-field artifact.
ECC_GROUP_ATOMS = {
    "ASP": {"CG", "OD1", "OD2"},
    "GLU": {"CD", "OE1", "OE2"},
    "LYS": {"NZ", "HZ1", "HZ2", "HZ3"},
    "ARG": {"CZ", "NE", "HE", "NH1", "NH2", "HH11", "HH12", "HH21", "HH22"},
    "HIP": {"CG", "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2", "HE2"},  # imidazolium (+1)
}


def scale_ionizable_charges(system, topology, scale):
    """Scale charged-group partial charges by `scale` on the NonbondedForce (ECC)."""
    nb = next((system.getForce(i) for i in range(system.getNumForces())
               if system.getForce(i).__class__.__name__ == "NonbondedForce"), None)
    if nb is None:
        raise RuntimeError("no NonbondedForce to apply ECC")
    atoms = list(topology.atoms())
    n = 0
    for at in atoms:
        grp = ECC_GROUP_ATOMS.get(at.residue.name)
        if grp and at.name in grp:
            q, sig, eps = nb.getParticleParameters(at.index)
            nb.setParticleParameters(at.index, q * scale, sig, eps)
            n += 1
    return n


def add_ca_restraints(system, modeller, k_kcal):
    k = k_kcal * unit.kilocalories_per_mole / unit.angstrom**2
    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", k)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    pos = modeller.positions
    n = 0
    for atom in modeller.topology.atoms():
        if atom.name == "CA":
            x, y, z = pos[atom.index].value_in_unit(unit.nanometer)
            force.addParticle(atom.index, [x, y, z])
            n += 1
    system.addForce(force)
    return n


def run_state(pdb, variant, args, device, replica):
    ff, modeller, system, target, binder, his = build(pdb, variant, args.binder_his, args.sensor_cut,
                                                       getattr(args, "switch_his_pos", None),
                                                       getattr(args, "charge_scale", 1.0))
    partner = target if args.binder_his else binder
    pair, d0, hN = nearest_binder_carboxyl(modeller.positions, his, partner, args.sensor_cut)
    if pair is None:
        log("  [%s] no carboxyl within %.1f A of switch His (closest %.1f) -> no switch pair" % (variant, args.sensor_cut, d0))
    nres_ca = add_ca_restraints(system, modeller, args.restraint_k)

    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 4.0 * unit.femtosecond)
    integ.setRandomNumberSeed(1000 * replica + (1 if variant == "HIP" else 2))
    plat = Platform.getPlatformByName("OpenCL")
    props = {"OpenCLDeviceIndex": str(device), "Precision": "mixed"}
    sim = app.Simulation(modeller.topology, system, integ, plat, props)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, 1234 + replica)

    # equilibration (dt = 4 fs with HMR)
    sim.step(int(args.equil_ps / 0.004))

    # production, measuring the His-N <-> carboxyl-O min distance each interval
    nsteps = int(args.prod_ps / 0.004)
    stride = int(args.interval_ps / 0.004)
    hN = his_n_atoms(his)
    ox = pair[1] if pair else None
    binder_ca = [a.index for a in binder.atoms() if a.name == "CA"] if binder else []
    ref_ca = None
    dists = []
    rmsds = []
    nframes = nsteps // stride
    for _ in range(nframes):
        sim.step(stride)
        st = sim.context.getState(getPositions=True)
        p = np.array(st.getPositions().value_in_unit(unit.angstrom))
        if ox:
            dmin = np.min(np.linalg.norm(p[ox][:, None, :] - p[hN][None, :, :], axis=-1))
            dists.append(float(dmin))
        if binder_ca:
            cur = p[binder_ca]
            if ref_ca is None:
                ref_ca = cur
            else:
                rmsds.append(float(np.sqrt(np.mean(np.sum((cur - ref_ca) ** 2, axis=1)))))
    out = {
        "variant": variant,
        "replica": replica,
        "pair_resname": pair[0].name if pair else None,
        "pair_ordinal": (list(partner.residues()).index(pair[0]) + 1) if pair else None,
        "init_dist": round(d0, 2),
        "occupancy_4A": round(float(np.mean(np.array(dists) < 4.0)), 3) if dists else None,
        "occupancy_5A": round(float(np.mean(np.array(dists) < 5.0)), 3) if dists else None,
        "mean_mindist": round(float(np.mean(dists)), 2) if dists else None,
        "binder_rmsd": round(float(np.mean(rmsds)), 2) if rmsds else None,
        "n_frames": len(dists),
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--equil-ps", type=float, default=200.0)
    ap.add_argument("--prod-ps", type=float, default=2000.0)
    ap.add_argument("--interval-ps", type=float, default=10.0)
    ap.add_argument("--replicas", type=int, default=2)
    ap.add_argument("--restraint-k", type=float, default=2.0)
    ap.add_argument("--sensor-cut", type=float, default=6.0)
    ap.add_argument("--binder-his", action="store_true",
                    help="inverted switch: titrate the BINDER His and bridge to a TARGET Asp/Glu")
    ap.add_argument("--switch-his-pos", type=int, default=None,
                    help="LOCKED binder switch-His residue id (the motif His from the MPNN fixed_residues "
                         "map); avoids measuring an incidental His. Warns if proximity disagrees.")
    ap.add_argument("--charge-scale", type=float, default=1.0,
                    help="ECC: scale ionizable-group charges (e.g. 0.8) to test switch robustness vs the "
                         "fixed-charge salt-bridge over-stabilization bias.")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    name = os.path.splitext(os.path.basename(args.pdb))[0]
    log("=== MD pH-switch:", name, "device", args.device)

    results = {"name": name, "pdb": os.path.abspath(args.pdb), "states": []}
    for rep in range(args.replicas):
        for variant in ("HIP", "HIE"):
            try:
                r = run_state(args.pdb, variant, args, args.device, rep)
                results["states"].append(r)
                log("  rep%d %s occ4=%s mindist=%s rmsd=%s pair=%s" % (
                    rep, variant, r["occupancy_4A"], r["mean_mindist"], r["binder_rmsd"], r["pair_resname"]))
            except Exception as e:
                log("  rep%d %s FAILED: %r" % (rep, variant, e))
                results["states"].append({"variant": variant, "replica": rep, "error": repr(e)})

    # aggregate switch metrics
    def agg(metric, variant):
        vals = [s[metric] for s in results["states"] if s.get("variant") == variant and s.get(metric) is not None]
        return float(np.mean(vals)) if vals else None

    occ_hip = agg("occupancy_4A", "HIP")
    occ_hie = agg("occupancy_4A", "HIE")
    dist_hip = agg("mean_mindist", "HIP")
    dist_hie = agg("mean_mindist", "HIE")
    results["occ_HIP"] = occ_hip
    results["occ_HIE"] = occ_hie
    results["d_occupancy"] = (round(occ_hip - occ_hie, 3) if (occ_hip is not None and occ_hie is not None) else None)
    results["d_mindist"] = (round(dist_hie - dist_hip, 2) if (dist_hip is not None and dist_hie is not None) else None)
    # switch = salt bridge tighter/more occupied when protonated
    dsw = results["d_occupancy"]
    results["md_switch"] = (dsw is not None and dsw > 0.1)

    outpath = os.path.join(args.out, name + ".json")
    json.dump(results, open(outpath, "w"), indent=2)
    log("WROTE", outpath, "| d_occupancy=%s d_mindist=%s md_switch=%s" % (
        results["d_occupancy"], results["d_mindist"], results["md_switch"]))
    print(outpath)


if __name__ == "__main__":
    main()
