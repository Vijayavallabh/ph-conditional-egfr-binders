#!/usr/bin/env python3
"""pH-dependent MM-GBSA binding free energy — a 2nd, functionally-direct Obj-1 signal.

Our four co-folders (AF2-ig/Boltz-2/Chai-1/AF3) are pH-BLIND: they score the
neutral-state fold only, so ALL pH-selectivity evidence has rested on the MD
salt-bridge *occupancy* (md_pH_switch.py) + the PROPKA/pKAI pKa legs. This script
adds an orthogonal, thermodynamic readout: the *binding free energy itself*,
computed in the two protonation states of the conserved sensor His433.

For one binder-target complex we build two systems differing ONLY in His433:
  * HIP = protonated imidazole (+1)  -> tumor pH 6.5
  * HIE = neutral imidazole   ( 0)   -> blood  pH 7.4
run the SAME restrained implicit-solvent MD as md_pH_switch.py, and on each
production snapshot compute a single-trajectory MM-GBSA interaction free energy
  dG_bind = <E(complex) - E(binder) - E(target)>      (Amber14 + GBn2, kcal/mol)
(no configurational entropy — standard single-trajectory MM-GBSA; it cancels to
first order in the difference below). The pH-selective BINDING signal is
  dd_bind = dG_bind(HIP) - dG_bind(HIE)   (< 0 = binds more tightly when
  His433 is protonated = tumor-ON), the functional complement of d_occupancy.

Because this is force-field based and independent of PROPKA's Wyman linkage and
of the co-folders, dd_bind is a genuinely orthogonal third line of pH evidence.
Trust the SIGN (consistent with our sign-consensus policy), not the magnitude:
single-trajectory GBSA magnitudes are noisy, but the sign is driven by the
His(+)...Glu(-) electrostatic term GB handles well.

Usage:
  md_pH_bind.py --pdb complex.pdb --out outdir [--device 0]
                [--equil-ps 200 --prod-ps 2000 --interval-ps 10 --replicas 2
                 --restraint-k 2.0 --sensor-cut 6.0]
Writes outdir/<name>.json with dG_bind_HIP/HIE, dd_bind_kcal, and (cross-check)
the occupancy switch metrics.
"""
import argparse
import json
import os
import sys

import numpy as np
from openmm import app, unit, Platform, LangevinMiddleIntegrator, VerletIntegrator
import openmm as mm

# reuse the validated build + sensor helpers (single source of truth)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from md_pH_switch import (  # noqa: E402
    build, add_ca_restraints, detect_chains, his_n_atoms, nearest_binder_carboxyl,
)

KJ2KCAL = 1.0 / 4.184


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def create_energy_system(ff, topology):
    """Unrestrained GBn2 system for single-trajectory MM-GBSA energy evaluation.
    NoCutoff is the correct MM-GBSA setting: the 191-residue EGFR target spans >2 nm, so a
    finite nonbonded cutoff would truncate long-range Coulomb/GB pairs from the complex energy
    without a matching term in the separated binder-only / target-only subsystems, biasing the
    interaction energy. (Sampling MD keeps its cutoff; only the energy eval uses NoCutoff.)"""
    return ff.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        hydrogenMass=None,
    )


def make_subsystem(ff, complex_modeller, keep):
    """Copy the (H-added, His-variant) complex modeller, delete the other chain,
    return (energy_system, context-less) for the kept chain only."""
    m = app.Modeller(complex_modeller.topology, complex_modeller.positions)
    tgt, bnd = detect_chains(m.topology)
    drop = tgt if keep == "binder" else bnd
    m.delete(list(drop.residues()))
    sysE = create_energy_system(ff, m.topology)
    return sysE, m.topology.getNumAtoms()


def energy_context(system, device):
    plat = Platform.getPlatformByName("OpenCL")
    props = {"OpenCLDeviceIndex": str(device), "Precision": "mixed"}
    return mm.Context(system, VerletIntegrator(1.0 * unit.femtosecond), plat, props)


def run_state(pdb, variant, args, device, replica):
    # NOTE: the SAMPLING MD below (build + CA restraints + LangevinMiddle 310 K / 4 fs-HMR +
    # seed scheme + minimize + equil + production stride) is an INTENTIONAL MIRROR of
    # md_pH_switch.run_state so the occupancy axis and this binding axis sample the SAME restrained
    # trajectory. Any change to that protocol MUST be mirrored in both files (it reuses build /
    # add_ca_restraints / nearest_binder_carboxyl from md_pH_switch). The only new logic here is the
    # per-frame single-trajectory MM-GBSA energy evaluation.
    ff, modeller, sampling_system, target, binder, his = build(pdb, variant)
    pair, d0, _ = nearest_binder_carboxyl(modeller.positions, his, binder, args.sensor_cut)

    # atom-index maps from the complex (order preserved when we strip a chain)
    tgt0, bnd0 = detect_chains(modeller.topology)
    binder_cidx = np.array([a.index for a in bnd0.atoms()], dtype=int)
    target_cidx = np.array([a.index for a in tgt0.atoms()], dtype=int)

    # energy systems (unrestrained): complex / binder-only / target-only
    sysE_cpx = create_energy_system(ff, modeller.topology)
    sysE_bnd, n_bnd = make_subsystem(ff, modeller, "binder")
    sysE_tgt, n_tgt = make_subsystem(ff, modeller, "target")
    assert n_bnd == len(binder_cidx), (n_bnd, len(binder_cidx))
    assert n_tgt == len(target_cidx), (n_tgt, len(target_cidx))
    ctx_cpx = energy_context(sysE_cpx, device)
    ctx_bnd = energy_context(sysE_bnd, device)
    ctx_tgt = energy_context(sysE_tgt, device)

    def pot(ctx, coords_nm):
        ctx.setPositions(unit.Quantity(coords_nm, unit.nanometer))
        return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)

    # sampling MD (restrained, identical protocol to md_pH_switch.py)
    add_ca_restraints(sampling_system, modeller, args.restraint_k)
    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 4.0 * unit.femtosecond)
    integ.setRandomNumberSeed(1000 * replica + (1 if variant == "HIP" else 2))
    plat = Platform.getPlatformByName("OpenCL")
    props = {"OpenCLDeviceIndex": str(device), "Precision": "mixed"}
    sim = app.Simulation(modeller.topology, sampling_system, integ, plat, props)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, 1234 + replica)
    sim.step(int(args.equil_ps / 0.004))

    nsteps = int(args.prod_ps / 0.004)
    stride = int(args.interval_ps / 0.004)
    hN = his_n_atoms(his)
    ox = pair[1] if pair else None
    binder_ca = [a.index for a in binder.atoms() if a.name == "CA"] if binder else []
    ref_ca = None
    dists, rmsds, dG = [], [], []
    for _ in range(nsteps // stride):
        sim.step(stride)
        p = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        p = np.asarray(p)
        # MM-GBSA interaction energy on this snapshot. Skip non-finite frames: an exploded/
        # overlapping snapshot -> inf, inf-inf=nan, which np.mean would propagate and json.dump
        # would write as a bare `NaN` token (invalid strict JSON + silently scored as leg=0 by the
        # consumer). A non-finite frame is dropped (and if ALL frames are bad, dG stays empty -> None).
        e_cpx = pot(ctx_cpx, p)
        e_bnd = pot(ctx_bnd, p[binder_cidx])
        e_tgt = pot(ctx_tgt, p[target_cidx])
        de = (e_cpx - e_bnd - e_tgt) * KJ2KCAL
        if np.isfinite(de):
            dG.append(de)
        # occupancy cross-check (angstrom)
        if ox:
            pa = p * 10.0
            dmin = float(np.min(np.linalg.norm(pa[ox][:, None, :] - pa[hN][None, :, :], axis=-1)))
            dists.append(dmin)
        if binder_ca:
            cur = p[binder_ca] * 10.0
            if ref_ca is None:
                ref_ca = cur
            else:
                rmsds.append(float(np.sqrt(np.mean(np.sum((cur - ref_ca) ** 2, axis=1)))))

    # free the 4 live OpenCL contexts now (1 Simulation + 3 energy) rather than at scope-exit GC,
    # to avoid ~2x transient peak GPU memory on the shared <=3-GPU box (a `for c in (...): del c`
    # loop only rebinds the loop var and frees nothing).
    del ctx_cpx, ctx_bnd, ctx_tgt, sim
    return {
        "variant": variant,
        "replica": replica,
        "pair_resname": pair[0].name if pair else None,
        "has_sensor_pair": pair is not None,
        "init_dist": round(d0, 2),
        "dG_bind_kcal": round(float(np.mean(dG)), 2) if dG else None,
        "dG_bind_sd": round(float(np.std(dG)), 2) if dG else None,
        "occupancy_4A": round(float(np.mean(np.array(dists) < 4.0)), 3) if dists else None,
        "mean_mindist": round(float(np.mean(dists)), 2) if dists else None,
        "binder_rmsd": round(float(np.mean(rmsds)), 2) if rmsds else None,
        "n_frames": len(dG),
        # per-frame His433-carboxylate min-distance series (A); kept so a switch that FLICKERS/STICKS
        # (bimodal distance distribution) is distinguishable from a stable one -- a single occupancy
        # number hides that. Empty when there is no sensor pair.
        "mindist_series": [round(float(x), 2) for x in dists],
    }


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
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    name = os.path.splitext(os.path.basename(args.pdb))[0]
    log("=== MD pH-bind (MM-GBSA):", name, "device", args.device)

    results = {"name": name, "pdb": os.path.abspath(args.pdb), "states": []}
    for rep in range(args.replicas):
        for variant in ("HIP", "HIE"):
            try:
                r = run_state(args.pdb, variant, args, args.device, rep)
                results["states"].append(r)
                log("  rep%d %s dG_bind=%s occ4=%s rmsd=%s" % (
                    rep, variant, r["dG_bind_kcal"], r["occupancy_4A"], r["binder_rmsd"]))
            except Exception as e:
                log("  rep%d %s FAILED: %r" % (rep, variant, e))
                results["states"].append({"variant": variant, "replica": rep, "error": repr(e)})

    def agg(metric, variant):
        vals = [s[metric] for s in results["states"]
                if s.get("variant") == variant and s.get(metric) is not None]
        return float(np.mean(vals)) if vals else None

    dg_hip, dg_hie = agg("dG_bind_kcal", "HIP"), agg("dG_bind_kcal", "HIE")
    occ_hip, occ_hie = agg("occupancy_4A", "HIP"), agg("occupancy_4A", "HIE")
    results["dG_bind_HIP"] = (round(dg_hip, 2) if dg_hip is not None else None)
    results["dG_bind_HIE"] = (round(dg_hie, 2) if dg_hie is not None else None)
    results["dd_bind_kcal"] = (round(dg_hip - dg_hie, 2)
                               if (dg_hip is not None and dg_hie is not None) else None)
    results["occ_HIP"] = occ_hip
    results["occ_HIE"] = occ_hie
    results["d_occupancy"] = (round(occ_hip - occ_hie, 3)
                              if (occ_hip is not None and occ_hie is not None) else None)
    # pH-selective binding: tighter (more negative dG) when protonated
    ddb = results["dd_bind_kcal"]
    results["bind_switch"] = (ddb is not None and ddb < 0.0)
    results["bind_switch_strong"] = (ddb is not None and ddb < -1.0)
    # no His433-carboxylate sensor pair -> dd_bind is a whole-interface electrostatic difference
    # (not a mechanistic switch readout); the consumer must consult this before counting the leg.
    results["has_sensor_pair"] = any(s.get("has_sensor_pair") for s in results["states"])

    outpath = os.path.join(args.out, name + ".json")
    json.dump(results, open(outpath, "w"), indent=2)
    log("WROTE", outpath, "| dG_HIP=%s dG_HIE=%s dd_bind=%s bind_switch=%s (d_occ=%s)" % (
        results["dG_bind_HIP"], results["dG_bind_HIE"], results["dd_bind_kcal"],
        results["bind_switch"], results["d_occupancy"]))
    print(outpath)


if __name__ == "__main__":
    main()
