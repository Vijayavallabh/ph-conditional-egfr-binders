#!/usr/bin/env python3
"""FAIR pH-switch head-to-head: a GFN2-xTB single-point interaction energy that treats cation-pi and salt
bridge on the SAME electronic footing -- the honest alternative to fixed-charge MD d_occ, which under-counts
cation-pi (missing polarization/charge-transfer) and over-counts salt bridges (ionic over-stabilization).

For a contact (binder His + target partner Phe/Asp), we build the folded-pose cluster in the His HIP
(protonated, tumor pH 6.5) and HIE (neutral, blood 7.4) states -- same heavy-atom geometry, protonation via
pdbfixer -- cap each sidechain to a small molecule (His->4-methylimidazole(ium), Phe->toluene, Asp->acetate),
and compute E_int = E(complex) - E(His) - E(partner) with xtb --gfn 2 --alpb water. The protonation-driven
switch strength is dE = E_int(HIP) - E_int(HIE) (more negative = protonation strengthens the contact). xtb's
self-consistent charges capture cation-pi polarization AND the aqueous screening of salt bridges, so
dE(cation-pi) vs dE(salt-bridge) is a fair comparison that the fixed-charge d_occ cannot give.
"""
import argparse, json, os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, "./scripts")

REPO = "."
XTB = f"{REPO}/envs/xtb/bin/xtb"
H2KCAL = 627.509474
BACKBONE = {"N", "C", "O", "CA", "H", "HA", "HA2", "HA3", "H1", "H2", "H3", "OXT"}
SC = {"HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"], "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
      "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"], "ASP": ["CB", "CG", "OD1", "OD2"],
      "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
      "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # indole (WAVE27 cation-pi)
      # NEUTRAL acceptors (charge 0) -- the WAVE 25 premise: do these raise the His pKa (stack) like a salt bridge?
      "SER": ["CB", "OG"], "THR": ["CB", "OG1", "CG2"], "ASN": ["CB", "CG", "OD1", "ND2"],
      "GLN": ["CB", "CG", "CD", "OE1", "NE2"]}
CHRG = {("HIS", "HIP"): 1, ("HIS", "HIE"): 0, ("PHE", None): 0, ("TYR", None): 0, ("ASP", None): -1,
        ("GLU", None): -1, ("SER", None): 0, ("THR", None): 0, ("ASN", None): 0, ("GLN", None): 0,
        ("TRP", None): 0}
ELEM = lambda name: "".join(c for c in name if c.isalpha())[0]


def protonate(pdb, his_chain_idx, his_resid, variant):
    from pdbfixer import PDBFixer
    from openmm.app import Modeller, ForceField
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    residues = list(modeller.topology.residues())
    variants = [None] * len(residues)
    for i, r in enumerate(residues):
        if r.chain.index == his_chain_idx and r.name.startswith("HI") and str(r.id) == str(his_resid):
            variants[i] = variant
    modeller.addHydrogens(ff, pH=7.0, variants=variants)
    from openmm import unit
    pos = np.array(modeller.positions.value_in_unit(unit.angstrom))
    return modeller.topology, pos


def extract_cluster(topology, pos, chain_idx, resid, kind):
    """return (elements, coords) of the capped sidechain small-molecule for residue (chain_idx, resid)."""
    res = None
    for r in topology.residues():
        if r.chain.index == chain_idx and str(r.id) == str(resid):
            res = r; break
    if res is None:
        raise RuntimeError(f"residue {chain_idx}/{resid} not found")
    atoms = {a.name: a.index for a in res.atoms()}
    keep = []
    ca = pos[atoms["CA"]] if "CA" in atoms else None
    cb = pos[atoms["CB"]] if "CB" in atoms else None
    for a in res.atoms():
        if a.name in BACKBONE:
            continue
        keep.append((ELEM(a.name), pos[a.index]))
    # cap: add an H where CA was, along CB->CA, so CB becomes a methyl carbon
    if ca is not None and cb is not None:
        d = ca - cb; d = d / (np.linalg.norm(d) + 1e-9)
        keep.append(("H", cb + 1.09 * d))
    return keep


def write_xyz(path, cluster):
    with open(path, "w") as f:
        f.write(f"{len(cluster)}\n\n")
        for el, xyz in cluster:
            f.write(f"{el} {xyz[0]:.4f} {xyz[1]:.4f} {xyz[2]:.4f}\n")


OPT = False  # module-level; set by main() from --opt
SOLVENT = "water"  # module-level; set by main() from --solvent. "vacuum" = no implicit solvent (eps=1, buried proxy)


def xtb_energy(xyz, chrg, wd):
    mode = "--opt" if OPT else "--sp"
    cmd = [XTB, os.path.basename(xyz), "--gfn", "2", "--chrg", str(chrg), mode]
    if SOLVENT != "vacuum":
        cmd[4:4] = ["--alpb", SOLVENT]  # insert after --gfn 2
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=wd)
    e = None
    for ln in r.stdout.splitlines():
        if "TOTAL ENERGY" in ln:
            e = float(ln.split()[-3])   # Eh (last one printed = final/optimized)
    if e is None:
        raise RuntimeError(f"xtb failed (chrg {chrg}):\n{r.stdout[-400:]}\n{r.stderr[-200:]}")
    return e


def e_int(pdb, his_chain_idx, his_resid, partner_chain_idx, partner_resid, partner_kind, variant, wd):
    topo, pos = protonate(pdb, his_chain_idx, his_resid, variant)
    his = extract_cluster(topo, pos, his_chain_idx, his_resid, "HIS")
    par = extract_cluster(topo, pos, partner_chain_idx, partner_resid, partner_kind)
    qh = CHRG[("HIS", variant)]; qp = CHRG[(partner_kind, None)]
    write_xyz(os.path.join(wd, "cplx.xyz"), his + par)
    write_xyz(os.path.join(wd, "his.xyz"), his)
    write_xyz(os.path.join(wd, "par.xyz"), par)
    Ec = xtb_energy("cplx.xyz", qh + qp, wd); Eh = xtb_energy("his.xyz", qh, wd); Ep = xtb_energy("par.xyz", qp, wd)
    return (Ec - Eh - Ep) * H2KCAL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--his", type=int, required=True); ap.add_argument("--partner", type=int, required=True)
    ap.add_argument("--partner-kind", default="PHE",
                    choices=["PHE", "TYR", "TRP", "ASP", "GLU", "SER", "THR", "ASN", "GLN"])
    ap.add_argument("--binder-chain", default="A"); ap.add_argument("--target-chain", default="B")
    ap.add_argument("--name", default=None)
    ap.add_argument("--opt", action="store_true", help="geometry-optimize each fragment/complex (fairer E_int)")
    ap.add_argument("--solvent", default="water",
                    help="ALPB solvent (water, chcl3, ...) or 'vacuum' (eps=1) as a DESOLVATED/BURIED proxy")
    a = ap.parse_args()
    global OPT, SOLVENT; OPT = a.opt; SOLVENT = a.solvent
    from pdbfixer import PDBFixer
    chain_ids = [c.id for c in PDBFixer(filename=a.pdb).topology.chains()]
    bi, ti = chain_ids.index(a.binder_chain), chain_ids.index(a.target_chain)
    with tempfile.TemporaryDirectory(dir=f"{REPO}/runs") as wd:
        e_hip = e_int(a.pdb, bi, a.his, ti, a.partner, a.partner_kind, "HIP", wd)
        e_hie = e_int(a.pdb, bi, a.his, ti, a.partner, a.partner_kind, "HIE", wd)
    res = dict(name=a.name or os.path.basename(a.pdb), his=a.his, partner=a.partner, kind=a.partner_kind,
               solvent=SOLVENT, opt=OPT,
               Eint_HIP=round(e_hip, 2), Eint_HIE=round(e_hie, 2), dE_switch=round(e_hip - e_hie, 2))
    print(json.dumps(res))


if __name__ == "__main__":
    main()
