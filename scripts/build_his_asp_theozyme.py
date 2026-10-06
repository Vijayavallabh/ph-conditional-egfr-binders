#!/usr/bin/env python3
"""MOONSHOT (WAVE 22): invert the pH switch to the BINDER side.

Our whole panel relies on the ANTIGEN His433, whose buried pKa carries a ~2 pH-unit prediction
error (the deepest uncertainty). This builds the opposite motif: a BINDER-borne histidine that
protonates at tumor pH (~6.5) and salt-bridges a CONSERVED (human=mouse) antigen aspartate on the
Domain III face. The titrating group is now one we design (controllable pKa); the antigen partner
is just a permanently-negative, conserved Asp, so the switch is cross-reactive by construction and
has the correct tumor-ON polarity (binder-His+ ... Asp- forms at 6.5, breaks at 7.4).

Output = an RFd3 theozyme PDB (chain A = DIII target 310-500; chain C res1 = the binder His placed
in H-bond geometry to the chosen Asp) + the inverted RFd3 spec JSON (fix the His ring, dock+accept
at the Asp). Mirrors the WAVE-9 His433-carboxylate motif (specs/rfd3/egfr_his433_switch.json) that
produced the triple-confirmed rfd3_g2_5, but with donor/acceptor swapped.

Conserved exposed Asp candidates on the His433 face (author / p00533 / O->His433N / relSASA):
  A323 / Asp347 / 16.3A / 0.71  (most exposed)   <- primary
  A436 / Asp460 / 12.0A / 0.19  (closest, central/EGF-blocking face)
  A434 / Asp458 / 15.4A / 0.25
"""
import json, sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select

REPO = "."
SRC = f"{REPO}/specs/rfd3/theozyme.pdb"   # has chain A = DIII target (author numbering 310-500)
OUTDIR = f"{REPO}/specs/rfd3"

HBOND = 2.8  # N(donor)...O(acceptor) target distance


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v


def build_his_for_asp(struct, asp_author, template_his_author=394):
    """Place a HIS (as chain C res 1) so its NE2 H-bonds the Asp carboxylate from the solvent side.
    Returns a dict atomname->coord for the His, plus diagnostics."""
    chA = struct[0]["A"]
    asp = chA[asp_author]
    assert asp.resname == "ASP", f"A{asp_author} is {asp.resname}, not ASP"
    od1, od2 = asp["OD1"].coord, asp["OD2"].coord
    o_mid = (od1 + od2) / 2.0
    # outward (solvent) direction: from DIII centroid toward the carboxylate
    centroid = np.mean([a.coord for r in chA for a in r], axis=0)
    u_out = unit(o_mid - centroid)
    # pick the carboxylate O whose outward exposure is greatest as the acceptor
    o_acc = od1 if np.dot(od1 - centroid, u_out) >= np.dot(od2 - centroid, u_out) else od2
    donor_pos = o_acc + HBOND * u_out  # where the His NE2 goes

    # template His sidechain (real internal geometry) from the target
    this = chA[template_his_author]
    names = ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"]
    T = {n: this[n].coord.copy() for n in names}
    # local frame: translate so NE2 at origin
    ne2 = T["NE2"].copy()
    for n in names:
        T[n] = T[n] - ne2
    ring_cen = np.mean([T[a] for a in ("CG", "ND1", "CD2", "CE1", "NE2")], axis=0)
    # rotation that maps the template's NE2->ring_centroid direction onto +u_out
    a_vec = unit(ring_cen)       # from NE2 toward ring (template frame)
    b_vec = u_out                # desired outward
    v = np.cross(a_vec, b_vec); s = np.linalg.norm(v); c = np.dot(a_vec, b_vec)
    if s < 1e-8:
        R0 = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R0 = np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))

    antigen_atoms = np.array([a.coord for r in chA for a in r if not (r.id[1] == asp_author)])

    # code-review: an additive score (cb_out + min(clash,3)) let an unbounded outward-CB term outweigh
    # the saturated clash term, so a sterically clashing rotamer could win. Gate on clash FIRST (reject
    # anything below a steric-clash floor), then maximize outward-CB among the clash-free rotamers; only
    # if NONE clear the gate do we fall back to the least-clashing one, with a loud warning.
    CLASH_MIN = 2.5  # A; min His-sidechain .. antigen distance to count as clash-free
    best_ok = None   # (cb_out, placed, clash, ne2_o): best clash-free rotamer
    best_any = None  # (clash, placed, clash, ne2_o, cb_out): least-bad fallback
    for deg in range(0, 360, 6):  # spin about u_out
        th = np.radians(deg)
        K = np.array([[0, -u_out[2], u_out[1]], [u_out[2], 0, -u_out[0]], [-u_out[1], u_out[0], 0]])
        Rspin = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
        R = Rspin @ R0
        placed = {n: R @ T[n] + donor_pos for n in names}
        # clash: min dist of His heavy sidechain atoms to antigen (excl target Asp)
        sc = np.array([placed[n] for n in ("CB", "CG", "ND1", "CD2", "CE1", "NE2")])
        d = np.linalg.norm(antigen_atoms[None, :, :] - sc[:, None, :], axis=2)
        clash = float(d.min())
        cb_out = float(np.dot(placed["CB"] - donor_pos, u_out))  # CB outward (binder side)
        ne2_o = float(np.linalg.norm(placed["NE2"] - o_acc))
        if clash >= CLASH_MIN and (best_ok is None or cb_out > best_ok[0]):
            best_ok = (cb_out, placed, clash, ne2_o)
        if best_any is None or clash > best_any[2]:
            best_any = (clash, placed, clash, ne2_o, cb_out)
    if best_ok is not None:
        cb_out, placed, clash, ne2_o = best_ok
    else:
        _, placed, clash, ne2_o, cb_out = best_any
        print("  WARN: no clash-free His rotamer (best clash %.2f A < %.1f) -> using least-clashing" % (clash, CLASH_MIN))
    return placed, dict(asp=asp_author, ne2_O=round(float(ne2_o), 2),
                        clash=round(float(clash), 2), cb_outward=round(float(cb_out), 2))


def write_theozyme(struct, placed, out_pdb):
    """Write chain A (target) + chain C (the binder His, res 1) to out_pdb via Bio.PDB (correct columns)."""
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom
    model = struct[0]
    if "C" in [c.id for c in model]:
        model.detach_child("C")   # drop the WAVE-9 Glu motif; we install our His
    elem = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C",
            "ND1": "N", "CD2": "C", "CE1": "C", "NE2": "N"}
    order = ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"]
    ch = Chain("C")
    res = Residue((" ", 1, " "), "HIS", "")
    for i, n in enumerate(order):
        res.add(Atom(n, placed[n].astype(float), 0.0, 1.0, " ", (" " + n).ljust(4), i + 1, elem[n]))
    ch.add(res)
    model.add(ch)
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_pdb)


def write_spec(name, out_pdb, asp_author, out_json, binder="25-50"):
    spec = {name: {
        "dialect": 2,
        "input": out_pdb,
        "contig": f"{binder},C1-1,{binder},/0,A310-500",
        "select_fixed_atoms": {"A310-500": "ALL", "C1": "CB,CG,ND1,CD2,CE1,NE2"},
        "select_hbond_donor": {"C1": "ND1,NE2"},                 # binder His (protonated) donates
        "select_hbond_acceptor": {f"A{asp_author}": "OD1,OD2"},  # conserved antigen Asp accepts
        "select_hotspots": {f"A{asp_author}": "OD1,OD2"},        # dock at the Asp
        "infer_ori_strategy": "hotspots",
        "is_non_loopy": True,
    }}
    json.dump(spec, open(out_json, "w"), indent=2)


def main():
    struct = PDBParser(QUIET=True).get_structure("t", SRC)
    targets = {"A323": 323, "A436": 436}  # Asp347 (most exposed), Asp460 (closest/central)
    for tag, au in targets.items():
        placed, diag = build_his_for_asp(struct, au)
        out_pdb = f"{OUTDIR}/theozyme_hisasp_{tag}.pdb"
        out_json = f"{OUTDIR}/egfr_hisasp_{tag}.json"
        write_theozyme(struct, placed, out_pdb)
        write_spec(f"egfr_hisasp_{tag}", out_pdb, au, out_json)
        print(f"{tag} (Asp{ {323:347,436:460}[au] }): NE2..O={diag['ne2_O']}A  "
              f"min-clash-to-antigen={diag['clash']}A  CB-outward={diag['cb_outward']}  "
              f"-> {out_pdb.split('/')[-1]}, {out_json.split('/')[-1]}")
        if diag["ne2_O"] > 3.3 or diag["clash"] < 2.0:
            print(f"  WARNING: geometry marginal (NE2..O {diag['ne2_O']}, clash {diag['clash']})")


if __name__ == "__main__":
    main()
