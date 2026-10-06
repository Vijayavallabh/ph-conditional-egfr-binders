#!/usr/bin/env python3
"""Rebuild AF2-ig inputs from MPNN-designed RFd complexes with a FULL-ATOM target.

AF2 initial-guess needs the TARGET with sidechains to template it (backbone-only target ->
plddt_target ~40, no interface, 0% yield). RFd keeps the target rigid but translates the system,
so we Kabsch-superpose the real full-atom 6ARU DIII onto each design's target backbone and emit:
  chain A = MPNN binder (designed seq, backbone; renumbered 1..Lb)
  chain B = full-atom target in the design's frame (renumbered 1..191, His433 = B/100)

Usage:
  build_fullatom_af2in.py --mpnn 'runs/rfd_mpnn_w1/gpu*/backbones/*.pdb' \
      --target data/targets/6ARU_domainIII.pdb --out runs/rfd_af2in_w1fa
"""
import argparse, glob, os
import numpy as np


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cP, cQ


def load_target(path, chain="A"):
    """full-atom target: ordered list of (resnum, line); CA dict."""
    lines, ca = [], {}
    for ln in open(path):
        if ln.startswith("ATOM") and ln[21] == chain:
            rn = int(ln[22:26]); an = ln[12:16].strip()
            lines.append((rn, ln))
            if an == "CA":
                ca[rn] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return lines, ca


def complex_chains(path):
    """return {chain: {'ca': {resnum:xyz}, 'lines':[...], 'n':int}} for ATOM records."""
    ch = {}
    for ln in open(path):
        if ln.startswith("ATOM"):
            c = ln[21]; an = ln[12:16].strip(); rn = int(ln[22:26])
            g = ch.setdefault(c, {"ca": {}, "res": set(), "lines": []})
            g["lines"].append(ln); g["res"].add(rn)
            if an == "CA":
                g["ca"][rn] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    for c in ch:
        ch[c]["n"] = len(ch[c]["res"])
    return ch


def xyz_of(ln):
    return np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])


def put_xyz(ln, xyz):
    return "%s%8.3f%8.3f%8.3f%s" % (ln[:30], xyz[0], xyz[1], xyz[2], ln[54:].rstrip("\n"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mpnn", required=True, help="glob of MPNN backbone complexes")
    ap.add_argument("--target", default="data/targets/6ARU_domainIII.pdb")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-rmsd", type=float, default=2.0, help="skip if target superpose RMSD exceeds this")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    tgt_lines, tgt_ca = load_target(a.target, a.target_chain)
    tgt_resnums = sorted(tgt_ca)
    tgt_all_xyz = np.array([xyz_of(ln) for _, ln in tgt_lines])

    n_ok = n_skip = 0
    worst = 0.0
    for comp in sorted(glob.glob(a.mpnn)):
        ch = complex_chains(comp)
        if len(ch) < 2:
            n_skip += 1; continue
        order = sorted(ch, key=lambda c: -ch[c]["n"])
        target_ch, binder_ch = order[0], order[1]
        comp_ca = ch[target_ch]["ca"]
        comp_resnums = sorted(comp_ca)
        # match by ORDER (same 191-residue domain), not resnum: PD targets are numbered 1..191,
        # de novo RFd targets 310..500 -> resnum matching fails for PD. Order is robust to both.
        m = min(len(tgt_resnums), len(comp_resnums))
        if m < 50:
            n_skip += 1; continue
        P = np.array([tgt_ca[tgt_resnums[i]] for i in range(m)])    # real target (by order)
        Q = np.array([comp_ca[comp_resnums[i]] for i in range(m)])  # design-frame target (by order)
        R, cP, cQ = kabsch(P, Q)
        rmsd = float(np.sqrt(np.mean(np.sum(((P - cP) @ R.T + cQ - Q) ** 2, axis=1))))
        worst = max(worst, rmsd)
        if rmsd > a.max_rmsd:
            n_skip += 1; continue
        # transform all target atoms into the design frame
        tgt_moved = (tgt_all_xyz - cP) @ R.T + cQ

        out_lines = []
        serial = 1
        # binder -> chain A, renumber 1..Lb
        resmap = {}; nr = 0
        for ln in ch[binder_ch]["lines"]:
            old = ln[22:27]
            if old not in resmap:
                nr += 1; resmap[old] = nr
            out_lines.append("%s%5d %sA%4d %s\n" % (ln[:6], serial, ln[12:21], resmap[old], ln[27:].rstrip("\n")))
            serial += 1
        out_lines.append("TER\n")
        # full-atom target -> chain B, renumber 1..191 (His433 author 409 -> B/100)
        resmap = {}; nr = 0
        for (rn, ln), xyz in zip(tgt_lines, tgt_moved):
            if rn not in resmap:
                nr += 1; resmap[rn] = nr
            nl = put_xyz(ln, xyz)
            out_lines.append("%s%5d %sB%4d %s\n" % (nl[:6], serial, nl[12:21], resmap[rn], nl[27:].rstrip("\n")))
            serial += 1
        out_lines.append("TER\nEND\n")
        open(os.path.join(a.out, os.path.basename(comp)), "w").write("".join(out_lines))
        n_ok += 1
    print(f"wrote {n_ok} af2ig inputs (full-atom target) -> {a.out} | skipped {n_skip} | worst superpose RMSD {worst:.2f} A")


if __name__ == "__main__":
    main()
