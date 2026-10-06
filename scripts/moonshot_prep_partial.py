#!/usr/bin/env python3
"""WAVE 22c downstream, stage 1: RFd1 PARTIAL-DIFFUSION backbones -> staged PDB (binder=A, target=B) +
per-design His-lock for the inverted switch.

RFd1 partial diffusion is BACKBONE-ONLY (binder = poly-GLY N,CA,C,O; target B keeps only backbone too,
so Asp460=B127 has no carboxylate). Partial diffusion PRESERVES the residue register (contig Lb-Lb, same
numbering) and, at low-moderate noise, the binding pose -- so the parent's switch position (moon01 His41)
stays the switch spot in each child. CA-proximity alone is the WRONG criterion (a residue 2 A closer in CA
may point its sidechain away; the parent His41 bridges at 2.93 A from CA-CA 10.4 A, while the nearest-CA
residue 37 does not). So we:
  1. lock the switch His at the PARENT position (--his-pos, default 41), register-preserved across diffusion,
  2. SET that residue to HIS in the staged PDB (backbone-only His is fine; soluble-MPNN fixes it by the
     PDB resname identity, independent of sidechain atoms) so His-locked MPNN keeps a His there,
  3. write binder=chain A + target=chain B and emit pdbs.json + fixed_residues_multi.json + manifest.

QC: drop children where that position's CA drifted beyond --max-ca of Asp460-CA (can no longer bridge);
the af2ig fold + 2-state MD are the real switch gates downstream (this just seeds a plausible His).
"""
import argparse, csv, glob, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder

ASP_ORDINAL = 127          # Asp460 in the B1-191 crop
TARGET_LEN = 191


def analyze(model, max_ca, his_pos):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders:
        return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    # require the hardcoded ordinal to actually be an Asp (guards a wrong parent/crop silently using
    # whatever residue sits at position ASP_ORDINAL)
    asp = [r for r in chains[tgt] if r.id[1] == ASP_ORDINAL and r.resname == "ASP"]
    pos = [r for r in chains[bnd] if r.id[1] == his_pos]
    if not asp or "CA" not in asp[0] or not pos or "CA" not in pos[0]:
        return None
    d = float(np.linalg.norm(pos[0]["CA"].coord - asp[0]["CA"].coord))
    if d > max_ca:
        return None
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), his=his_pos,
                ca_dist=round(d, 2), chains=chains)


def write_pdb(info, out_pdb):
    sb = StructureBuilder()
    sb.init_structure("x"); sb.init_model(0)
    for src_chain, new_id in ((info["bnd"], "A"), (info["tgt"], "B")):
        sb.init_chain(new_id)
        for r in info["chains"][src_chain]:
            # set the chosen binder residue to HIS (the switch His); keep everything else as-is
            resname = "HIS" if (new_id == "A" and r.id[1] == info["his"]) else r.resname
            sb.init_seg(" "); sb.init_residue(resname, " ", r.id[1], " ")
            for at in r:
                sb.init_atom(at.name, at.coord,
                             at.bfactor if at.bfactor is not None else 0.0,
                             at.occupancy if at.occupancy is not None else 1.0,
                             " ", at.fullname, at.serial_number, at.element)
    io = PDBIO(); io.set_structure(sb.get_structure()); io.save(out_pdb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/rfd/pd_inv02", help="dir of RFd1 gpu*/des_*.pdb")
    ap.add_argument("--out", default="runs/moonshot/pd_inv02/staged")
    ap.add_argument("--his-pos", type=int, default=41, help="parent switch-His position to lock (moon01=41)")
    ap.add_argument("--max-ca", type=float, default=12.0, help="max switch-His-CA .. Asp460-CA (A); parent 10.4")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    # FINAL designs only: a.gen/<shard>/des_<int>.pdb (one level). Exclude the per-step trajectory
    # PDBs (RFd1 writes des_N_*_traj.pdb into a traj/ subdir -- a recursive glob would ingest those).
    pdbs_in = sorted(p for p in glob.glob(os.path.join(a.gen, "*", "des_*.pdb"))
                     if "traj" not in os.path.basename(p))
    print(f"found {len(pdbs_in)} RFd1 partial backbones (final only)")
    pdbs, fixed, rows, kept, dropped = [], {}, [], 0, 0
    for cf in pdbs_in:
        # include the gpu/T subdir so des_0 from different shards don't collide
        sub = os.path.relpath(os.path.dirname(cf), a.gen).replace(os.sep, "_")
        name = f"{sub}_{os.path.basename(cf).replace('.pdb','')}"
        try:
            info = analyze(PDBParser(QUIET=True).get_structure("x", cf)[0], a.max_ca, a.his_pos)
        except Exception as e:
            print(f"  parse FAIL {name}: {repr(e)[:80]}"); continue
        if info is None:
            dropped += 1; continue
        out_pdb = os.path.abspath(os.path.join(a.out, f"{name}.pdb"))
        write_pdb(info, out_pdb)
        pdbs.append(out_pdb)
        fixed[out_pdb] = f"A{info['his']}"
        rows.append(dict(name=name, binder_len=info["nbind"], his_pos=info["his"], ca_dist=info["ca_dist"]))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(fixed, open(os.path.join(a.out, "fixed_residues_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "his_pos", "ca_dist"])
        w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (binder-CA..Asp460-CA <= {a.max_ca} A); dropped {dropped}")
    if kept == 0 and pdbs_in:
        import sys
        sys.exit(f"FATAL: kept 0 of {len(pdbs_in)} backbones -- check ASP_ORDINAL={ASP_ORDINAL} is an Asp "
                 f"in the target crop and --his-pos {a.his_pos} is within {a.max_ca} A of it.")
    if rows:
        import statistics as st
        print(f"  CA-dist median {st.median(r['ca_dist'] for r in rows):.2f} A; "
              f"binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}")
    print(f"-> {a.out}/pdbs.json, fixed_residues_multi.json, manifest.csv")


if __name__ == "__main__":
    main()
