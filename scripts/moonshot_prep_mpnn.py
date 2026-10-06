#!/usr/bin/env python3
"""WAVE 22 moonshot downstream, stage 1: RFd3 .cif.gz -> staged PDB (binder=A, target=B) +
per-design His-lock map for motif-locked soluble MPNN.

Each RFd3 design carries a binder histidine salt-bridging a CONSERVED antigen Asp (Asp347=B14 or
Asp460=B127). We must KEEP that His when ProteinMPNN redesigns the binder, so we:
  1. gunzip+parse each CIF, identify the target chain (the 191-res DIII) and the binder (the other),
  2. find the binder His whose ND1/NE2 is closest to a target Asp carboxylate (the motif bridge),
  3. write a PDB with binder=chain A (numbering preserved) + target=chain B,
  4. emit pdbs.json (list) + fixed_residues_multi.json ({pdb: "A<hispos>"}) + a manifest CSV.

QC: drop designs whose motif bridge is broken (His N .. Asp O > cutoff) -- RFd3 occasionally wanders.
"""
import argparse, csv, glob, gzip, json, os, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder

TARGET_LEN = 191            # DIII crop
BRIDGE_CUT = 4.0            # His N .. Asp O max for an intact motif
HIS = ("HIS", "HIP", "HIE", "HID")


def load_cif(path):
    # parse the decompressed CIF from an in-memory handle (code-review: /tmp is read-only / stale-serving
    # on this box, and a tmp file was also leaked on parser exceptions). Biopython's get_structure accepts
    # a file-like object, so nothing touches disk.
    import io
    with gzip.open(path, "rt") as g:
        txt = g.read()
    s = MMCIFParser(QUIET=True).get_structure("x", io.StringIO(txt))
    return s[0]


def analyze(model):
    chains = {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}
    # target = the chain closest to TARGET_LEN; binder = the other
    tgt = min(chains, key=lambda c: abs(len(chains[c]) - TARGET_LEN))
    binders = [c for c in chains if c != tgt]
    if not binders:
        return None
    bnd = max(binders, key=lambda c: len(chains[c]))
    # motif His in binder closest to a target Asp carboxylate
    asps = [r for r in chains[tgt] if r.resname == "ASP"]
    best = None
    for hr in chains[bnd]:
        if hr.resname not in HIS:
            continue
        for a in ("ND1", "NE2"):
            if a not in hr:
                continue
            for ar in asps:
                for o in ("OD1", "OD2"):
                    if o in ar:
                        d = np.linalg.norm(hr[a].coord - ar[o].coord)
                        if best is None or d < best[0]:
                            best = (d, hr.id[1], ar.id[1])
    if best is None:
        return None
    return dict(tgt=tgt, bnd=bnd, nbind=len(chains[bnd]), bridge=round(float(best[0]), 2),
                his=int(best[1]), asp=int(best[2]), chains=chains)


def write_pdb(info, out_pdb):
    """Write binder (as chain A, numbering preserved) + target (as chain B)."""
    sb = StructureBuilder()
    sb.init_structure("x"); sb.init_model(0)
    for src_chain, new_id in ((info["bnd"], "A"), (info["tgt"], "B")):
        sb.init_chain(new_id)
        for r in info["chains"][src_chain]:
            sb.init_seg(" ")
            sb.init_residue(r.resname if r.resname not in ("HIP", "HIE", "HID") else "HIS",
                            " ", r.id[1], " ")
            for at in r:
                sb.init_atom(at.name, at.coord,
                             at.bfactor if at.bfactor is not None else 0.0,
                             at.occupancy if at.occupancy is not None else 1.0, " ",
                             at.fullname, at.serial_number, at.element)
    io = PDBIO(); io.set_structure(sb.get_structure()); io.save(out_pdb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default="runs/moonshot/gen", help="dir of RFd3 gpu*/ outputs")
    ap.add_argument("--out", default="runs/moonshot/staged", help="staged PDBs + json dir")
    ap.add_argument("--bridge-cut", type=float, default=BRIDGE_CUT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cifs = sorted(glob.glob(os.path.join(a.gen, "gpu*", "*.cif.gz")))
    print(f"found {len(cifs)} RFd3 CIFs")
    pdbs, fixed, rows, kept, broken = [], {}, [], 0, 0
    for cf in cifs:
        # include the gpu subdir: RFd3 names by key+batch+model only (NOT seed), so all 8 GPU shards
        # produce identical basenames -> must disambiguate or distinct designs collide/overwrite.
        name = os.path.basename(os.path.dirname(cf)) + "_" + os.path.basename(cf).replace(".cif.gz", "")
        try:
            info = analyze(load_cif(cf))
        except Exception as e:
            print(f"  parse FAIL {name}: {repr(e)[:80]}"); continue
        if info is None:
            continue
        if info["bridge"] > a.bridge_cut:
            broken += 1; continue
        out_pdb = os.path.abspath(os.path.join(a.out, f"{name}.pdb"))
        write_pdb(info, out_pdb)
        pdbs.append(out_pdb)
        fixed[out_pdb] = f"A{info['his']}"
        rows.append(dict(name=name, binder_len=info["nbind"], his_pos=info["his"],
                         asp=info["asp"], bridge=info["bridge"]))
        kept += 1
    json.dump(pdbs, open(os.path.join(a.out, "pdbs.json"), "w"))
    json.dump(fixed, open(os.path.join(a.out, "fixed_residues_multi.json"), "w"))
    with open(os.path.join(a.out, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "his_pos", "asp", "bridge"])
        w.writeheader(); w.writerows(rows)
    print(f"kept {kept} (motif bridge <= {a.bridge_cut} A); dropped {broken} broken-bridge")
    if rows:
        import statistics as st
        print(f"  binder_len {min(r['binder_len'] for r in rows)}-{max(r['binder_len'] for r in rows)}; "
              f"bridge median {st.median(r['bridge'] for r in rows):.2f} A; "
              f"Asp347={sum(r['asp']==14 for r in rows)} Asp460={sum(r['asp']==127 for r in rows)}")
    print(f"-> {a.out}/pdbs.json, fixed_residues_multi.json, manifest.csv")


if __name__ == "__main__":
    main()
