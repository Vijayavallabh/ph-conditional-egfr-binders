#!/usr/bin/env python3
"""Stage 2 - build pH-switch variant complexes by point-mutating the binder (NO pH_mode).

Objective-1 engine (pH-selectivity). For each binder-target complex PDB:
  * identify the BINDER (smaller chain) and TARGET (EGFR Domain III, larger chain) by residue count;
  * locate the conserved target sensors His433 (target pdbnum 100, HIS) and Glu496 (pdbnum 163, GLU);
  * find the binder residue nearest His433's imidazole N (ND1/NE2)  -> `nh`
    and the binder residue nearest Glu496's carboxyl O (OE1/OE2)    -> `ng`;
  * emit four variant complexes (clone per variant):
        orig  - unchanged baseline
        hisD  - nh -> ASP  (binder carboxylate to salt-bridge protonated His433+ at low pH)
        hisE  - nh -> GLU  (same mechanism, longer carboxylate)
        gluH  - ng -> HIS  (binder His to partner the conserved acidic Glu496)
    using toolbox.mutate_residue (local repack, pack_radius 6.0).

Each variant is scored at pH 6.5 vs 7.4 downstream by scripts/stage4_ph_score.py; a more negative
pH_margin = ddg(6.5) - ddg(7.4) means stronger binding at tumor pH ("tumor-ON").

Run:
  python scripts/stage2_mutate.py --pdbs 'runs/.../*.pdb' --outdir out/variants --manifest out/variants.csv
  python scripts/stage2_mutate.py --pdbs @list.txt      --outdir out/variants --manifest out/variants.csv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import time
import traceback
from pathlib import Path

import pyrosetta
from pyrosetta.toolbox import mutate_residue

REPO = Path(__file__).resolve().parents[1]

HIS433_PDBNUM = 100          # target renumbered: pdbnum 100 == His433 (conserved pH sensor)
GLU496_PDBNUM = 163          # target renumbered: pdbnum 163 == Glu496 (conserved acidic sensor)
HIS_IMIDAZOLE_N = ("ND1", "NE2")    # His imidazole nitrogens (protonated at low pH)
GLU_CARBOXYL_O = ("OE1", "OE2")     # Glu carboxyl oxygens

FIELDS = ["design", "tag", "variant_pdb", "binder_chain", "target_chain",
          "mut_pose_resid", "mut_binder_pdbnum", "mut_to",
          "nearest_his433_dist_ang", "nearest_glu496_dist_ang", "error"]


def resolve_pdbs(spec: str) -> list[str]:
    """Expand --pdbs: a glob pattern, or @listfile (one path/glob per line, # comments ok)."""
    paths: list[str] = []
    if spec.startswith("@"):
        with open(spec[1:]) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                hits = sorted(glob.glob(line))
                paths += hits if hits else [line]          # keep literal so it fails recorded
    else:
        paths = sorted(glob.glob(spec))
    # de-dup, keep deterministic order
    seen, out = set(), []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def chain_residue_counts(pose) -> dict:
    pi = pose.pdb_info()
    counts: dict = {}
    for i in range(1, pose.total_residue() + 1):
        counts[pi.chain(i)] = counts.get(pi.chain(i), 0) + 1
    return counts


def find_target_sensor(pose, target_chain: str, pdbnum: int, name3: str) -> int:
    """Return the pose index of the target residue at `pdbnum` with the expected 3-letter name."""
    pi = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        if pi.chain(i) == target_chain and pi.number(i) == pdbnum:
            got = pose.residue(i).name3()
            if got != name3:
                raise ValueError(f"target {target_chain}/{pdbnum} is {got}, expected {name3}")
            return i
    raise ValueError(f"target sensor {target_chain}/{pdbnum} ({name3}) not found")


def sensor_atom_points(pose, idx: int, atom_names) -> list:
    res = pose.residue(idx)
    pts = [res.xyz(a) for a in atom_names if res.has(a)]
    if not pts:
        raise ValueError(f"residue {idx} ({res.name3()}) has none of {atom_names}")
    return pts


def nearest_binder_residue(pose, binder_idx: list[int], pts) -> tuple[int, float]:
    """Binder residue whose nearest atom (any atom) is closest to the given sensor points."""
    best_i, best_d = None, float("inf")
    for i in binder_idx:
        res = pose.residue(i)
        dmin = min((res.xyz(j) - p).norm()
                   for j in range(1, res.natoms() + 1) for p in pts)
        if dmin < best_d:
            best_d, best_i = dmin, i
    if best_i is None:
        raise ValueError("no binder residues to search")
    return best_i, best_d


def process_pdb(pdb: str, outdir: str, sf) -> list[dict]:
    stem = Path(pdb).stem
    pose = pyrosetta.pose_from_file(pdb)
    pi = pose.pdb_info()

    counts = chain_residue_counts(pose)
    if len(counts) < 2:
        raise ValueError(f"expected >=2 chains, found {sorted(counts)}")
    binder_chain = min(counts, key=lambda c: counts[c])
    target_chain = max(counts, key=lambda c: counts[c])
    if binder_chain == target_chain:
        raise ValueError(f"could not separate binder/target from counts {counts}")

    his_idx = find_target_sensor(pose, target_chain, HIS433_PDBNUM, "HIS")
    glu_idx = find_target_sensor(pose, target_chain, GLU496_PDBNUM, "GLU")
    his_pts = sensor_atom_points(pose, his_idx, HIS_IMIDAZOLE_N)
    glu_pts = sensor_atom_points(pose, glu_idx, GLU_CARBOXYL_O)

    binder_idx = [i for i in range(1, pose.total_residue() + 1) if pi.chain(i) == binder_chain]
    nh, dh = nearest_binder_residue(pose, binder_idx, his_pts)       # binder res nearest His433
    ng, dg = nearest_binder_residue(pose, binder_idx, glu_pts)       # binder res nearest Glu496
    nh_pdbnum = pi.number(nh)
    ng_pdbnum = pi.number(ng)
    dh_r, dg_r = round(dh, 3), round(dg, 3)

    os.makedirs(outdir, exist_ok=True)

    # (tag, mutate-resid or None, one-letter target or None)
    plan = [("orig", None, None),
            ("hisD", nh, "D"),
            ("hisE", nh, "E"),
            ("gluH", ng, "H")]

    rows = []
    for tag, resid, aa in plan:
        work = pose.clone()
        if resid is not None:
            mutate_residue(work, resid, aa, pack_radius=6.0, pack_scorefxn=sf)
        out_path = os.path.join(outdir, f"{stem}__{tag}.pdb")
        work.dump_pdb(out_path)
        rows.append({
            "design": stem, "tag": tag, "variant_pdb": out_path,
            "binder_chain": binder_chain, "target_chain": target_chain,
            "mut_pose_resid": "" if resid is None else resid,
            "mut_binder_pdbnum": "" if resid is None else (nh_pdbnum if resid == nh else ng_pdbnum),
            "mut_to": "" if aa is None else aa,
            "nearest_his433_dist_ang": dh_r, "nearest_glu496_dist_ang": dg_r,
            "error": "",
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Build pH-switch binder variants (stage 2, no pH_mode).")
    ap.add_argument("--pdbs", required=True, help="glob pattern OR @listfile of complex PDBs")
    ap.add_argument("--outdir", required=True, help="directory for variant PDBs")
    ap.add_argument("--manifest", required=True, help="output manifest CSV path")
    a = ap.parse_args()

    pdbs = resolve_pdbs(a.pdbs)
    print(f"[stage2] {len(pdbs)} input PDB(s) -> {a.outdir}")
    if not pdbs:
        raise SystemExit("no input PDBs matched --pdbs")

    pyrosetta.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")
    sf = pyrosetta.get_fa_scorefxn()

    os.makedirs(a.outdir, exist_ok=True)
    Path(a.manifest).parent.mkdir(parents=True, exist_ok=True)

    n_ok = n_fail = n_variants = 0
    with open(a.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for k, pdb in enumerate(pdbs, 1):
            t0 = time.time()
            try:
                rows = process_pdb(pdb, a.outdir, sf)
                for r in rows:
                    w.writerow(r)
                n_ok += 1
                n_variants += len(rows)
                print(f"  [{k}/{len(pdbs)}] {Path(pdb).name}: {len(rows)} variants "
                      f"({time.time() - t0:.2f}s)")
            except Exception as e:                              # noqa: BLE001 - robust at scale
                n_fail += 1
                w.writerow({"design": Path(pdb).stem, "tag": "FAILED", "variant_pdb": pdb,
                            "error": f"{type(e).__name__}: {e}"})
                print(f"  [{k}/{len(pdbs)}] {Path(pdb).name}: FAILED - {type(e).__name__}: {e}")
                traceback.print_exc()
            fh.flush()

    print(f"[stage2] done: {n_ok} ok, {n_fail} failed, {n_variants} variants -> {a.manifest}")


if __name__ == "__main__":
    main()
