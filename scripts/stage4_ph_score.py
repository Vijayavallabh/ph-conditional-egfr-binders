#!/usr/bin/env python3
"""Stage 4 - interface ddG of each complex at a fixed pH (pH_mode / Rosetta-pH titration).

Objective-1 scorer. Runs once per pH as a separate process (pH is baked into pyrosetta.init,
which cannot be re-initialised). For every complex PDB it computes the binder-target interface
binding energy with InterfaceAnalyzerMover across jump 1, with pH-dependent protonation enabled,
and records the protonation state (.name()) of the conserved target His433 sensor.

Pipeline usage (two runs, then merge by `pdb`):
  python scripts/stage4_ph_score.py --pH 6.5 --pdbs 'out/variants/*.pdb' --out out/ph65.csv
  python scripts/stage4_ph_score.py --pH 7.4 --pdbs 'out/variants/*.pdb' --out out/ph74.csv
  # pH_margin = ddg(6.5) - ddg(7.4);  MORE NEGATIVE = better at tumor pH ("tumor-ON").

ddg is interface dG (lower = more favorable binding).
"""
from __future__ import annotations

import argparse
import csv
import glob
import time
import traceback
from pathlib import Path

import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

REPO = Path(__file__).resolve().parents[1]

FIELDS = ["pdb", "path", "pH", "ddg", "dsasa", "his433_state", "error"]


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
                paths += hits if hits else [line]
    else:
        paths = sorted(glob.glob(spec))
    seen, out = set(), []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def target_chain_of(pose) -> str:
    """Target = larger chain (EGFR Domain III, 191 res); binder = smaller chain."""
    pi = pose.pdb_info()
    counts: dict = {}
    for i in range(1, pose.total_residue() + 1):
        counts[pi.chain(i)] = counts.get(pi.chain(i), 0) + 1
    if len(counts) < 2:
        raise ValueError(f"expected >=2 chains, found {sorted(counts)}")
    return max(counts, key=lambda c: counts[c])


def sensor_state(pose, target_chain: str, pdbnum: int) -> str:
    """Full residue-type name (incl. protonation variant) of the target sensor at `pdbnum`."""
    pi = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        if pi.chain(i) == target_chain and pi.number(i) == pdbnum:
            return pose.residue(i).name()
    return ""


def score_pdb(pdb: str, pH: float, sensor: int, sf) -> dict:
    pose = pyrosetta.pose_from_file(pdb)
    target_chain = target_chain_of(pose)

    ia = InterfaceAnalyzerMover(1)          # jump 1 = binder-target interface
    ia.set_scorefunction(sf)
    ia.set_pack_separated(True)
    ia.set_pack_input(True)
    ia.set_compute_packstat(False)
    ia.apply(pose)                          # pH-dependent repack happens here

    ddg = ia.get_interface_dG()
    dsasa = ia.get_interface_delta_sasa()
    his433_state = sensor_state(pose, target_chain, sensor)   # protonation after pH repack

    return {"pdb": Path(pdb).name, "path": pdb, "pH": pH,
            "ddg": round(ddg, 4), "dsasa": round(dsasa, 4),
            "his433_state": his433_state, "error": ""}


def main() -> None:
    ap = argparse.ArgumentParser(description="Interface ddG at a fixed pH (stage 4, pH_mode).")
    ap.add_argument("--pH", type=float, required=True, help="pH value baked into Rosetta pH_mode")
    ap.add_argument("--pdbs", required=True, help="glob pattern OR @listfile of complex PDBs")
    ap.add_argument("--out", required=True, help="output CSV path")
    ap.add_argument("--sensor", type=int, default=100,
                    help="target pdbnum of His433 sensor to report protonation for (default 100)")
    a = ap.parse_args()

    pdbs = resolve_pdbs(a.pdbs)
    print(f"[stage4] pH={a.pH} : {len(pdbs)} PDB(s) -> {a.out}")
    if not pdbs:
        raise SystemExit("no input PDBs matched --pdbs")

    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ignore_zero_occupancy false "
                   f"-pH_mode -value_pH {a.pH}")
    sf = pyrosetta.get_fa_scorefxn()

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)

    n_ok = n_fail = 0
    times = []
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for k, pdb in enumerate(pdbs, 1):
            t0 = time.time()
            try:
                row = score_pdb(pdb, a.pH, a.sensor, sf)
                w.writerow(row)
                dt = time.time() - t0
                times.append(dt)
                n_ok += 1
                print(f"  [{k}/{len(pdbs)}] {row['pdb']}: ddg={row['ddg']} "
                      f"dsasa={row['dsasa']} his433={row['his433_state']} ({dt:.2f}s)")
            except Exception as e:                              # noqa: BLE001 - robust at scale
                n_fail += 1
                w.writerow({"pdb": Path(pdb).name, "path": pdb, "pH": a.pH,
                            "ddg": "", "dsasa": "", "his433_state": "",
                            "error": f"{type(e).__name__}: {e}"})
                print(f"  [{k}/{len(pdbs)}] {Path(pdb).name}: FAILED - {type(e).__name__}: {e}")
                traceback.print_exc()
            fh.flush()

    mean_t = sum(times) / len(times) if times else 0.0
    print(f"[stage4] pH={a.pH} done: {n_ok} ok, {n_fail} failed, "
          f"mean {mean_t:.2f}s/complex -> {a.out}")


if __name__ == "__main__":
    main()
