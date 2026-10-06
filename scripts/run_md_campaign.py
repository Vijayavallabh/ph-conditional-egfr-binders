#!/usr/bin/env python3
"""Shard the MD pH-switch validation (scripts/md_pH_switch.py) across GPUs.

Takes a file listing complex PDB paths (one per line), splits them across the given GPUs
(OpenCL device index == physical GPU), runs one md_pH_switch.py per design, then summarizes
the switch metrics (d_occupancy>0 and d_mindist>0 = salt bridge forms when His433 protonated).

Usage:
  python scripts/run_md_campaign.py --pdb-list designs.txt --out runs/md_submission \
      --gpus 0-2 --equil-ps 200 --prod-ps 2000 --replicas 3
"""
from __future__ import annotations
import argparse, glob, json, os, subprocess, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MD = REPO / "scripts" / "md_pH_switch.py"
PY = REPO / "envs" / "md" / "bin" / "python"


def parse_gpus(s):
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-"); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-list", required=True)
    ap.add_argument("--out", default=str(REPO / "runs" / "md_submission"))
    ap.add_argument("--gpus", default="0-2")  # user ≤3-GPU cap (2026-10-02, default 0,1,2); pass --gpus to override
    ap.add_argument("--equil-ps", default="200")
    ap.add_argument("--prod-ps", default="2000")
    ap.add_argument("--interval-ps", default="10")
    ap.add_argument("--replicas", default="3")
    ap.add_argument("--restraint-k", default="2.0")
    ap.add_argument("--sensor-cut", default="6.0")
    ap.add_argument("--binder-his", action="store_true",
                    help="inverted switch: titrate the BINDER His, bridge to a TARGET Asp/Glu")
    ap.add_argument("--switch-his-pos", type=int, default=None,
                    help="LOCKED binder switch-His id, forwarded to md_pH_switch.py (avoids measuring "
                         "an incidental His); use when every design in the list locks the His at one position")
    ap.add_argument("--charge-scale", type=float, default=1.0,
                    help="ECC charge scaling, forwarded to md_pH_switch.py (e.g. 0.8)")
    a = ap.parse_args()

    pdbs = [l.strip() for l in open(a.pdb_list) if l.strip() and not l.startswith("#")]
    pdbs = [p for p in pdbs if os.path.exists(p)]
    if not pdbs:
        sys.exit("no existing PDBs in list")
    gpus = parse_gpus(a.gpus)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    done = {Path(j).stem for j in glob.glob(str(out / "*.json"))}
    before = len(pdbs)
    pdbs = [p for p in pdbs if Path(p).stem not in done]
    if before != len(pdbs):
        print(f"skipping {before - len(pdbs)} already-done designs", flush=True)
    if not pdbs:
        print("all designs already done"); return

    # round-robin assign designs to GPUs
    buckets = {g: [] for g in gpus}
    for i, p in enumerate(pdbs):
        buckets[gpus[i % len(gpus)]].append(p)

    def launch_bucket(gpu, paths):
        """sequential shell running md_pH_switch.py for each design on one GPU."""
        inner = []
        for p in paths:
            inner.append(
                "%s %s --pdb %s --out %s --device %d --equil-ps %s --prod-ps %s "
                "--interval-ps %s --replicas %s --restraint-k %s --sensor-cut %s%s%s%s" % (
                    PY, MD, p, out, gpu, a.equil_ps, a.prod_ps, a.interval_ps,
                    a.replicas, a.restraint_k, a.sensor_cut,
                    " --binder-his" if a.binder_his else "",
                    (" --switch-his-pos %d" % a.switch_his_pos) if a.switch_his_pos is not None else "",
                    (" --charge-scale %g" % a.charge_scale) if a.charge_scale != 1.0 else ""))
        script = " ; ".join(inner)
        logf = open(out / f"md_gpu{gpu}.log", "w")
        return subprocess.Popen(["bash", "-c", script], stdout=logf, stderr=subprocess.STDOUT), logf

    procs = []
    for g in gpus:
        if buckets[g]:
            procs.append(launch_bucket(g, buckets[g]))
            time.sleep(0.5)
    print(f"launched MD on {len(procs)} gpus over {len(pdbs)} designs", flush=True)
    for pr, logf in procs:
        pr.wait(); logf.close()

    # summarize
    rows = []
    for j in sorted(glob.glob(str(out / "*.json"))):
        try:
            d = json.load(open(j))
            rows.append(d)
        except Exception:
            pass
    nsw = sum(1 for d in rows if d.get("md_switch"))
    print("\n%-52s %7s %7s %6s" % ("design", "d_occ", "d_dist", "switch"))
    for d in sorted(rows, key=lambda d: -(d.get("d_occupancy") or -9)):
        print("%-52s %7s %7s %6s" % (d["name"][:52], d.get("d_occupancy"), d.get("d_mindist"), d.get("md_switch")))
    print(f"\nMD switch (d_occupancy>0.1): {nsw}/{len(rows)}")


if __name__ == "__main__":
    main()
