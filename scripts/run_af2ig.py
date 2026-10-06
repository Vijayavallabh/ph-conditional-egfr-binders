#!/usr/bin/env python3
"""Fold-validate a directory of binder-target complex PDBs with AF2 initial-guess across GPUs.

Shards the PDBs across the given GPUs, runs one scripts/fold_af2ig.sh per GPU (CUDA_VISIBLE_DEVICES
pinned), then merges every shard's out.sc into one scorefile. Self-consistency / binder filter:
  binder_aligned_rmsd < 2.0  (the sequence folds to the designed backbone)
  pae_interaction    < 10     (confident interface)
  plddt_binder       > 80     (confident binder fold)
pH-AGNOSTIC: this validates the binder, not the pH switch.

Usage:
  python scripts/run_af2ig.py --pdbdir runs/obj1_foldin --out runs/af2ig --gpus 0-2 --mode exact
"""
from __future__ import annotations
import argparse, glob, os, shutil, subprocess, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SH = REPO / "scripts" / "fold_af2ig.sh"


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
    ap.add_argument("--pdbdir", required=True)
    ap.add_argument("--out", default=str(REPO / "runs" / "af2ig"))
    ap.add_argument("--gpus", default="0-2")  # user ≤3-GPU cap (2026-10-02, default 0,1,2); pass --gpus to override
    ap.add_argument("--mode", default="exact")
    a = ap.parse_args()

    pdbs = sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb")))
    if not pdbs:
        sys.exit(f"no PDBs in {a.pdbdir}")
    gpus = parse_gpus(a.gpus)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    shard_root = out / "shardin"
    if shard_root.exists():
        shutil.rmtree(shard_root)
    shard_root.mkdir(parents=True)

    # Sort by binder (chain A) length and block-shard so each GPU sees few distinct lengths
    # (af2ig compiles per length; round-robin would force ~every length to recompile on every GPU).
    def chainA_len(pdb):
        n = 0
        try:
            for l in open(pdb):
                if l[:4] == "ATOM" and l[21] == "A" and l[12:16].strip() == "CA":
                    n += 1
        except OSError:
            pass
        return n
    pdbs = sorted(pdbs, key=chainA_len)
    k = len(gpus)
    per = (len(pdbs) + k - 1) // k
    shards = [pdbs[i * per:(i + 1) * per] for i in range(k)]

    procs = []
    for k, (gpu, sh) in enumerate(zip(gpus, shards)):
        if not sh:
            continue
        sin = shard_root / f"gpu{gpu}"; sin.mkdir(parents=True, exist_ok=True)
        for p in sh:
            os.symlink(os.path.abspath(p), sin / os.path.basename(p))
        sout = out / f"gpu{gpu}"
        logf = open(out / f"gpu{gpu}.log", "w")
        # ABSOLUTE paths: fold_af2ig.sh cd's into the kit dir, which breaks relative paths.
        cmd = ["bash", str(SH), str(gpu), os.path.abspath(sin), os.path.abspath(sout), a.mode]
        procs.append((subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT), gpu, sout, logf))
        time.sleep(1.0)
    print(f"launched {len(procs)} af2ig workers over {len(pdbs)} PDBs | gpus={gpus} mode={a.mode}", flush=True)
    failed = []
    for pr, gpu, sout, logf in procs:
        rc = pr.wait(); logf.close()
        if rc != 0 or not (sout / "out.sc").exists():
            failed.append((gpu, rc))
    if failed:
        # A length-block-sharded GPU crash (OOM on the longest binders) else drops that whole length
        # range silently -- surface it so the caller knows the scorefile is INCOMPLETE, not clean.
        print(f"WARNING: {len(failed)} af2ig shard(s) FAILED or produced no out.sc: "
              f"{', '.join(f'gpu{g}(rc={rc})' for g, rc in failed)} -- merged scores are INCOMPLETE",
              flush=True)

    # merge scorefiles
    merged = out / "out.sc"
    header = None
    lines = []
    for _, gpu, sout, _ in procs:
        sc = sout / "out.sc"
        if not sc.exists():
            continue
        for ln in open(sc):
            if ln.startswith("SCORE:") and "description" in ln:
                header = ln
            elif ln.startswith("SCORE:"):
                lines.append(ln)
    with open(merged, "w") as fh:
        if header:
            fh.write(header)
        fh.writelines(lines)
    print(f"\nmerged {len(lines)} scored -> {merged}", flush=True)
    # quick pass/fail summary
    if header:
        cols = header.split()[1:]
        def gi(name):
            return cols.index(name) if name in cols else None
        ri, pi, pb = gi("binder_aligned_rmsd"), gi("pae_interaction"), gi("plddt_binder")
        npass = 0
        for ln in lines:
            v = ln.split()[1:]
            try:
                if float(v[ri]) < 2.0 and float(v[pi]) < 10.0 and float(v[pb]) > 80.0:
                    npass += 1
            except (ValueError, IndexError, TypeError):
                pass
        print(f"self-consistent (rmsd<2 & pae_int<10 & plddt_binder>80): {npass}/{len(lines)}")


if __name__ == "__main__":
    main()
