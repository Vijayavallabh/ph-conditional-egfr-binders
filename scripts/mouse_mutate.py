#!/usr/bin/env python3
"""Build FULL-ATOM mouse-DIII complexes for a direct cross-reactivity cofold.

af2ig needs the target's side chains (a backbone-only target gives a false "no interface" for
every binder). So we mutate chain B (human EGFR DIII) to the mouse sequence at every divergent
position and repack chain B + the interface with PyRosetta, keeping the binder (chain A) and the
binding pose. Human/mouse DIII share the backbone and all epitope anchors (His433 incl.), so a
truly cross-reactive binder should fold-validate against the mouse complex with pae like the human.

Usage: python scripts/mouse_mutate.py --indir runs/final20 --outdir runs/final20_mouse --workers N
"""
from __future__ import annotations
import argparse, glob, json, os, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def divergent_map():
    m = json.load(open(REPO / "data" / "targets" / "mouse_domainIII_map.json"))
    out = {}
    for p in m["positions"]:
        if not p.get("same", True):
            out[int(p["ordinal"])] = p["mouse"]
    return out


def run_worker(shard, outdir, control=False):
    import pyrosetta
    from pyrosetta.toolbox import mutate_residue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ChainSelector, NeighborhoodResidueSelector, NotResidueSelector, OrResidueSelector)
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, PreventRepackingRLT, OperateOnResidueSubset)
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    pyrosetta.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false -no_optH false -ex1 -ex2")
    sf = pyrosetta.get_fa_scorefxn()
    div = {} if control else divergent_map()   # control: repack chain B only (human), for a fair pae baseline
    suffix = "_humanrepack" if control else "_mouse"
    os.makedirs(outdir, exist_ok=True)
    for pdb in shard:
        name = Path(pdb).stem
        try:
            pose = pyrosetta.pose_from_file(pdb)
        except Exception as e:
            print(f"{name}: load fail {repr(e)[:60]}", flush=True); continue
        pi = pose.pdb_info()
        # normalise chain-B pdb numbering to crop ordinal: some complexes number B continuously after
        # the binder (e.g. egfr_phsw_19 is 62..252, not 1..191), which would otherwise mutate the WRONG
        # positions. offset = 0 for the standard 1..191 case (no change to those complexes).
        bnums = [pi.number(i) for i in range(1, pose.total_residue() + 1) if pi.chain(i) == "B"]
        boff = (min(bnums) - 1) if bnums else 0
        mutated = 0
        for i in range(1, pose.total_residue() + 1):
            ordv = pi.number(i) - boff
            if pi.chain(i) == "B" and ordv in div:
                aa = div[ordv]
                if pose.residue(i).name1() != aa:
                    try:
                        mutate_residue(pose, i, aa, pack_radius=0.0, pack_scorefxn=sf)
                        mutated += 1
                    except Exception:
                        pass
        # repack chain B + interface neighborhood
        chB = ChainSelector("B")
        nbr = NeighborhoodResidueSelector(chB, 6.0, True)
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(nbr)))
        task = tf.create_task_and_apply_taskoperations(pose)
        PackRotamersMover(sf, task).apply(pose)
        out = os.path.join(outdir, f"{name}{suffix}.pdb")
        pose.dump_pdb(out)
        print(f"{name}: mutated {mutated} -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--shard")
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--control", action="store_true", help="human repack baseline (no mutations)")
    a = ap.parse_args()
    if a.worker:
        shard = [l.strip() for l in open(a.shard) if l.strip()]
        run_worker(shard, a.outdir, control=a.control)
        return
    pdbs = sorted(glob.glob(os.path.join(a.indir, "*.pdb")))
    if not pdbs:
        sys.exit("no input PDBs")
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    shard_dir = outdir / "shards"; shard_dir.mkdir(exist_ok=True)
    k = max(1, min(a.workers, len(pdbs)))
    shards = [pdbs[i::k] for i in range(k)]
    procs = []
    for j, sh in enumerate(shards):
        if not sh:
            continue
        sf = shard_dir / f"s{j}.txt"; sf.write_text("\n".join(sh) + "\n")
        log = open(shard_dir / f"s{j}.log", "w")
        cmd = [sys.executable, __file__, "--worker", "--shard", str(sf), "--indir", a.indir, "--outdir", str(outdir)]
        if a.control:
            cmd.append("--control")
        procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT))
    for p in procs:
        p.wait()
    print(f"mouse complexes: {len(glob.glob(os.path.join(a.outdir, '*_mouse.pdb')))} -> {a.outdir}")


if __name__ == "__main__":
    main()
