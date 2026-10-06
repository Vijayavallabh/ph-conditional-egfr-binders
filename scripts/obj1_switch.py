#!/usr/bin/env python3
"""Objective-1 engine (CORRECTED): engineer + score a pH switch on each filtered binder.

WHY THIS REPLACES stage2_mutate.py + stage4_ph_score.py
-------------------------------------------------------
PyRosetta `-pH_mode -value_pH` + InterfaceAnalyzer on ref2015 does NOT titrate: ddG is
byte-identical at pH 2 vs 12 and His stays neutral at all pH (the pH score term is absent in
2026.39 and the protonated-His variant is not in the ResidueTypeSet). And a naive
nearest-residue -> Asp point mutation both clashes (+180 REU on an unrelaxed stub) and
CRASHES His433's pKa (desolvation, no salt bridge) -> it destroys the switch.

So we instead:
  1. SCAN binder positions near the pH sensor and install a partner, then FastRelax:
       - binder Asp/Glu near target His433 (B/100): salt bridge forms when His433 is
         protonated at low pH -> stabilises the protonated state -> tumour-ON.
       - binder His near target Glu496 (B/163, conserved acidic): the binder His protonates
         at low pH and salt-bridges the always-negative Glu -> tumour-ON.
  2. SCORE each variant by the Wyman pH-linkage computed from PROPKA pKas:
       n_H(pH) = sum_groups 1/(1+10^(pH-pKa));  dn_H = n_H(complex) - n_H(target) - n_H(binder)
       ddG_pH  = RT*ln10 * integral_{6.5}^{7.4} dn_H dpH        (RT*ln10 = 1.419 kcal/mol @310K)
     ddG_pH > 0  => complex binds more protons than the free parts => binds tighter at pH 6.5
     => tumour-ON. (ddG_pH = dG_bind(7.4) - dG_bind(6.5).)
  3. Keep the best-scoring variant per design (+ a best-acid x best-His stacking combo), and
     also report a Rosetta interface ddG (affinity, Obj3) for the kept variant.

Validated on the pipetest design: metric spans -0.198..+0.589 kcal/mol across variants; the
best switch position is NOT the closest contact, and some positions invert the switch -- so
every candidate is scored, never assumed.

Usage
-----
  # driver: shard the manifest's PDBs across worker processes
  python scripts/obj1_switch.py --manifest results/filtered_designs.csv \
      --outdir runs/obj1 --out results/obj1_switch.csv --workers 48
  # (internal) one worker over a shard file
  python scripts/obj1_switch.py --worker --shard <file> --outdir runs/obj1 --out <partial.csv>

Must run with envs/pyrosetta activated (PyRosetta 2026.39 + propka 3.5.1).
"""
from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RT_LN10 = 1.419            # kcal/mol at 310 K
HIS433_NUM = 100           # target pdbnum of His433 (crop ordinal)
GLU496_NUM = 163           # target pdbnum of Glu496
D1 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "H": "HIS"}
AA3to1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V", "MSE": "M"}

# ------------------------------------------------------------------ PROPKA linkage
def _parse_pka(pkafile):
    out, inss = [], False
    try:
        fh = open(pkafile)
    except OSError:
        return out
    for line in fh:
        if "SUMMARY OF THIS PREDICTION" in line:
            inss = True
            continue
        if inss:
            t = line.split()
            if len(t) >= 4 and t[0].isalpha() and len(t[0]) in (2, 3):
                try:
                    out.append((t[0], int(t[1]), t[2], float(t[3])))
                except ValueError:
                    pass
            elif line.strip().startswith("-") and out:
                break
    fh.close()
    return out


def _f(pka, pH):
    return 1.0 / (1.0 + 10 ** (pH - pka))


def _run_propka(pdb_basename, workdir):
    r = subprocess.run(["propka3", pdb_basename], cwd=workdir,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pka = os.path.join(workdir, pdb_basename[:-4] + ".pka")
    if r.returncode != 0:
        # propka failed (possibly after writing a PARTIAL .pka) -> discard it so _parse_pka returns []
        # and ddg_ph_from_pdb returns None (a clean "no measurement"), never a truncated Wyman linkage
        # recorded as a valid ddg_pH -- that score drives all Objective-1 ranking.
        if os.path.exists(pka):
            os.remove(pka)
    return pka


# Physical pKa clamps (Cys deliberately excluded: unreliable + unwanted as a switch). These kill
# PROPKA desolvation artifacts (buried acids at pKa ~9, buried Cys at pKa ~18) that otherwise
# dominate the linkage with non-physical, non-robust "switches".
PKA_CLAMP = {"ASP": (1.0, 7.0), "GLU": (1.5, 7.0), "HIS": (3.0, 8.5),
             "LYS": (6.0, 12.5), "ARG": (10.0, 14.0), "TYR": (8.0, 13.0)}


def _cb_coords(pdb):
    """{(resn,num,chain): (x,y,z)} by CB (fallback CA) for distance-to-sensor tests."""
    cb, ca = {}, {}
    for l in open(pdb):
        if l[:6] in ("ATOM  ", "HETATM"):
            an = l[12:16].strip()
            try:
                key = (l[17:20].strip(), int(l[22:26]), l[21])
                xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                continue
            if an == "CB":
                cb[key] = xyz
            elif an == "CA":
                ca.setdefault(key, xyz)
    for k, v in ca.items():
        cb.setdefault(k, v)
    return cb


def ddg_ph_from_pdb(pdb, binder_chain, target_chain, sensor_cut=12.0):
    """Wyman pH-linkage from PROPKA pKas, per-group matched (non-shifting groups cancel), with
    pKas clamped to physical ranges and Cys dropped. Returns the linkage RESTRICTED to groups
    near the conserved pH sensors (His433 / Glu496) as the primary score, plus the total and the
    His433-only term:
        contrib_g = RT*ln10 * integral_{6.5}^{7.4} [f(pKa_cx,g) - f(pKa_free,g)] dpH
        ddg_pH (sensor) = sum of contrib_g for g with CB within sensor_cut of His433 or Glu496
    ddg_pH > 0 => complex holds more protons at the epitope => binds tighter at pH 6.5 (tumour-ON),
    via a CONSERVED-anchor mechanism (cross-reactive + wet-lab-defensible)."""
    wd = tempfile.mkdtemp(prefix="phlink_")
    try:
        cx = os.path.join(wd, "cx.pdb")
        shutil.copyfile(pdb, cx)
        tg = os.path.join(wd, "tg.pdb")
        bd = os.path.join(wd, "bd.pdb")
        with open(tg, "w") as ot, open(bd, "w") as ob:
            for line in open(pdb):
                if line[:6] in ("ATOM  ", "HETATM"):
                    ch = line[21]
                    if ch == target_chain:
                        ot.write(line)
                    elif ch == binder_chain:
                        ob.write(line)
            ot.write("END\n")
            ob.write("END\n")

        def clamp(groups):
            out = {}
            for rn, n, c, p in groups:
                if rn in PKA_CLAMP:
                    lo, hi = PKA_CLAMP[rn]
                    out[(rn, n, c)] = min(max(p, lo), hi)
            return out

        gcx = clamp(_parse_pka(_run_propka("cx.pdb", wd)))
        gtg = clamp(_parse_pka(_run_propka("tg.pdb", wd)))
        gbd = clamp(_parse_pka(_run_propka("bd.pdb", wd)))
        if not gcx or not gtg:
            return None
        free = {}
        free.update(gtg)
        free.update(gbd)
        pairs = [(k, gcx[k], free[k]) for k in gcx if k in free]

        cb = _cb_coords(pdb)
        his_cb = cb.get(("HIS", HIS433_NUM, target_chain))
        glu_cb = cb.get(("GLU", GLU496_NUM, target_chain))

        def near_sensor(key):
            p = cb.get(key)
            if p is None:
                return False
            if his_cb and math.dist(p, his_cb) <= sensor_cut:
                return True
            if glu_cb and math.dist(p, glu_cb) <= sensor_cut:
                return True
            return False

        grid = [6.5 + 0.05 * k for k in range(19)]      # 6.5 .. 7.4

        def trapz(vals):
            return sum((vals[i] + vals[i + 1]) / 2 * (grid[i + 1] - grid[i]) for i in range(len(vals) - 1))

        def contrib(pc, pf):
            return RT_LN10 * trapz([_f(pc, p) - _f(pf, p) for p in grid])

        total = sensor = his433 = 0.0
        dn_s65 = dn_s74 = 0.0
        local = []
        his_key = ("HIS", HIS433_NUM, target_chain)
        for k, pc, pf in pairs:
            gi = contrib(pc, pf)
            total += gi
            if k == his_key:
                his433 = gi
            if near_sensor(k):
                sensor += gi
                dn_s65 += _f(pc, 6.5) - _f(pf, 6.5)
                dn_s74 += _f(pc, 7.4) - _f(pf, 7.4)
                if abs(gi) > 0.02:
                    local.append((k, gi, pc, pf))
        local.sort(key=lambda x: x[1], reverse=True)
        drivers = ";".join(f"{rn}{n}{c}:{gi:+.2f}(pKa{pf:.1f}->{pc:.1f})"
                           for (rn, n, c), gi, pc, pf in local[:3])
        his_cx = gcx.get(his_key)
        his_tg = gtg.get(his_key)
        return {
            "ddg_pH": round(sensor, 3),              # PRIMARY: conserved-sensor linkage
            "ddg_pH_total": round(total, 3),
            "ddg_pH_his433": round(his433, 3),
            "dn_6.5": round(dn_s65, 3), "dn_7.4": round(dn_s74, 3),
            "his433_cx_pka": his_cx, "his433_tg_pka": his_tg,
            "his433_dpKa": round(his_cx - his_tg, 2) if (his_cx and his_tg) else None,
            "ph_drivers": drivers,
        }
    finally:
        shutil.rmtree(wd, ignore_errors=True)


# ------------------------------------------------------------------ PyRosetta pieces
_PYR = None


def init_pyrosetta():
    global _PYR
    if _PYR is None:
        import pyrosetta
        pyrosetta.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false "
                       "-no_optH false -ex1 -ex2")
        _PYR = pyrosetta
    return _PYR


def clean_pdb(src, dst):
    """Drop placeholder atoms at the origin. BoltzGen inverse-folded outputs carry only real
    backbone coords for the binder; every side-chain atom sits at (0,0,0). Stripping them lets
    Rosetta rebuild the side chains from ideal geometry on load."""
    with open(dst, "w") as o:
        for l in open(src):
            if l[:6] in ("ATOM  ", "HETATM"):
                try:
                    x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
                except ValueError:
                    continue
                if abs(x) < 1e-3 and abs(y) < 1e-3 and abs(z) < 1e-3:
                    continue
            o.write(l)


def load_and_relax(pyr, sf, pdb):
    """Clean -> load (Rosetta builds the missing side chains) -> full FastRelax. Raw stubs load
    at ~3e4 REU (idealised side chains clash); a constrained relax settles them to ~ -800 REU
    while holding the designed backbone / binding mode (constrain_relax_to_start_coords)."""
    from pyrosetta.rosetta.protocols.relax import FastRelax
    tmp = tempfile.mktemp(suffix=".pdb")
    clean_pdb(pdb, tmp)
    try:
        pose = pyr.pose_from_file(tmp)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    fr = FastRelax(sf, 1)
    fr.constrain_relax_to_start_coords(True)
    fr.apply(pose)
    return pose


def chain_sizes(pose):
    pi = pose.pdb_info()
    from collections import Counter
    c = Counter(pi.chain(i) for i in range(1, pose.total_residue() + 1))
    return c


def find_res(pose, chain, pdbnum):
    pi = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        if pi.chain(i) == chain and pi.number(i) == pdbnum:
            return i
    return None


def near_positions(pose, binder_chain, center_xyz, cutoff, limit):
    """binder residues whose CB (or CA) is within cutoff of any center atom; sorted by distance."""
    import numpy as np
    pi = pose.pdb_info()
    cen = [np.array(x) for x in center_xyz]
    out = []
    for i in range(1, pose.total_residue() + 1):
        if pi.chain(i) != binder_chain:
            continue
        r = pose.residue(i)
        if r.name3() in ("GLY", "PRO"):      # skip; backbone-risky mutation targets
            continue
        cb = np.array(r.xyz("CB")) if r.has("CB") else np.array(r.xyz("CA"))
        d = min(float(np.linalg.norm(cb - c)) for c in cen)
        if d < cutoff:
            out.append((i, pi.number(i), r.name3(), round(d, 2)))
    out.sort(key=lambda x: x[3])
    return out[:limit]


def local_relax(pyr, pose, sf, focus_ids, radius=8.0):
    from pyrosetta.rosetta.core.select.residue_selector import (
        NeighborhoodResidueSelector, ResidueIndexSelector, NotResidueSelector)
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, PreventRepackingRLT, OperateOnResidueSubset)
    from pyrosetta.rosetta.protocols.relax import FastRelax
    sel = ResidueIndexSelector(",".join(str(x) for x in focus_ids))
    nbr = NeighborhoodResidueSelector(sel, radius, True)
    mv = nbr.apply(pose)
    mm = pyr.rosetta.core.kinematics.MoveMap()
    for i in range(1, pose.total_residue() + 1):
        mm.set_bb(i, bool(mv[i]))
        mm.set_chi(i, bool(mv[i]))
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(nbr)))
    fr = FastRelax(sf, 1)
    fr.constrain_relax_to_start_coords(True)
    fr.set_movemap(mm)
    fr.set_task_factory(tf)
    fr.apply(pose)


def interface_ddg(pyr, pose):
    try:
        ia = pyr.rosetta.protocols.analysis.InterfaceAnalyzerMover(1)
        ia.set_pack_separated(True)
        ia.set_pack_input(True)
        ia.set_compute_packstat(False)
        p = pose.clone()
        ia.apply(p)
        return round(ia.get_interface_dG(), 2)
    except Exception:
        return None


def binder_seq(pose, binder_chain):
    pi = pose.pdb_info()
    s = []
    for i in range(1, pose.total_residue() + 1):
        if pi.chain(i) == binder_chain:
            s.append(AA3to1.get(pose.residue(i).name3(), "X"))
    return "".join(s)


# ------------------------------------------------------------------ per-design driver
def make_variant(pyr, base, sf, muts, focus_extra):
    """Clone base, apply [(resid, one-letter)] mutations, local-relax, return pose."""
    from pyrosetta.toolbox import mutate_residue
    p = base.clone()
    foci = list(focus_extra)
    for resid, aa1 in muts:
        mutate_residue(p, resid, aa1, pack_radius=6.0, pack_scorefxn=sf)
        foci.append(resid)
    local_relax(pyr, p, sf, foci)
    return p


def process_design(pyr, sf, pdb, outdir, acid_cut=9.0, his_cut=9.0, max_acid=8, max_his=4):
    import numpy as np
    name = Path(pdb).stem
    try:
        pose = load_and_relax(pyr, sf, pdb)      # strip origin stubs, rebuild side chains, relax
    except Exception as e:
        return {"design": name, "error": f"load:{repr(e)[:80]}"}
    sizes = chain_sizes(pose)
    if len(sizes) < 2:
        return {"design": name, "error": "fewer than 2 chains"}
    target_chain = max(sizes, key=sizes.get)
    binder_chain = min(sizes, key=sizes.get)
    his = find_res(pose, target_chain, HIS433_NUM)
    glu = find_res(pose, target_chain, GLU496_NUM)
    if his is None:
        return {"design": name, "error": "His433 not found"}
    hr = pose.residue(his)
    his_centers = [hr.xyz(a) for a in ("ND1", "NE2") if hr.has(a)]
    glu_centers = []
    if glu is not None:
        gr = pose.residue(glu)
        glu_centers = [gr.xyz(a) for a in ("OE1", "OE2") if gr.has(a)]

    acid_pos = near_positions(pose, binder_chain, his_centers, acid_cut, max_acid)
    his_pos = near_positions(pose, binder_chain, glu_centers, his_cut, max_his) if glu_centers else []

    # baseline: relax around His433, no mutation
    variants = [("orig", [], [his])]
    for (i, pnum, aa, d) in acid_pos:
        variants.append((f"D@{pnum}", [(i, "D")], [his]))
        variants.append((f"E@{pnum}", [(i, "E")], [his]))
    for (i, pnum, aa, d) in his_pos:
        variants.append((f"H@{pnum}", [(i, "H")], [his, glu]))

    scored = []
    for tag, muts, foci in variants:
        try:
            p = make_variant(pyr, pose, sf, muts, foci)
        except Exception as e:
            continue
        tmp = os.path.join(outdir, f".{name}__{tag}.pdb".replace("@", "at"))
        p.dump_pdb(tmp)
        sc = ddg_ph_from_pdb(tmp, binder_chain, target_chain)
        os.remove(tmp)
        if sc is None:
            continue
        scored.append((tag, muts, foci, sc))

    if not scored:
        return {"design": name, "error": "no scored variants"}

    best = max(scored, key=lambda s: s[3]["ddg_pH"])

    # stacking combo: best single acid + best single His (if both exist and improve)
    best_acid = max((s for s in scored if s[0].startswith(("D@", "E@"))),
                    key=lambda s: s[3]["ddg_pH"], default=None)
    best_his = max((s for s in scored if s[0].startswith("H@")),
                   key=lambda s: s[3]["ddg_pH"], default=None)
    if best_acid and best_his:
        muts = best_acid[1] + best_his[1]
        if len({m[0] for m in muts}) == 2:          # distinct positions
            try:
                p = make_variant(pyr, pose, sf, muts, [his, glu])
                tmp = os.path.join(outdir, f".{name}__combo.pdb")
                p.dump_pdb(tmp)
                sc = ddg_ph_from_pdb(tmp, binder_chain, target_chain)
                os.remove(tmp)
                if sc is not None:
                    combo = (f"combo[{best_acid[0]}+{best_his[0]}]", muts, [his, glu], sc)
                    scored.append(combo)
                    if sc["ddg_pH"] > best[3]["ddg_pH"]:
                        best = combo
            except Exception:
                pass

    # rebuild + persist the best variant, add interface ddG + sequence
    tag, muts, foci, sc = best
    p = make_variant(pyr, pose, sf, muts, foci) if muts or tag == "orig" else pose
    outpdb = os.path.join(outdir, f"{name}__best.pdb")
    p.dump_pdb(outpdb)
    iddg = interface_ddg(pyr, p)
    seq = binder_seq(p, binder_chain)
    return {
        "design": name, "best_tag": tag,
        "ddg_pH": sc["ddg_pH"], "ddg_pH_total": sc.get("ddg_pH_total"),
        "ddg_pH_his433": sc.get("ddg_pH_his433"),
        "dn_6.5": sc["dn_6.5"], "dn_7.4": sc["dn_7.4"],
        "his433_cx_pka": sc["his433_cx_pka"], "his433_dpKa": sc["his433_dpKa"],
        "ph_drivers": sc.get("ph_drivers", ""),
        "interface_ddg": iddg, "n_variants": len(scored),
        "orig_ddg_pH": next((s[3]["ddg_pH"] for s in scored if s[0] == "orig"), None),
        "binder_chain": binder_chain, "target_chain": target_chain,
        "binder_seq": seq, "binder_len": len(seq),
        "best_pdb": outpdb, "src_pdb": pdb, "error": "",
    }


FIELDS = ["design", "best_tag", "ddg_pH", "ddg_pH_total", "ddg_pH_his433",
          "dn_6.5", "dn_7.4", "his433_cx_pka", "his433_dpKa", "ph_drivers",
          "interface_ddg", "n_variants", "orig_ddg_pH",
          "binder_chain", "target_chain", "binder_len", "binder_seq",
          "best_pdb", "src_pdb", "error"]


def run_worker(shard_file, outdir, out_csv):
    pyr = init_pyrosetta()
    sf = pyr.get_fa_scorefxn()
    os.makedirs(outdir, exist_ok=True)
    pdbs = [l.strip() for l in open(shard_file) if l.strip()]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for pdb in pdbs:
            t0 = time.time()
            try:
                row = process_design(pyr, sf, pdb, outdir)
            except Exception as e:
                row = {"design": Path(pdb).stem, "error": f"proc:{repr(e)[:100]}"}
            w.writerow({k: row.get(k, "") for k in FIELDS})
            fh.flush()
            print(f"[{Path(pdb).stem}] tag={row.get('best_tag')} ddg_pH={row.get('ddg_pH')} "
                  f"iddg={row.get('interface_ddg')} err={row.get('error')} {time.time()-t0:.0f}s",
                  flush=True)


def run_driver(pdbs, outdir, out_csv, workers):
    os.makedirs(outdir, exist_ok=True)
    shard_dir = os.path.join(outdir, "shards")
    os.makedirs(shard_dir, exist_ok=True)
    workers = max(1, min(workers, len(pdbs)))
    shards = [[] for _ in range(workers)]
    for i, p in enumerate(pdbs):
        shards[i % workers].append(p)
    procs = []
    for k, sh in enumerate(shards):
        if not sh:
            continue
        sf_file = os.path.join(shard_dir, f"shard_{k:03d}.txt")
        Path(sf_file).write_text("\n".join(sh) + "\n")
        pcsv = os.path.join(shard_dir, f"part_{k:03d}.csv")
        log = open(os.path.join(shard_dir, f"shard_{k:03d}.log"), "w")
        cmd = [sys.executable, __file__, "--worker", "--shard", sf_file,
               "--outdir", outdir, "--out", pcsv]
        procs.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), pcsv, log))
        time.sleep(0.3)       # stagger PyRosetta inits
    print(f"launched {len(procs)} workers over {len(pdbs)} designs -> {outdir}", flush=True)
    for pr, _, _ in procs:
        pr.wait()
    # merge
    rows = []
    for _, pcsv, log in procs:
        log.close()
        if os.path.exists(pcsv):
            rows += list(csv.DictReader(open(pcsv)))
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    ok = [r for r in rows if not r.get("error")]
    print(f"\nDONE {len(ok)}/{len(rows)} designs scored -> {out_csv}", flush=True)
    pos = sorted((r for r in ok if r.get("ddg_pH") not in ("", None)),
                 key=lambda r: float(r["ddg_pH"]), reverse=True)[:15]
    print("top-15 by ddg_pH (tumour-ON):")
    for r in pos:
        print(f"  {r['design'][:40]:40s} tag={r['best_tag']:16s} ddg_pH={r['ddg_pH']} "
              f"dpKa={r['his433_dpKa']} iddg={r['interface_ddg']} len={r['binder_len']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--shard")
    ap.add_argument("--manifest")
    ap.add_argument("--pdbs", nargs="*")
    ap.add_argument("--outdir", default=str(REPO / "runs" / "obj1"))
    ap.add_argument("--out", default=str(REPO / "results" / "obj1_switch.csv"))
    ap.add_argument("--workers", type=int, default=48)
    ap.add_argument("--limit", type=int, default=0, help="cap #designs (0=all)")
    a = ap.parse_args()

    if a.worker:
        run_worker(a.shard, a.outdir, a.out)
        return

    pdbs = []
    if a.manifest:
        for r in csv.DictReader(open(a.manifest)):
            if r.get("pdb"):
                pdbs.append(r["pdb"])
    if a.pdbs:
        for p in a.pdbs:
            pdbs += sorted(glob.glob(p)) if any(c in p for c in "*?[") else [p]
    pdbs = [p for p in pdbs if os.path.exists(p)]
    if a.limit:
        pdbs = pdbs[:a.limit]
    if not pdbs:
        sys.exit("no input PDBs")
    run_driver(pdbs, a.outdir, a.out, a.workers)


if __name__ == "__main__":
    main()
