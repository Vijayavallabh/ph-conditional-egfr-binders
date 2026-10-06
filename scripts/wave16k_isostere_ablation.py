#!/usr/bin/env python3
"""WAVE 16k / deep-research M4: isostere ablation -- is the pH switch the designed Asp/Glu <-> His433
salt bridge, or bulk polyanion electrostatics?

For each frozen design: relax the complex, identify the binder Asp/Glu whose carboxylate salt-bridges
His433, then mutate ONLY that residue to its charge-removed ISOSTERE (Asp->Asn, Glu->Gln; size + H-bond
geometry preserved, charge removed), local-relax, and recompute the His433 Wyman linkage. If the switch
is the designed salt bridge, removing the partner's charge should DROP His433's pKa and collapse the
His433-specific linkage (ddg_pH_his433 -> ~0). If the switch survives the isostere, it is bulk
electrostatics, not the bridge -- a mechanism red flag.

Reuses obj1_switch.py verbatim (load_and_relax / make_variant / ddg_ph_from_pdb / HIS433_NUM=100).
Runs on the SAME complex PDBs M1 used (runs/ens_propka/*.json 'pdb'). CPU (PyRosetta + propka3);
leaves all GPUs to the pmx-NES FEP (M2'). Validation only; frozen set unchanged.

Run: PATH=envs/pyrosetta/bin:$PATH envs/pyrosetta/bin/python scripts/wave16k_isostere_ablation.py
"""
import sys, os, csv, glob, json
# propka3 (called by obj1_switch with cwd=tempdir) must be found by ABSOLUTE path
os.environ["PATH"] = "./envs/pyrosetta/bin:" + os.environ.get("PATH", "")
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import obj1_switch as O

ISO = {"ASP": "N", "GLU": "Q"}   # charge-removed isostere (one-letter for mutate_residue)

def partner_in_pose(pose, binder_chain, his_resid):
    """binder ASP/GLU with sidechain-O nearest His433 ND1/NE2 -> (resid, dist, resn)."""
    import numpy as np
    pi = pose.pdb_info()
    hr = pose.residue(his_resid)
    hn = [np.array(hr.xyz(a)) for a in ("ND1", "NE2") if hr.has(a)]
    best = (None, 1e9, None)
    for i in range(1, pose.total_residue() + 1):
        if pi.chain(i) != binder_chain:
            continue
        r = pose.residue(i)
        if r.name3() not in ISO:
            continue
        os_ = [np.array(r.xyz(a)) for a in ("OD1", "OD2", "OE1", "OE2") if r.has(a)]
        for o in os_:
            for n in hn:
                dd = float(np.linalg.norm(o - n))
                if dd < best[1]:
                    best = (i, dd, r.name3())
    return best

def main():
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    glob_arg = sys.argv[sys.argv.index("--glob") + 1] if "--glob" in sys.argv else os.path.join(O.REPO, "runs/ens_propka/egfr_phsw_*.json")
    out_arg = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(O.REPO, "results/wave16k_isostere_ablation.csv")
    pyr = O.init_pyrosetta(); sf = pyr.get_fa_scorefxn()
    outdir = os.path.join(O.REPO, "runs/wave16k_isostere"); os.makedirs(outdir, exist_ok=True)
    rows = []
    for jf in sorted(glob.glob(glob_arg)):
        d = json.load(open(jf)); name = d["name"]; pdb = d["pdb"]
        if only and name != only:
            continue
        rec = {"name": name, "err": ""}
        if not os.path.exists(pdb):
            rec["err"] = "no pdb"; rows.append(rec); print(name, "no pdb", flush=True); continue
        try:
            pose = O.load_and_relax(pyr, sf, pdb)
            sizes = O.chain_sizes(pose); tc = max(sizes, key=sizes.get); bc = min(sizes, key=sizes.get)
            his = O.find_res(pose, tc, O.HIS433_NUM)
            if his is None:
                rec["err"] = "no His433"; rows.append(rec); print(name, "no His433", flush=True); continue
            pr, pd, prn = partner_in_pose(pose, bc, his)
            if pr is None:
                rec["err"] = "no binder Asp/Glu"; rows.append(rec); print(name, "no acidic partner", flush=True); continue
            wt_pdb = os.path.join(outdir, f"{name}__wt.pdb"); pose.dump_pdb(wt_pdb)
            wt = O.ddg_ph_from_pdb(wt_pdb, bc, tc)
            p2 = O.make_variant(pyr, pose, sf, [(pr, ISO[prn])], [his])
            iso_pdb = os.path.join(outdir, f"{name}__iso.pdb"); p2.dump_pdb(iso_pdb)
            iso = O.ddg_ph_from_pdb(iso_pdb, bc, tc)
            rec.update(
                partner=f"{prn}{pose.pdb_info().number(pr)}{bc}", partner_dist=round(pd, 2),
                wt_his_cx_pka=(wt or {}).get("his433_cx_pka"), iso_his_cx_pka=(iso or {}).get("his433_cx_pka"),
                wt_ddg_his433=(wt or {}).get("ddg_pH_his433"), iso_ddg_his433=(iso or {}).get("ddg_pH_his433"),
                wt_ddg_sensor=(wt or {}).get("ddg_pH"), iso_ddg_sensor=(iso or {}).get("ddg_pH"))
            print(f"{name}: {rec['partner']} | His433 pKa {rec['wt_his_cx_pka']}->{rec['iso_his_cx_pka']} "
                  f"| ddg_His433 {rec['wt_ddg_his433']}->{rec['iso_ddg_his433']} "
                  f"| ddg_sensor {rec['wt_ddg_sensor']}->{rec['iso_ddg_sensor']}", flush=True)
        except Exception as e:
            rec["err"] = "proc:" + repr(e)[:80]; print(name, "ERR", rec["err"], flush=True)
        rows.append(rec)
    cols = ["name", "partner", "partner_dist", "wt_his_cx_pka", "iso_his_cx_pka",
            "wt_ddg_his433", "iso_ddg_his433", "wt_ddg_sensor", "iso_ddg_sensor", "err"]
    out = out_arg
    if only:
        return
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({k: r.get(k, "") for k in cols})
    # summary: did the isostere collapse the His433 linkage?
    ok = [r for r in rows if not r["err"] and r.get("wt_ddg_his433") is not None and r.get("iso_ddg_his433") is not None]
    print(f"\nDONE {len(ok)}/{len(rows)} scored -> {out}")
    print(f"{'design':15s}{'partner':12s}{'His pKa WT->iso':>18}{'ddg_His433 WT->iso':>22}{'collapsed?':>12}")
    for r in ok:
        dpka = f"{r['wt_his_cx_pka']}->{r['iso_his_cx_pka']}"
        dddg = f"{r['wt_ddg_his433']}->{r['iso_ddg_his433']}"
        collapsed = (r["iso_ddg_his433"] is not None and r["wt_ddg_his433"] not in (None, 0)
                     and r["iso_ddg_his433"] < 0.5 * r["wt_ddg_his433"])
        print(f"{r['name']:15s}{r.get('partner',''):12s}{dpka:>18}{dddg:>22}{('YES' if collapsed else 'no'):>12}")

if __name__ == "__main__":
    main()
