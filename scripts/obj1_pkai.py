#!/usr/bin/env python3
"""Orthogonal (PROPKA-independent) pKa confirmation of the His433 pH switch via pKAI.

pKAI is a graph neural network trained to reproduce Poisson-Boltzmann (PypKa) pKas --
a methodology distinct from PROPKA's empirical scheme AND from force-field MD, so
agreement across all three = genuine triple-confirmation of the decisive objective.

For each final complex (binder=A, target=B, His433=B/100, Glu496=B/163), using
bound-geometry references (to match scripts/obj1_switch.py's PROPKA protocol):
  - pKAI on complex, isolated target (chain B), isolated binder (chain A)
  - dpKa(His433) = pKa_complex - pKa_free_target   (>0 => binder stabilizes protonated His => tumor-ON)
  - Wyman linkage ddG_pH = RT*ln10 * integral_{6.5}^{7.4} d(n_H) dpH, restricted to the conserved
    sensors His433/Glu496 + binder carboxylates within 8 A of His433; n_H = 1/(1+10^(pH-pKa)).
Writes results/pkai_switch.csv; prints sign-agreement vs the PROPKA engine.

Usage: envs/pkai/bin/python scripts/obj1_pkai.py [final20_pdbs.tsv] [--model pKAI|pKAI+] [--jobs N] [--out CSV]
"""
import csv, os, math, argparse
from pathlib import Path
os.environ.setdefault("OMP_NUM_THREADS", "4")

REPO = Path(__file__).resolve().parents[1]   # derive from script location (portable across sessions)
# repo-relative default so a later session with no tsv arg doesn't point at a vanished session scratchpad
# (code-review 2026-10-02; cf. the repo's known stale-/tmp failure mode). Pass the tsv arg for other inputs.
DEFAULT_TSV = REPO/"results"/"final20_pdbs.tsv"
TMP = REPO/"runs/tmp/pkai"; TMP.mkdir(parents=True, exist_ok=True)
# KEEP IN SYNC with scripts/obj1_switch.py (the canonical PROPKA engine): the Wyman-linkage integral below
# (RT_LN10, HIS/GLU residue nums, the pH grid and trapezoid, nH==_f) duplicates obj1_switch's _f/RT_LN10/
# HIS433_NUM/GLU496_NUM/trapz. Only the pKa SOURCE differs (pKAI vs PROPKA). A change to one (e.g. the
# 298->310 K correction) MUST be mirrored here or the two methods drift silently. (code-review 2026-10-02)
RT_LN10 = 1.419          # kcal/mol at 310 K (body temp; matches scripts/obj1_switch.py PROPKA engine)
HIS_RES, GLU_RES, HIS37_RES = 100, 163, 37  # His433=B100, Glu496=B163, His370(2nd conserved sensor)=B37
HIS_VARIANTS = {"HIP", "HIE", "HID", "HSP", "HSE", "HSD"}  # protonation-state names pKAI does not titrate

def extract_chain(pdb, chain, out):
    with open(out, "w") as fh:
        for l in open(pdb):
            if l.startswith("ATOM") and l[21] == chain:
                fh.write(l)
        fh.write("TER\nEND\n")

def normalize_his(pdb):
    """pKAI titrates only residues named exactly HIS; rename HIP/HIE/HID/HSP/HSE/HSD -> HIS so a
    protonation-state-named His433 (e.g. from an Amber/OpenMM-exported complex) is not silently dropped."""
    txt = open(pdb).read()
    if not any(v in txt for v in HIS_VARIANTS):
        return pdb
    out = str(TMP/(Path(pdb).stem + "_hisnorm.pdb"))
    with open(out, "w") as fh:
        for l in txt.splitlines(keepends=True):
            if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() in HIS_VARIANTS:
                l = l[:17] + "HIS" + l[20:]
            fh.write(l)
    return out

def atoms(pdb, chain, resnum=None, names=None):
    out = []
    for l in open(pdb):
        if not l.startswith("ATOM") or l[21] != chain:
            continue
        rn = int(l[22:26]); an = l[12:16].strip()
        if resnum is not None and rn != resnum:
            continue
        if names and an not in names:
            continue
        out.append((rn, an, (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out

def dist(a, b): return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))
def nH(pka, pH): return 1.0/(1.0+10**(pH-pka))

def run_pkai(pdb, model):
    import contextlib
    from pkai.pKAI import pKAI
    norm = normalize_his(pdb)
    with open(os.devnull, "w") as dn, contextlib.redirect_stdout(dn):
        res = pKAI(norm, model_name=model, device="cpu")
    return {(str(c), int(rn)): (resn, float(pk)) for (c, rn, resn, pk) in res}

def score_one(args):
    name, pdb, model = args          # model passed explicitly (spawn workers re-import this module)
    try:
        tgt = str(TMP/f"{name}_B.pdb"); bnd = str(TMP/f"{name}_A.pdb")
        extract_chain(pdb, "B", tgt); extract_chain(pdb, "A", bnd)
        cx, ft, fb = run_pkai(pdb, model), run_pkai(tgt, model), run_pkai(bnd, model)
        his_c, his_f = cx.get(("B", HIS_RES)), ft.get(("B", HIS_RES))
        if not his_c or not his_f:
            return {"name": name, "error": "no His433 pKa"}
        if his_c[0] != "HIS":        # a titratable non-His at B100 (crop/renumber mismatch) must not be mislabelled
            return {"name": name, "error": f"B{HIS_RES} is {his_c[0]}, not His"}
        glu_c, glu_f = cx.get(("B", GLU_RES)), ft.get(("B", GLU_RES))
        h37_c, h37_f = cx.get(("B", HIS37_RES)), ft.get(("B", HIS37_RES))  # 2nd conserved epitope His
        groups = [("B", HIS_RES, his_c[1], his_f[1])]
        if glu_c and glu_f:
            groups.append(("B", GLU_RES, glu_c[1], glu_f[1]))
        hisN = [a[2] for a in atoms(pdb, "B", HIS_RES, {"ND1", "NE2"})]
        near = sorted({rn for rn, an, co in atoms(pdb, "A", None, {"OD1", "OD2", "OE1", "OE2"})
                       if hisN and min(dist(co, h) for h in hisN) < 8.0})
        for rn in near:
            cpk, fpk = cx.get(("A", rn)), fb.get(("A", rn))
            if cpk and fpk:
                groups.append(("A", rn, cpk[1], fpk[1]))
        pHs = [6.5+0.05*i for i in range(19)]
        integ = sum(0.5*(sum(nH(c, pHs[k])-nH(f, pHs[k]) for *_, c, f in groups)
                          + sum(nH(c, pHs[k+1])-nH(f, pHs[k+1]) for *_, c, f in groups))
                    * (pHs[k+1]-pHs[k]) for k in range(len(pHs)-1))
        return {"name": name, "his433_pka_free": round(his_f[1], 2),
                "his433_pka_complex": round(his_c[1], 2),
                "his433_dpka_pkai": round(his_c[1]-his_f[1], 2),
                "_dpka_raw": his_c[1]-his_f[1],          # unrounded, for exact sign comparison
                "glu496_dpka_pkai": round(glu_c[1]-glu_f[1], 2) if (glu_c and glu_f) else "",
                "n_near_binder_acids": len(near),
                "ddg_pH_pkai": round(RT_LN10*integ, 3),
                "his433_driven_pkai": (his_c[1]-his_f[1]) > 0.3,
                "his37_dpka_pkai": round(h37_c[1]-h37_f[1], 2) if (h37_c and h37_f) else "",
                "his37_driven_pkai": bool(h37_c and h37_f and (h37_c[1]-h37_f[1]) > 0.3)}
    except Exception as e:
        return {"name": name, "error": str(e)[:140]}

def main():
    ap = argparse.ArgumentParser(description="pKAI orthogonal His433 pH-switch confirmation")
    ap.add_argument("tsv", nargs="?", default=str(DEFAULT_TSV),
                    help="TSV with name,pdb columns (default: scratchpad final20_pdbs.tsv)")
    ap.add_argument("--model", default="pKAI", choices=["pKAI", "pKAI+"])
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.tsv), delimiter="\t"))
    jobs = [(r["name"], r["pdb"], a.model) for r in rows]
    man = {r["name"]: r for r in csv.DictReader(open(REPO/"results/submission_manifest_v2.csv"))}
    rank_of = {r["name"]: (int(man[r["name"]]["rank"]) if r["name"] in man else 999) for r in rows}
    if jobs and a.jobs > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(min(a.jobs, len(jobs))) as p:
            out = p.map(score_one, jobs)
    else:
        out = [score_one(j) for j in jobs]
    res = {r["name"]: r for r in out}
    cols = ["rank", "name", "his433_pka_free", "his433_pka_complex", "his433_dpka_pkai",
            "his37_dpka_pkai", "his37_driven_pkai", "glu496_dpka_pkai", "n_near_binder_acids",
            "ddg_pH_pkai", "his433_driven_pkai", "propka_his433_dpKa", "propka_ddg_pH",
            "sign_agree_his433", "error"]
    outf = Path(a.out) if a.out else REPO/"results/pkai_switch.csv"
    outf.parent.mkdir(parents=True, exist_ok=True)
    agree = tot = drv = 0
    with open(outf, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        print(f"{'rk':>2} {'name':<13}{'Hfree':>6}{'Hcplx':>6}{'dpKa':>6}{'ddgPH':>7}{'drv':>5} | "
              f"{'P_dpKa':>7}{'P_ddg':>7}{'AGR':>5}")
        for r in sorted(rows, key=lambda x: rank_of.get(x["name"], 999)):
            nm = r["name"]; d = res.get(nm, {}); m = man.get(nm, {})
            row = {"rank": rank_of.get(nm, ""), "name": nm}
            for k in ["his433_pka_free", "his433_pka_complex", "his433_dpka_pkai", "his37_dpka_pkai",
                      "his37_driven_pkai", "glu496_dpka_pkai", "n_near_binder_acids", "ddg_pH_pkai",
                      "his433_driven_pkai", "error"]:
                row[k] = d.get(k, "")
            pdpka, pddg = m.get("his433_dpKa", ""), m.get("ddg_pH", "")
            row["propka_his433_dpKa"], row["propka_ddg_pH"] = pdpka, pddg
            sa = ""
            dpka_raw = d.get("_dpka_raw")           # exact sign (not the 2-dp-rounded value)
            if dpka_raw is not None and pdpka not in ("", None):
                try: sa = (float(dpka_raw) > 0) == (float(pdpka) > 0)
                except Exception: sa = ""
            row["sign_agree_his433"] = sa
            if sa is True: agree += 1
            if sa in (True, False): tot += 1
            if d.get("his433_driven_pkai") is True: drv += 1
            w.writerow(row)
            if d.get("error"):
                print(f"{str(rank_of.get(nm, '')):>2} {nm:<13}  ERROR: {d['error']}")
            else:
                print(f"{str(rank_of.get(nm, '')):>2} {nm:<13}{str(d.get('his433_pka_free','')):>6}{str(d.get('his433_pka_complex','')):>6}"
                      f"{str(d.get('his433_dpka_pkai','')):>6}{str(d.get('ddg_pH_pkai','')):>7}{str(d.get('his433_driven_pkai',''))[:5]:>5} | "
                      f"{str(pdpka):>7}{str(pddg):>7}{str(sa):>5}")
    print(f"\npKAI His433-driven: {drv}/{len(rows)} | His433 pKa-shift sign agreement pKAI-vs-PROPKA: {agree}/{tot}")
    print(f"wrote {outf}")

if __name__ == "__main__":
    main()
