#!/usr/bin/env python3
"""WAVE 18 finalize: merge every validation leg for the MD-selective designs, rank by pH-selectivity
(MD occupancy) + binding consensus (af2ig + Chai + Boltz) + causal isostere, fold-cap for diversity,
pick the final new set, and emit the resubmit CSV (8 locked + new) + a per-design metrics table.
Robust to legs still in progress (isostere / MM-GBSA columns filled when present)."""
import csv, glob, json, os

REPO = "."; os.chdir(REPO)

def ld(p, k="name"):
    return {r[k]: r for r in csv.DictReader(open(p))} if os.path.exists(p) else {}

def ff(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d

def iso_verdict(wt, isod):
    """Tri-state so 'uninformative' is never reported as 'the control refuted the designed bridge'."""
    if wt is None or isod is None or wt <= 0:
        return "incon."   # not run, or near-zero/negative static His433 linkage = nothing to collapse
    return "yes" if isod < 0.5 * wt else "no"   # ran: bridge collapsed (causal) vs survived (not)

sel = ld("results/md_selective_w18.csv")          # hip, hie, docc, tier, fold, lin, pae, tm
cof = ld("results/cofold_consensus_w18.csv")       # chai_iptm_human, chai_mouse_iptm, boltz2_iptm
full = ld("results/metrics_w18_full.csv")          # sequence, ddg_pH, pkai_*, xreact_delta_pae, novelty_level, free_cys
iso = ld("results/isostere_sel11.csv")             # wt_ddg_his433, iso_ddg_his433
mmg = {}
# jsons land either flat (<name>.json) or nested (<name>.json/<name>__best.json, when --out was a dir)
for jf in glob.glob("runs/mmgbsa_w18/**/*__best.json", recursive=True) + glob.glob("runs/mmgbsa_w18/*.json"):
    if not os.path.isfile(jf):
        continue
    try:
        d = json.load(open(jf)); n = os.path.basename(jf)[:-5].replace("__best", "")
        if d.get("dd_bind_kcal") is not None:
            mmg[n] = ff(d.get("dd_bind_kcal"))
    except Exception:
        pass

rows = []
for name, s in sel.items():
    f, c, i = full.get(name, {}), cof.get(name, {}), iso.get(name, {})
    chai_h, boltz = ff(c.get("chai_iptm_human")), ff(c.get("boltz2_iptm"))
    wt, isod = ff(i.get("wt_ddg_his433")), ff(i.get("iso_ddg_his433"))
    rows.append(dict(
        name=name, sequence=f.get("sequence", ""), fold=s["fold"], lineage=s["lin"],
        d_occ=ff(s["docc"]), occ_HIP=ff(s["hip"]), occ_HIE=ff(s["hie"]),
        ddg_pH=ff(f.get("ddg_pH")), pkai_ddg_pH=ff(f.get("pkai_ddg_pH")),
        xreact_delta_pae=ff(f.get("xreact_delta_pae")), af2_pae=ff(s["pae"]),
        chai_iptm_human=chai_h, chai_mouse_iptm=ff(c.get("chai_mouse_iptm")), boltz2_iptm=boltz,
        struct_tm=ff(s["tm"]), novelty_level=f.get("novelty_level", ""), free_cys=f.get("free_cys", "0"),
        iso_causal=iso_verdict(wt, isod),
        mmgbsa_dd_bind=mmg.get(name),
        bind_consensus=round(min(chai_h or 0.0, boltz or 0.0), 3),
        md_tier=s["tier"],
    ))

# explicit MD-selectivity floor, fail-loud: d_occ >= 0.2 AND occ_HIP > occ_HIE (the locked-8 calibration,
# whose accepted floor reaches d_occ 0.155 / occ_HIE 0.69). md_selective_w18 is already the gated 11/32,
# so this re-asserts the gate in-script (an ungated input can't silently ship a constitutive bridge) and
# logs any drop rather than capping silently.
FLOOR_DOCC = 0.2
def passes_floor(r):
    return (r["d_occ"] or 0) >= FLOOR_DOCC and (r["occ_HIP"] or 0) > (r["occ_HIE"] or 0)
dropped_floor = [r for r in rows if not passes_floor(r)]
if dropped_floor:
    print("floor-gate dropped %d: %s" % (len(dropped_floor),
          ", ".join("%s(d_occ=%s,occ_HIE=%s)" % (r["name"], r["d_occ"], r["occ_HIE"]) for r in dropped_floor)))
rows = [r for r in rows if passes_floor(r)]

# rank: selectivity (d_occ) primary; binding consensus secondary; penalize residual neutral-pH leak (occ_HIE)
def score(r):
    return round((r["d_occ"] or 0) + 0.5 * (r["bind_consensus"] or 0) - 0.4 * max(0.0, (r["occ_HIE"] or 0) - 0.2), 3)
for r in rows:
    r["score"] = score(r)
rows.sort(key=lambda r: -r["score"])

# lineage-cap <=2: cap true cousins (same partial-diffusion parent lineage), the diversity axis the
# user flagged ("distinct designs, not cousins of the same few"). All 11 MD-selective pass (<=2/lineage),
# so none of the gate-clearing designs is dropped just for sharing a coarse foldseek structural cluster.
picked, linct = [], {}
for r in rows:
    if linct.get(r["lineage"], 0) >= 2:
        r["selected"] = ""
        continue
    picked.append(r); linct[r["lineage"]] = linct.get(r["lineage"], 0) + 1
    r["selected"] = "yes"
for i, r in enumerate(picked, start=21):
    r["new_name"] = f"egfr_phsw_{i}"

cols = ["new_name", "name", "fold", "lineage", "score", "md_tier", "d_occ", "occ_HIP", "occ_HIE",
        "ddg_pH", "pkai_ddg_pH", "iso_causal", "mmgbsa_dd_bind", "xreact_delta_pae", "af2_pae",
        "chai_iptm_human", "chai_mouse_iptm", "boltz2_iptm", "bind_consensus", "struct_tm",
        "novelty_level", "free_cys", "selected", "sequence"]
with open("results/metrics_final_w18.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
    for r in rows:
        r.setdefault("new_name", "")
        w.writerow(r)

# paired isostere-knockout mechanism control: egfr_phsw_22's binder is 06__65; its lone His433
# salt-bridging Glu (Glu66, 2.66 A) -> Gln removes the charge. Verified single E66Q extraction from
# runs/wave16k_isostere/..._iso.pdb (wt==design, exactly one diff E66Q). Isostere ablation predicts the
# switch collapses (ddg_His433 0.99 -> -0.22, full sign-flip) while the fold/binding is preserved: so it
# should BIND but NOT be pH-selective -- an in-vitro negative control that isolates the designed bridge.
KO_ANCHOR = "pdnov_egfr_phsw_06__des_65_2"
KO_POS = 66   # 1-indexed binder position of the His433-salt-bridging Glu (isostere run: partner GLU66A, 2.66 A)
anchor = next((r for r in picked if r["name"] == KO_ANCHOR), None)
ko_seq = None
if anchor is not None:
    s = anchor["sequence"]
    assert len(s) >= KO_POS and s[KO_POS - 1] == "E", \
        "KO anchor pos%d is %r, expected E (Glu); the isostere partner moved -- re-derive before shipping" \
        % (KO_POS, s[KO_POS - 1:KO_POS])
    ko_seq = s[:KO_POS - 1] + "Q" + s[KO_POS:]   # E66Q charge-removed isostere, derived from the shipped anchor
    assert sum(a != b for a, b in zip(s, ko_seq)) == 1 and len(ko_seq) == len(s), \
        "KO must differ from its anchor by exactly one residue"

# resubmit CSV: 8 locked + picked new (+ knockout control when its anchor is shipped)
locked = list(csv.DictReader(open("results/designs_8.csv")))
up = [{"name": r["name"], "sequence": r["sequence"], "molecule_class": "single_chain"} for r in locked]
for r in picked:
    up.append({"name": r["new_name"], "sequence": r["sequence"], "molecule_class": "single_chain"})
if anchor is not None:
    ko_name = f"egfr_phsw_{20 + len(picked) + 1}"   # next index after the picked block
    up.append({"name": ko_name, "sequence": ko_seq, "molecule_class": "single_chain"})
    print(f"+ knockout control {ko_name} = {anchor['new_name']} (06__65) E66Q [mechanism negative control]")
else:
    print(f"WARNING: knockout anchor {KO_ANCHOR} not in picked set -> control omitted (would be orphaned)")
with open("results/adaptyv_challenge1_resubmit_FINAL.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["name", "sequence", "molecule_class"]); w.writeheader(); w.writerows(up)

print(f"picked {len(picked)} new (fold-capped) across {len({r['fold'] for r in picked})} folds, "
      f"{len({r['lineage'] for r in picked})} lineages | resubmit total = {len(up)}")
print(f"{'new':12}{'orig':30}{'docc':>6}{'HIE':>6}{'chai_h':>7}{'boltz':>7}{'iso':>5}{'mmg':>7}  fold")
for r in picked:
    print(f"{r['new_name']:12}{r['name']:30}{(r['d_occ'] or 0):6.2f}{(r['occ_HIE'] or 0):6.2f}"
          f"{(r['chai_iptm_human'] or 0):7.2f}{(r['boltz2_iptm'] or 0):7.2f}{'Y' if r['iso_causal'] else '-':>5}"
          f"{(r['mmgbsa_dd_bind'] if r['mmgbsa_dd_bind'] is not None else 0):7.2f}  {r['fold'].split('__')[1] if '__' in r['fold'] else r['fold']}")
