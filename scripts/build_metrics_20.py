#!/usr/bin/env python3
"""Build the submission's per-design metrics table (results/metrics_20.csv) and the human-readable
docs/DESIGN_TABLE_20.md from the validated sources, so the 20 rows are never hand-transcribed.

Tier A  (8 locked)       <- results/metrics_8.csv        (PROPKA+pKAI+MD d_occ, +JustHISpKa where positive)
Tier B  (11 MD-selective) <- results/metrics_final_w18.csv (selected==yes; gate = explicit-protonation MD
                                                            occupancy occ_HIP high AND occ_HIE low + isostere)
Tier C  (1 control)       <- the paired E66Q knockout of egfr_phsw_22 (06__65); binds, predicted NOT selective.

Blanks are honest: a cell is empty when that leg was not run on that design (e.g. isostere/occ_HIE were
not part of the locked-8 workflow; PRODIGY/JustHISpKa were not re-run on the new 11, which instead carry
MD occupancy + isostere + MM-GBSA). Nothing is invented to fill a column.
"""
import csv, os

REPO = "."; os.chdir(REPO)
COLS = ["name", "tier", "length", "lineage", "d_occ", "occ_HIE", "ddg_pH_propka", "pkai_ddg_pH",
        "iso_causal", "mmgbsa_dd_bind", "chai_human", "chai_mouse", "boltz", "af2ig_iface_pae",
        "prodigy_dG_kcal", "free_cys", "novelty", "notes", "sequence"]


def rd(p):
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def g(d, k, default=""):
    v = d.get(k, default)
    return "" if v is None else v


# AF2 initial-guess INTERFACE PAE (pae_interaction), the quantity Adaptyv ranks on.
# Keyed on the parent binder name (the pdnov_..._des_... design that was AF2-ig folded);
# general: any row whose parent appears here gets the correct interface PAE, else blank.
iface_pae = {r["design"]: r["pae_interaction"] for r in rd("results/rfd_af2_passers_w18.csv")}


def jhp_note(jhp):
    """Tier A confirmation note: quad only where JustHISpKa agrees in DIRECTION (wyman > 0).
    egfr_phsw_08 is negative (wyman -0.74) and egfr_phsw_19 was not scored -> triple-confirmed."""
    try:
        v = float(jhp)
    except (TypeError, ValueError):
        v = None
    if v is not None and v > 0:
        return f"triple+JustHISpKa (PROPKA+pKAI+MD d_occ+JustHISpKa); JustHISpKa wyman={jhp}"
    if v is not None:
        return f"triple-confirmed (PROPKA+pKAI+MD d_occ); JustHISpKa wyman={jhp} negative, not counted"
    return "triple-confirmed (PROPKA+pKAI+MD d_occ); JustHISpKa not scored"


rows = []

# Tier A: locked 8 (resubmit order 02,04,06,08,09,12,16,19 as they appear in metrics_8)
for r in rd("results/metrics_8.csv"):
    rows.append({
        "name": r["name"], "tier": "A_locked", "length": g(r, "length"), "lineage": r["name"],
        "d_occ": g(r, "md_docc"), "occ_HIE": "", "ddg_pH_propka": g(r, "ddg_pH"),
        "pkai_ddg_pH": g(r, "pkai_ddg_pH"), "iso_causal": "", "mmgbsa_dd_bind": g(r, "dd_bind_kcal"),
        "chai_human": g(r, "chai_iptm_human"), "chai_mouse": g(r, "chai_mouse_iptm"),
        "boltz": g(r, "boltz2"), "af2ig_iface_pae": g(r, "af2_pae"),
        "prodigy_dG_kcal": g(r, "prodigy_dG_kcal"),
        "free_cys": g(r, "free_cys"), "novelty": "portal-passed",
        "notes": jhp_note(g(r, "justhispka_wyman_shift")),
        "sequence": g(r, "sequence"),
    })

# Tier B: the 11 MD-selective new designs (selected==yes)
anchor_newname = ""
for r in rd("results/metrics_final_w18.csv"):
    if r.get("selected") != "yes":
        continue
    if r["name"] == "pdnov_egfr_phsw_06__des_65_2":
        anchor_newname = r["new_name"]
    seq = g(r, "sequence")
    iso = g(r, "iso_causal") or "incon."   # finalizer writes yes/no/incon. directly (tri-state)
    note = f"MD-{g(r, 'md_tier')}; isostere-causal={iso}; parent {r['name']}"
    rows.append({
        "name": r["new_name"], "tier": "B_md_selective", "length": str(len(seq)), "lineage": g(r, "lineage"),
        "d_occ": g(r, "d_occ"), "occ_HIE": g(r, "occ_HIE"), "ddg_pH_propka": g(r, "ddg_pH"),
        "pkai_ddg_pH": g(r, "pkai_ddg_pH"), "iso_causal": iso, "mmgbsa_dd_bind": g(r, "mmgbsa_dd_bind"),
        "chai_human": g(r, "chai_iptm_human"), "chai_mouse": g(r, "chai_mouse_iptm"),
        "boltz": g(r, "boltz2_iptm"), "af2ig_iface_pae": iface_pae.get(r["name"], ""),
        "prodigy_dG_kcal": "",
        "free_cys": g(r, "free_cys"), "novelty": g(r, "novelty_level") or "portal-gated",
        "notes": note, "sequence": seq,
    })

# Tier C: the paired knockout control (read it back from the resubmit CSV so the sequence is the one shipped)
sub = {r["name"]: r for r in rd("results/adaptyv_challenge1_resubmit_FINAL.csv")}
# the control is the one resubmit row that is neither a locked-8 design nor a Tier-B pick (robust to count)
known = {r["name"] for r in rows}
ko = next((n for n in sub if n not in known), None)
if ko:
    seq = sub[ko]["sequence"]
    rows.append({
        "name": ko, "tier": "C_control", "length": str(len(seq)), "lineage": "egfr_phsw_06",
        "d_occ": "", "occ_HIE": "", "ddg_pH_propka": "", "pkai_ddg_pH": "", "iso_causal": "",
        "mmgbsa_dd_bind": "", "chai_human": "", "chai_mouse": "", "boltz": "", "af2ig_iface_pae": "",
        "prodigy_dG_kcal": "", "free_cys": "0", "novelty": "1-mut of portal-passed",
        "notes": f"paired mechanism control = {anchor_newname or 'egfr_phsw_22'} (06__65) E66Q; removes the "
                 "His433-salt-bridging Glu charge; isostere predicts switch collapses (0.99->-0.22) so it "
                 "should BIND but NOT be pH-selective (in-vitro negative control)",
        "sequence": seq,
    })

with open("results/metrics_20.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(rows)

print(f"wrote results/metrics_20.csv ({len(rows)} rows: "
      f"{sum(r['tier']=='A_locked' for r in rows)} locked + "
      f"{sum(r['tier']=='B_md_selective' for r in rows)} MD-selective + "
      f"{sum(r['tier']=='C_control' for r in rows)} control)")
assert len(rows) == 20, f"expected 20 rows, got {len(rows)}"


# ---- markdown table block for DESIGN_TABLE_20.md (compact, submission-facing) ----
def fmt(v, nd=2):
    try:
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return v if v not in (None, "") else "-"


def tbl(tier_rows, cols, hdr):
    out = ["| " + " | ".join(hdr) + " |", "|" + "|".join(["---"] * len(hdr)) + "|"]
    for r in tier_rows:
        out.append("| " + " | ".join(str(c(r)) for c in cols) + " |")
    return "\n".join(out)


A = [r for r in rows if r["tier"] == "A_locked"]
B = [r for r in rows if r["tier"] == "B_md_selective"]
C = [r for r in rows if r["tier"] == "C_control"]
colsA = [lambda r: r["name"], lambda r: r["length"], lambda r: fmt(r["d_occ"]),
         lambda r: fmt(r["ddg_pH_propka"]), lambda r: fmt(r["pkai_ddg_pH"]),
         lambda r: fmt(r["mmgbsa_dd_bind"]), lambda r: fmt(r["chai_human"]), lambda r: fmt(r["boltz"]),
         lambda r: fmt(r["prodigy_dG_kcal"], 1), lambda r: r["free_cys"]]
hdrA = ["design", "len", "MD d_occ", "PROPKA ddG_pH", "pKAI", "MM-GBSA dd_bind", "Chai iPTM", "Boltz iPTM", "PRODIGY dG", "free-Cys"]
colsB = [lambda r: r["name"], lambda r: r["lineage"], lambda r: r["length"], lambda r: fmt(r["d_occ"]),
         lambda r: fmt(r["occ_HIE"]), lambda r: fmt(r["ddg_pH_propka"]), lambda r: fmt(r["pkai_ddg_pH"]),
         lambda r: r["iso_causal"] or "incon.", lambda r: fmt(r["mmgbsa_dd_bind"]),
         lambda r: fmt(r["chai_human"]), lambda r: fmt(r["boltz"]), lambda r: fmt(r["af2ig_iface_pae"], 2)]
hdrB = ["design", "lineage", "len", "MD d_occ", "occ_HIE", "PROPKA ddG_pH", "pKAI ddG_pH", "iso-causal",
        "MM-GBSA dd_bind", "Chai iPTM", "Boltz iPTM", "AF2-ig iface PAE"]

blocks = {
    "TIER_A": tbl(A, colsA, hdrA),
    "TIER_B": tbl(sorted(B, key=lambda r: int(r["name"].split("_")[-1])), colsB, hdrB),
    "KO_LINE": (f"`{C[0]['name']}` is {anchor_newname or 'egfr_phsw_22'} (binder 06__65) with a single "
                f"**E66Q** mutation ({C[0]['length']} aa). It removes the His433-salt-bridging Glu charge; "
                f"the isostere analysis predicts the switch collapses (0.99 to -0.22) while the fold and the "
                f"rest of the interface are intact, so it should bind but lose the pH preference, a built-in "
                f"in-vitro negative control.") if C else "(no control row)",
}
DOC = f"""# Design table: twenty de novo pH-conditional EGFR domain III binders

Anthropic x Adaptyv 2026, Challenge 1 (Track 3). Single-chain, in-silico only. Objective priority:
pH-selectivity first, then human/mouse cross-reactivity, then affinity, with novelty as a gate.

The panel is 8 locked designs that were already validated and submitted (Tier A), 11 new designs chosen
on explicit-protonation molecular dynamics selectivity (Tier B), and 1 paired knockout that acts as a
built-in wet-lab mechanism control (Tier C). All 20 pass the official validator (single chain, length,
ANARCI, novelty) and carry zero free cysteines.

## How to read the primary metric

Every design puts one aspartate or glutamate so it salt-bridges His433 when His433 is protonated. A
positive static pH-switch score (PROPKA Wyman linkage, pKAI His433 shift) only confirms the switch points
the right way: tighter at tumor pH. It does not prove selectivity, because a bridge that is also present
when His433 is neutral scores positive yet binds at both pH. The one measurement that separates the two
is explicit-protonation MD: the fraction of the trajectory the designed bridge is formed with His433
protonated (occ_HIP) must be high AND the fraction with His433 neutral (occ_HIE) must be low. The gap,
d_occ = occ_HIP - occ_HIE, is the selectivity metric; low occ_HIE means little binding at blood pH.

## Tier A: 8 locked designs (already submitted)

Each locked design is confirmed by PROPKA, pKAI, and the MD d_occ; JustHISpKa agrees in direction for six
of the eight (egfr_phsw_08 is negative on JustHISpKa and egfr_phsw_19 was not scored by it), so we treat
those two as triple-confirmed. Cross-reactivity is checked by four methods and binding is measured by
PRODIGY and co-folders. d_occ here is from the earlier MD; per-state neutral-pH occupancy (occ_HIE) was not
logged for this tier, so occ_HIE is blank.

{blocks["TIER_A"]}

## Tier B: 11 new designs (selected on MD selectivity)

Partial-diffusion backbones (the novelty route), soluble ProteinMPNN (0-Cys), AF2 initial-guess interface
filter, then MD-screened on the two protonation states. All eleven clear the 0.2 selectivity floor
(d_occ >= 0.2 with occ_HIP > occ_HIE). That floor is applied to the eleven new designs only; the eight
locked designs were selected earlier under a related but not identical MD protocol, and two of them
(egfr_phsw_02 at d_occ 0.155 and egfr_phsw_16 at 0.195) sit just below it. Per-state neutral-pH occupancy
(occ_HIE) was not logged for the locked tier, so we do not compare occ_HIE across tiers. The strongest new
designs have occ_HIE near zero; egfr_phsw_28/30/31 sit at the floor with more neutral-pH leak (occ_HIE 0.45
to 0.52) and are the weaker end, kept for binding diversity. iso-causal = yes means the isostere control
(acid to its neutral amide) collapses the His433 linkage, so the switch is the designed bridge and not bulk
charge; "incon." means the static His433 linkage was already near zero so the control is uninformative there
(these still cleared the MD occupancy gate). The eleven span nine partial-diffusion lineages, capped at two
per lineage; note egfr_phsw_22/23/24 fall in one foldseek structural cluster (distinct lineages 06/06/09),
so three strong designs share that backbone family.

{blocks["TIER_B"]}

AF2-ig iface PAE is the AF2 initial-guess interface PAE (pae_interaction, from the AF2 initial-guess fold;
lower is better; the same quantity Adaptyv ranks on), sourced per design from its parent binder's AF2-ig
fold. All eleven fall in a tight 5.7 to 7.7 band, well under the 11 cutoff.

Chai and Boltz disagree on egfr_phsw_02 (Tier A, Chai iPTM 0.21 vs Boltz 0.86) and egfr_phsw_26 (Chai iPTM
0.21 vs Boltz 0.59, with AF2-ig interface PAE 6.73); both engines are shown so neither call is hidden. We
treat a confident call from any two of AF2, Boltz, and Chai as binding support and note that Chai is the
least reliable of the three on these interfaces.

## Tier C: 1 paired knockout (built-in mechanism control)

{blocks["KO_LINE"]}

If the wet lab sees egfr_phsw_22 bind more tightly at pH 6.5 than 7.4 while this single-point mutant binds
but loses that pH preference, the designed Glu66:His433 bridge is confirmed as the switch in vitro. The
mutant is still de novo (a one-residue variant of a design that passed the novelty gate).

## Provenance

- Per-design numeric table with sequences: `results/metrics_20.csv`
- Upload CSV (name, sequence, molecule_class): `results/adaptyv_challenge1_resubmit_FINAL.csv`
- Full methods and reproduce commands: `docs/METHODS_20.md`
- Engine positive control (FcRn:Fc): `scripts/obj1_fcrn_control.py` (both legs PASS)
"""
open("docs/DESIGN_TABLE_20.md", "w").write(DOC)
print("wrote docs/DESIGN_TABLE_20.md")
