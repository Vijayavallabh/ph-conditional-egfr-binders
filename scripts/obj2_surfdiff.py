#!/usr/bin/env python
"""
WAVE 16b - SurfDiff binder-free cross-species (human<->mouse) surface check for objective 2.

An ADDITIVE 4th cross-reactivity signal, computed from the TARGETS ALONE (no binder, no co-folder) with
SurfDiff (Jedlitzke, ..., Sormanni 2025, bioRxiv 10.64898/2025.12.23.696170; PMID 41509378; GPL/CC-BY-NC).
SurfDiff `compare` scores each human-DIII residue's Residue Uniqueness Score (RUS) vs the mouse DIII:
HIGH RUS = species-divergent surface neighbourhood (5 A), LOW RUS = conserved. For each design we aggregate
RUS over the DIII residues it contacts (its epitope patch): LOW mean-RUS = conserved epitope = cross-reactive.

POSITIVE CONTROL (gates trust): RUS must be elevated at the 24 known human<->mouse divergent ordinals (esp.
the cetuximab species-escape cluster 134/158/159/164 that makes cetuximab human-specific) vs the conserved
residues (one-sided Mann-Whitney). If the control fails, discard (do NOT use the signal).

ADDITIVE check; does NOT override the validated WAVE-15 3-signal obj-2 consensus or re-rank.
NOTE (honest): SurfDiff's RUS is a NEIGHBOURHOOD/patch signal, so a conserved anchor (His433) can read as
high-RUS when its 5 A surroundings are divergent -- a more pessimistic, patch-level view than our
residue-level contact-severity leg; they are complementary, not interchangeable.

Contact parsing mirrors scripts/xreact_mechanism.py (ATOM-only, altLoc filter, chain-A/B any-atom
min-distance, ordinal = chain-B resnum - (min chain-B resnum - 1)); kept fail-loud rather than silent.

Inputs : SurfDiff compare output (human_query_info.csv) + the 20 complexes (contacts) + the divergence map.
Output : results/wave16b_surfdiff_obj2.csv
Run    : third_party/surfdiff/.venv/bin/python scripts/obj2_surfdiff.py <human_query_info.csv>
         (needs numpy + scipy; the SurfDiff venv has both. No SurfDiff package import, so `uv run python`
          also works if the core env provides scipy.)
"""
from __future__ import annotations
import os, sys, csv, glob, json
import numpy as np

REPO = os.environ.get("ADAPTYV_REPO") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QINFO = sys.argv[1] if len(sys.argv) > 1 else None
IN_GLOB = f"{REPO}/runs/md_bind_gplus/in/egfr_phsw_*.pdb"
MAP = f"{REPO}/data/targets/mouse_domainIII_map.json"
OUT = f"{REPO}/results/wave16b_surfdiff_obj2.csv"
CET = [134, 158, 159, 164]          # cetuximab species-escape cluster (human-specific)
CROP_LEN = 191                      # DIII crop length (ordinals 1..191); asserted against the map
CONTACT_CUT, CORE_CUT = 5.0, 4.0


def load_rus(qinfo):
    with open(qinfo, encoding="utf-8-sig") as fh:      # tolerate a UTF-8 BOM
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"[FATAL] {qinfo} has no data rows.")
    rus_cols = [c for c in rows[0] if c.endswith("_RUS")]
    if len(rus_cols) != 1:
        sys.exit(f"[FATAL] expected exactly one *_RUS column in {qinfo}, found {rus_cols}.")
    rus_c = rus_cols[0]
    if "residue_number" not in rows[0]:
        sys.exit(f"[FATAL] no 'residue_number' column in {qinfo} (cols: {list(rows[0])}).")
    out = {}
    for r in rows:
        rn, rv = r.get("residue_number"), r.get(rus_c)
        if not rn or rv is None or str(rv).strip() == "":
            continue
        try:
            v = float(rv)
        except ValueError:
            continue
        if v != v:         # NaN
            continue
        out[int(float(rn))] = v
    return out, rus_c


def divergent_ordinals():
    m = json.load(open(MAP))
    pos = m["positions"]
    if len(pos) != CROP_LEN:
        sys.exit(f"[FATAL] map has {len(pos)} positions, expected CROP_LEN={CROP_LEN}.")
    return sorted(int(p["ordinal"]) for p in pos if not p.get("same", True))


def contacts_of(pdb):
    """Ordinals of chain-B residues with any atom within CONTACT_CUT / CORE_CUT of chain A (binder).
    ATOM-only, altLoc in {' ','A'}; fail-loud on empty chains or non-contiguous chain-B numbering.
    KEEP IN SYNC with scripts/xreact_mechanism.py's contact/ordinal parser (same boff = chain-B
    resnum - (min-1) mapping; note the footgun that some designs number chain B
    62..252): a change to the mapping must be mirrored there, or the severity leg and this RUS leg
    silently map the same contact to different DIII residues."""
    A, B = [], {}
    for ln in open(pdb):
        if ln[:6] != "ATOM  ":
            continue
        if ln[16] not in (" ", "A"):           # altLoc filter (match canonical)
            continue
        ch = ln[21]
        if ch not in ("A", "B"):
            continue
        try:
            resnum = int(ln[22:26])            # parsed inside the guard (match xreact_mechanism)
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        if ch == "A":
            A.append(xyz)
        else:
            B.setdefault(resnum, []).append(xyz)
    if not A or not B:
        sys.exit(f"[FATAL] {pdb}: empty chain A ({len(A)} atoms) or B ({len(B)} residues).")
    bnums = sorted(B)
    if bnums[-1] - bnums[0] + 1 != len(bnums):
        sys.exit(f"[FATAL] {pdb}: chain B numbering not contiguous "
                 f"({bnums[0]}..{bnums[-1]}, {len(bnums)} residues) -> ordinal map unreliable.")
    if len(bnums) != CROP_LEN:
        sys.exit(f"[FATAL] {pdb}: chain B has {len(bnums)} residues, expected {CROP_LEN}.")
    A = np.asarray(A)
    boff = bnums[0] - 1
    cont, core = set(), set()
    for resid, atoms in B.items():
        d = np.sqrt(((np.asarray(atoms)[:, None, :] - A[None, :, :]) ** 2).sum(-1)).min()
        if d < CONTACT_CUT:
            cont.add(resid - boff)
        if d < CORE_CUT:
            core.add(resid - boff)
    return cont, core


def main():
    if not QINFO or not os.path.exists(QINFO):
        sys.exit("[FATAL] pass the SurfDiff human_query_info.csv path")
    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        sys.exit("[FATAL] scipy required (use third_party/surfdiff/.venv/bin/python, or a core env with scipy).")

    rus, rus_c = load_rus(QINFO)
    div = divergent_ordinals()
    if set(CET) - set(div):
        sys.exit(f"[FATAL] cetuximab cluster {CET} not a subset of the map's divergent set.")
    # fail-loud frame alignment: the RUS CSV must be in the same 1..191 ordinal frame as the map
    if 100 not in rus:
        sys.exit("[FATAL] His433 (ordinal 100) absent from RUS CSV -> wrong numbering frame; aborting.")
    present_div = [o for o in div if o in rus]
    if len(present_div) < 20:
        sys.exit(f"[FATAL] only {len(present_div)}/24 divergent ordinals present in RUS CSV -> frame mismatch.")
    missing_cet = [o for o in CET if o not in rus]
    if missing_cet:
        sys.exit(f"[FATAL] cetuximab cluster ordinals {missing_cet} absent from RUS CSV -> the positive "
                 "control would average an empty set (NaN) and silently fail; aborting.")

    cons = [o for o in range(1, CROP_LEN + 1) if o not in div]
    rus_div = [rus[o] for o in div if o in rus]
    rus_cons = [rus[o] for o in cons if o in rus]

    # ---- POSITIVE CONTROL (one-sided, tie-corrected via scipy) ----
    U, p = mannwhitneyu(rus_div, rus_cons, alternative="greater")
    md, mc = float(np.mean(rus_div)), float(np.mean(rus_cons))
    cet_mean = float(np.mean([rus[o] for o in CET if o in rus]))
    print("==== WAVE16b SurfDiff obj-2 positive control (human vs mouse DIII) ====")
    print(f"RUS divergent({len(rus_div)}) mean={md:.3f} median={np.median(rus_div):.3f} | "
          f"conserved({len(rus_cons)}) mean={mc:.3f} median={np.median(rus_cons):.3f}")
    print(f"cetuximab-escape cluster {CET} mean RUS={cet_mean:.3f}")
    print(f"Mann-Whitney (one-sided, tie-corrected) divergent>conserved: U={U:.0f} p={p:.2e}")
    ctrl_pass = (md > mc) and (p < 0.05) and (cet_mean > mc)
    print(f"CONTROL: {'PASS (RUS elevated at divergent/cetuximab surface)' if ctrl_pass else 'FAIL -> discard signal'}")
    if not ctrl_pass:
        print("Control failed -> SurfDiff signal discarded (not used).")
        return 2

    # ---- per-design aggregation over epitope contacts ----
    pdbs = sorted(glob.glob(IN_GLOB))
    if not pdbs:                       # fail loud BEFORE opening OUT in 'w' (else we'd clobber it to 0 bytes)
        sys.exit(f"[FATAL] no complexes matched {IN_GLOB} -> nothing to score "
                 f"(refusing to clobber {OUT}).")
    rows = []
    for pdb in pdbs:
        name = os.path.splitext(os.path.basename(pdb))[0]
        cont, core = contacts_of(pdb)
        cont_r = [rus[o] for o in cont if o in rus]
        core_r = [rus[o] for o in core if o in rus]
        if not cont_r:
            sys.exit(f"[FATAL] {name}: no contact ordinal mapped into the RUS CSV -> frame/contact bug.")
        div_contacts = sorted(o for o in cont if o in div)
        rows.append({
            "name": name,
            "n_contacts": len(cont), "n_core": len(core),
            "mean_rus_contacts": round(float(np.mean(cont_r)), 3),
            "max_rus_contacts": round(float(np.max(cont_r)), 3),
            "mean_rus_core": round(float(np.mean(core_r)), 3) if core_r else float("nan"),
            "n_divergent_contacts": len(div_contacts),
            "divergent_contacts": ";".join(map(str, div_contacts)),
        })
        print(f"{name}: contacts {len(cont)} (core {len(core)}) mean_RUS={rows[-1]['mean_rus_contacts']} "
              f"max_RUS={rows[-1]['max_rus_contacts']} div_contacts={div_contacts}")
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    mr = np.array([r["mean_rus_contacts"] for r in rows], float)   # all finite (fail-loud above)
    print(f"\nper-design mean_rus_contacts: min={mr.min():.3f} median={np.median(mr):.3f} max={mr.max():.3f} "
          f"(lower = more conserved epitope = cross-reactive)")
    print(f"vs DIII-wide: divergent mean {md:.3f}, conserved mean {mc:.3f} "
          f"(epitopes below the divergent mean favour conservation)")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
