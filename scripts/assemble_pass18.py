#!/usr/bin/env python3
"""Portal-confirmed 18-design submission (2026-10-05 novelty check of the W21 panel).

The W21 panel (20 designs) was run through the Adaptyv portal novelty check. 18 scored 3/4
(PASS); exactly two scored 2/4 (FAIL): egfr_phsw_35 (lineage 09) and egfr_phsw_37 (lineage 01
double-partial-diffusion). The portal rule is "remove designs below 3/4 to submit", so this file is
the W21 panel minus those two, in the same best-of order. Every submitted design is portal-verified
3/4 -> the "all submitted designs >= 3/4" gate is satisfied with certainty, zero further iteration.

Each sequence and its molecule_class are carried over from the W21 CSV unchanged (only surrounding
whitespace is stripped; the W21 sequences are already clean uppercase), so the 18 stay byte-identical
to the W21 designs the portal scored 3/4. Track 3 allows up to 20; 18 all-confirmed is valid and keeps
the full strong front tier: the top of the order is the highest-d_occ switches (0.927-0.978), the two
severe-expression designs are parked in the tail, and the knockout control is last.
"""
import csv, sys, subprocess

REPO = "."
W21 = f"{REPO}/results/adaptyv_challenge1_resubmit_W21.csv"
OUT = f"{REPO}/results/adaptyv_challenge1_resubmit_PASS18.csv"

# Portal 2026-10-05: these two W21 designs scored 2/4 and are removed.
PORTAL_FAIL = {"egfr_phsw_35", "egfr_phsw_37"}

# The 18 portal-confirmed (3/4) designs, best-of order (= W21 row order minus the 2 fails). The
# explicit list documents the rationale per design; main() asserts it equals the derived W21-row-order
# -minus-fails, so it cannot silently drift if W21 is ever regenerated in a different order.
ORDER = [
    "egfr_phsw_29",  # d_occ 0.927, CLEAN, confirmed 3/4
    "egfr_phsw_23",  # d_occ 0.957; PARENT of the egfr_phsw_32 knockout (E66Q) -> keep in top tier
    "egfr_phsw_33",  # refill, d_occ 0.932, CLEAN -> now portal-confirmed 3/4
    "egfr_phsw_27",  # d_occ 0.978 (top switch)
    "egfr_phsw_22",  # d_occ 0.962 (strong lin01 switch)
    "egfr_phsw_34",  # refill, d_occ 0.927 -> confirmed 3/4
    "egfr_phsw_36",  # refill, d_occ 0.952, near-twin of 29, CLEAN -> confirmed 3/4
    "egfr_phsw_04",  # locked de novo, d_occ 0.89
    "egfr_phsw_38",  # refill, d_occ 0.55, CLEAN -> confirmed 3/4
    "egfr_phsw_08",  # locked de novo, d_occ 0.57, CLEAN
    "egfr_phsw_39",  # refill, d_occ 0.485, CLEAN -> confirmed 3/4
    "egfr_phsw_19",  # locked de novo, d_occ 0.48, CLEAN
    "egfr_phsw_12",  # locked de novo, d_occ 0.335, CLEAN
    "egfr_phsw_16",  # locked de novo, d_occ 0.195, CLEAN
    "egfr_phsw_02",  # locked de novo, d_occ 0.155, CLEAN
    "egfr_phsw_06",  # locked de novo, SEVERE expr -> tail
    "egfr_phsw_09",  # locked de novo (lin09 SEED; the seed passed novelty, its pd-children 24/35 failed). SEVERE expr -> tail
    "egfr_phsw_32",  # knockout control (E66Q of egfr_phsw_23) -> last
]


def main():
    with open(W21) as fh:
        w21_rows = list(csv.DictReader(fh))
    w21 = {r["name"]: r for r in w21_rows}
    w21_order = [r["name"] for r in w21_rows]

    # hard invariants
    assert set(ORDER) <= set(w21), "a confirmed design is missing from W21"
    assert len(ORDER) == 18 and len(set(ORDER)) == 18, "expected 18 distinct designs"
    assert PORTAL_FAIL.isdisjoint(ORDER), "a portal-failed design leaked into the order"
    assert set(w21) - set(ORDER) == PORTAL_FAIL, "dropped set != the 2 portal fails"
    # order cannot silently diverge from 'W21 row order minus the 2 fails' (best-of order drives selection)
    expected = [n for n in w21_order if n not in PORTAL_FAIL]
    assert ORDER == expected, f"ORDER diverged from W21 row order minus fails:\n have {ORDER}\n want {expected}"

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for n in ORDER:
            row = w21[n]
            seq = row["sequence"].strip()  # carry over verbatim (already clean uppercase)
            assert "C" not in seq.upper(), f"{n} has a cysteine"
            w.writerow([n, seq, (row.get("molecule_class") or "single_chain").strip()])
    print(f"wrote {OUT}: 18 portal-confirmed (3/4) designs, all 0-Cys, guaranteed all-pass")
    print(f"dropped (portal 2/4): {sorted(PORTAL_FAIL)}")

    # validation gate: fail loud if the official validator rejects (don't claim all-pass on exit 0)
    print("\n=== validate_submission.py (track 3) ===")
    rc = subprocess.run([f"{REPO}/envs/md/bin/python", f"{REPO}/scripts/validate_submission.py",
                         OUT, "--track", "3"]).returncode
    if rc != 0:
        sys.exit(f"[FATAL] validate_submission.py FAILED (exit {rc}) -- submission is NOT valid")
    print("validator PASS")


if __name__ == "__main__":
    main()
