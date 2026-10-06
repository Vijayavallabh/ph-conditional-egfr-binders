#!/usr/bin/env python3
"""Definition-of-done gate: validate a candidate submission CSV against the
Anthropic x Adaptyv 2026 EGFR challenge contract. Exit 0 = PASS, 1 = FAIL.

Usage:
    uv run python scripts/validate_submission.py submissions/track3.csv --track 3
    uv run python scripts/validate_submission.py file.csv --strict     # warnings -> failure
    uv run python scripts/validate_submission.py file.csv --no-anarci   # skip Ab numbering
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Work even without an editable install of the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from adaptyv_egfr.submission import format_report, validate_csv  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Validate an EGFR submission CSV.")
    ap.add_argument("csv", type=Path, help="candidate submission .csv")
    ap.add_argument("--track", type=int, choices=[1, 2, 3], default=3,
                    help="per-track design cap (Track 1: 40, Tracks 2/3: 20); this repo defaults to Track 3")
    ap.add_argument("--strict", action="store_true", help="treat warnings as failure")
    ap.add_argument("--no-anarci", action="store_true", help="skip antibody numbering checks")
    a = ap.parse_args(argv)

    rep = validate_csv(a.csv, track=a.track, use_anarci=not a.no_anarci)
    print(format_report(rep, a.csv))
    failed = (not rep.ok) or (a.strict and rep.n_warnings > 0)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
