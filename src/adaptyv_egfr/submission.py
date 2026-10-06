"""Submission contract + validator for the Anthropic x Adaptyv 2026 EGFR challenge.

Single source of truth for the submission CSV format, used by
``scripts/validate_submission.py``.

CONTRACT (from the Adaptyv Challenge 1 submission page). The live /submit form is a
client-rendered SPA, so re-verify exact column names/casing on the logged-in form
before the final upload.

  * one CSV; rows ordered best-first (top row = highest-ranked design)
  * required columns: name, sequence, molecule_class
  * molecule_class in {single_chain, nanobody, scfv, fab_kappa, fab_lambda} (lowercase;
    Adaptyv renamed the single-chain class 'protein' -> 'single_chain' on the live portal,
    2026-10 — both are accepted here, 'single_chain' is the current portal value)
  * sequence: natural amino acids only; Fab formatted "VH:VL" (one ':' separator)
  * length 10-250 aa inclusive (per chain for a Fab)
  * names unique; sequences unique (de-novo, de-duplicated)
  * nanobody/scfv/fab must pass ANARCI-style numbering; fab_kappa/fab_lambda VL
    isotype (K/L) must match the declared class
  * per-track caps: Track 1 <= 40 designs; Tracks 2 & 3 <= 20 designs each
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

REQUIRED_COLUMNS: tuple[str, ...] = ("name", "sequence", "molecule_class")
MOLECULE_CLASSES: tuple[str, ...] = ("single_chain", "protein", "nanobody", "scfv", "fab_kappa", "fab_lambda")
# Adaptyv renamed the single-chain class 'protein' -> 'single_chain' on the live portal (2026-10).
# 'single_chain' is the current portal value; 'protein' kept as a legacy alias so prior CSVs still pass.
SINGLE_CHAIN_CLASSES: frozenset[str] = frozenset({"single_chain", "protein"})
FAB_CLASSES: frozenset[str] = frozenset({"fab_kappa", "fab_lambda"})
ANTIBODY_CLASSES: frozenset[str] = frozenset({"nanobody", "scfv", "fab_kappa", "fab_lambda"})
AA: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")
LEN_MIN: int = 10
LEN_MAX: int = 250
TRACK_CAPS: dict[int, int] = {1: 40, 2: 20, 3: 20}
# Soft warning window for a plausible immunoglobulin V-domain length.
VDOMAIN_MIN, VDOMAIN_MAX = 90, 150


def length_category(n: int) -> str:
    """Scoring sub-category by single-chain length (per the challenge)."""
    if n < 40:
        return "protein_microbinder"
    if n <= 100:
        return "protein_minibinder"
    return "large_protein"


# --------------------------------------------------------------------------- #
# Reports
# --------------------------------------------------------------------------- #
@dataclass
class RowReport:
    index: int
    name: str
    molecule_class: str
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    info: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not self.errors


@dataclass
class Report:
    rows: list[RowReport] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)  # file-level
    warnings: list[str] = field(default_factory=list)  # file-level
    n_rows: int = 0
    anarci_available: bool = False

    @property
    def ok(self) -> bool:
        return not self.errors and all(r.ok for r in self.rows)

    @property
    def n_errors(self) -> int:
        return len(self.errors) + sum(len(r.errors) for r in self.rows)

    @property
    def n_warnings(self) -> int:
        return len(self.warnings) + sum(len(r.warnings) for r in self.rows)


# --------------------------------------------------------------------------- #
# anarcism (ANARCI-equivalent) integration - imported lazily & defensively
# --------------------------------------------------------------------------- #
def anarci_available() -> bool:
    try:
        import anarcism  # noqa: F401

        return True
    except Exception:
        return False


def _domains(seq: str):
    """List of anarcism DomainResult for ``seq`` (empty list on parse-miss, None on failure)."""
    try:
        import anarcism

        return list(anarcism.number_sequence(seq).domains)
    except Exception:
        return None


def _chain_type(dom) -> str:
    return getattr(dom, "chain_type", "?")


def _clean_seq(s: str) -> str:
    return "".join((s or "").split()).upper()


def _bad_chars(seq: str) -> set[str]:
    return {c for c in seq if c not in AA}


# --------------------------------------------------------------------------- #
# Row-level validation
# --------------------------------------------------------------------------- #
def validate_row(row: dict, *, use_anarci: bool = True) -> RowReport:
    idx = int(row.get("_index", 0))
    name = (row.get("name") or "").strip()
    raw_class = (row.get("molecule_class") or "").strip()
    rep = RowReport(index=idx, name=name, molecule_class=raw_class)

    if not name:
        rep.errors.append("empty 'name'")
    if "@" in name:
        rep.warnings.append("'name' contains '@' - names are published verbatim; avoid PII")

    mclass = raw_class
    if mclass not in MOLECULE_CLASSES:
        if mclass.lower() in MOLECULE_CLASSES:
            rep.warnings.append(f"molecule_class casing: use lowercase '{mclass.lower()}'")
            mclass = mclass.lower()
        else:
            rep.errors.append(f"molecule_class '{raw_class}' not in {list(MOLECULE_CLASSES)}")
            return rep
    rep.molecule_class = mclass

    seq = _clean_seq(row.get("sequence"))
    if not seq:
        rep.errors.append("empty 'sequence'")
        return rep

    is_fab = mclass in FAB_CLASSES
    ncolon = seq.count(":")
    if is_fab:
        if ncolon != 1:
            rep.errors.append(f"Fab must be 'VH:VL' with exactly one ':' (found {ncolon})")
            return rep
        vh, vl = seq.split(":")
        chains = {"VH": vh, "VL": vl}
    else:
        if ncolon != 0:
            rep.errors.append(f"'{mclass}' is single-chain; ':' not allowed")
            return rep
        chains = {"chain": seq}

    for label, chain in chains.items():
        bad = _bad_chars(chain)
        if bad:
            rep.errors.append(f"{label}: non-standard residues {sorted(bad)}")
        n = len(chain)
        if not (LEN_MIN <= n <= LEN_MAX):
            rep.errors.append(f"{label}: length {n} outside {LEN_MIN}-{LEN_MAX} aa")
        elif is_fab and not (VDOMAIN_MIN <= n <= VDOMAIN_MAX):
            rep.warnings.append(f"{label}: length {n} unusual for a V-domain (~{VDOMAIN_MIN}-{VDOMAIN_MAX})")

    if not is_fab:
        rep.info["length"] = len(seq)
        if mclass in SINGLE_CHAIN_CLASSES:
            rep.info["category"] = length_category(len(seq))

    if use_anarci and mclass in ANTIBODY_CLASSES:
        if not anarci_available():
            rep.warnings.append("anarcism not installed - antibody numbering/isotype unchecked")
        elif not rep.errors:  # only number well-formed chains
            _check_antibody(rep, mclass, chains)
    return rep


def _check_antibody(rep: RowReport, mclass: str, chains: dict) -> None:
    import anarcism

    if mclass == "nanobody":
        doms = _domains(chains["chain"]) or []
        cts = [_chain_type(d) for d in doms]
        rep.info["ig_domains"] = cts
        if cts != ["H"]:
            rep.errors.append(f"nanobody must be a single VHH (H) domain; ANARCI found {cts or 'none'}")
    elif mclass == "scfv":
        doms = _domains(chains["chain"]) or []
        cts = [_chain_type(d) for d in doms]
        rep.info["ig_domains"] = cts
        if not ("H" in cts and ("K" in cts or "L" in cts)):
            rep.errors.append(f"scFv must contain one H and one K/L domain; ANARCI found {cts or 'none'}")
    else:  # fab_kappa / fab_lambda
        try:
            pr = anarcism.validate_antibody_pair(chains["VH"], chains["VL"])
        except Exception as e:  # pragma: no cover - defensive
            rep.errors.append(f"ANARCI pair validation errored: {e}")
            return
        rep.info["ig_domains"] = [_chain_type(pr.vh), _chain_type(pr.vl)]
        if not pr.ok:
            rep.errors.append(f"ANARCI: VH/VL pair invalid: {list(pr.errors)}")
        if _chain_type(pr.vh) != "H":
            rep.errors.append(f"Fab VH is not a heavy domain (got {_chain_type(pr.vh)})")
        want = "K" if mclass == "fab_kappa" else "L"
        got = _chain_type(pr.vl)
        if got != want:
            rep.errors.append(f"molecule_class '{mclass}' but VL isotype is '{got}' (expected '{want}')")


# --------------------------------------------------------------------------- #
# Set-level validation
# --------------------------------------------------------------------------- #
def validate_rows(rows: Iterable[dict], *, track: int | None = None, use_anarci: bool = True) -> Report:
    rep = Report(anarci_available=anarci_available())
    rows = list(rows)
    rep.n_rows = len(rows)
    if not rows:
        rep.errors.append("no data rows")
        return rep

    for i, r in enumerate(rows, start=1):
        rr = dict(r)
        rr["_index"] = i
        rep.rows.append(validate_row(rr, use_anarci=use_anarci))

    # duplicate names
    by_name: dict[str, list[int]] = {}
    for rr in rep.rows:
        by_name.setdefault(rr.name, []).append(rr.index)
    for nm, idxs in by_name.items():
        if nm and len(idxs) > 1:
            rep.errors.append(f"duplicate name '{nm}' at rows {idxs}")

    # duplicate sequences (uniqueness hard filter)
    by_seq: dict[str, list[int]] = {}
    for i, r in enumerate(rows, start=1):
        by_seq.setdefault(_clean_seq(r.get("sequence")), []).append(i)
    for s, idxs in by_seq.items():
        if s and len(idxs) > 1:
            rep.errors.append(f"duplicate sequence at rows {idxs} (de-duplicate before submitting)")

    # per-track cap
    cap = TRACK_CAPS.get(track) if track else None
    if cap is not None and len(rows) > cap:
        rep.errors.append(f"{len(rows)} designs exceeds Track {track} cap of {cap}")
    elif cap is None and len(rows) > max(TRACK_CAPS.values()):
        rep.warnings.append(
            f"{len(rows)} designs exceeds the largest cap ({max(TRACK_CAPS.values())}); pass --track"
        )

    # optional ordering sanity
    if rows and "rank" in rows[0]:
        try:
            ranks = [float(r["rank"]) for r in rows]
            if ranks != sorted(ranks):
                rep.warnings.append("'rank' column is not ascending; rows must be best-first")
        except Exception:
            rep.warnings.append("'rank' column present but not numeric")
    return rep


def validate_csv(path, *, track: int | None = None, use_anarci: bool = True) -> Report:
    p = Path(path)
    if not p.exists():
        rep = Report()
        rep.errors.append(f"file not found: {p}")
        return rep
    with p.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]
    missing = [c for c in REQUIRED_COLUMNS if c not in headers]
    rep = validate_rows(rows, track=track, use_anarci=use_anarci)
    if missing:
        rep.errors.insert(0, f"missing required column(s): {missing}; found {headers}")
    return rep


def write_submission(rows: Iterable[dict], path, *, extra_columns: Iterable[str] | None = None) -> Path:
    """Write a ranked submission CSV (required columns first, in given best-first order)."""
    cols = list(REQUIRED_COLUMNS) + [c for c in (extra_columns or []) if c not in REQUIRED_COLUMNS]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    return p


def format_report(rep: Report, path=None) -> str:
    lines = [f"Submission validation: {path}" if path else "Submission validation"]
    lines.append(
        f"rows={rep.n_rows}  errors={rep.n_errors}  warnings={rep.n_warnings}  "
        f"anarci={'yes' if rep.anarci_available else 'no'}"
    )
    for e in rep.errors:
        lines.append(f"  [FILE ERROR] {e}")
    for w in rep.warnings:
        lines.append(f"  [FILE WARN ] {w}")
    for rr in rep.rows:
        if rr.errors or rr.warnings:
            lines.append(f"  row {rr.index} '{rr.name}' ({rr.molecule_class}):")
            for e in rr.errors:
                lines.append(f"      ERROR  {e}")
            for w in rr.warnings:
                lines.append(f"      warn   {w}")
    lines.append("RESULT: " + ("PASS" if rep.ok else "FAIL"))
    return "\n".join(lines)
