"""adaptyv_egfr - submission tooling for the Anthropic x Adaptyv 2026 EGFR
conditional-binder challenge (Challenge 1).

Objective priority (this order IS the scoring priority):
    1. pH-selective binding        bind @ pH 6.5, not @ pH 7.4   (highest)
    2. human/mouse cross-reactive   the same sequence binds both
    3. high affinity to human EGFR extracellular region

This package ships the submission-format validator (``submission.py``). The design
and validation pipeline lives in ``scripts/``; see ``docs/METHODS.md``.
"""
from __future__ import annotations

__version__ = "1.0.0"

# Objective keys in priority order.
OBJECTIVES = ("ph_selective", "cross_reactive", "human_affinity")
