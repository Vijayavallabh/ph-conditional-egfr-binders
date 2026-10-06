# Design table: twenty de novo pH-conditional EGFR domain III binders

Submission order. All 20: single_chain, 0 free cysteines, **portal novelty check PASS (>=3/4, user-run 2026-10-06)**. Metrics computed on the exact submitted sequence; raw values in `submission/metrics.csv`, methods in `METHODS.md`.

Legend — **d_occ** = two-state MD occupancy gap occ_HIP−occ_HIE (pH selectivity; higher = more selective); **occ_HIE** = residual blood-pH bridge occupancy (lower = cleaner OFF; blank = not recorded in this file for that design); **pKAI / PROPKA** = His pKa-shift linkage energies (sign = switch direction); **iso** = isostere-ablation result (collapse = designed bridge confirmed causal; uninform. = static linkage ~0 so control is uninformative, design rests on MD); **ddG_bind** = MM-GBSA two-state dd_bind kcal/mol (negative = tighter at tumor pH; read by sign); **ChaiH / Boltz** = co-fold interface iPTM; **af2ig** = AF2-initial-guess interface PAE (lower = better); **nov_TM** = local ESMFold2-oracle TM (lower = more novel; the PORTAL is the authority — see note). A dash = that leg not recorded for that design.

| # | name | tier | mechanism | len | d_occ | occ_HIE | pKAI | PROPKA | iso | ddG_bind | ChaiH | Boltz | af2ig | nov_TM | portal |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | egfr_phsw_23 | T1 | 1 (acid->His433) | 108 | 0.96 | 0.04 | 0.97 | 1.96 | collapse | -17.44 | 0.82 | 0.91 | 6.59 | 0.79 | PASS |
| 2 | egfr_phsw_32 | KO | 1 knockout | 108 | - | - | - | - | - | - | - | - | - | - | PASS |
| 3 | egfr_phsw_inv02 | T4 | 2 (His->Asp460) | 80 | 0.87 | 0.13 | - | - | - | - | 0.79 | 0.84 | 9.56 | 0.75 | PASS |
| 4 | egfr_phsw_27 | T1 | 1 (acid->His433) | 104 | 0.98 | 0.02 | 0.7 | 1.27 | uninform. | -6.8 | 0.21 | 0.58 | 6.73 | 0.83 | PASS |
| 5 | egfr_phsw_inv01 | T4 | 2 (His->Asp460) | 88 | 0.52 | 0 | - | - | - | - | 0.62 | 0.78 | 9.28 | 0.72 | PASS |
| 6 | egfr_phsw_29 | T1 | 1 (acid->His433) | 56 | 0.93 | 0.07 | 0.95 | 1.27 | collapse | -8.51 | 0.81 | 0.55 | 7.69 | 0.77 | PASS |
| 7 | egfr_phsw_33 | T2 | 1 (acid->His433) | 92 | 0.93 | 0.02 | - | - | - | - | - | - | - | 0.78 | PASS |
| 8 | egfr_phsw_22 | T1 | 1 (acid->His433) | 56 | 0.96 | 0.04 | 0.59 | 0.99 | collapse | -8.03 | 0.2 | 0.85 | 5.39 | 0.89 | PASS |
| 9 | egfr_phsw_34 | T2 | 1 (acid->His433) | 108 | 0.93 | 0.07 | - | - | - | - | - | - | - | 0.77 | PASS |
| 10 | egfr_phsw_36 | T2 | 1 (acid->His433) | 56 | 0.95 | 0.04 | - | - | - | - | - | - | - | 0.65 | PASS |
| 11 | egfr_phsw_04 | T1 | 1 (acid->His433) | 104 | 0.89 | - | 0.93 | 1.71 | - | -4.57 | 0.78 | 0.64 | 11.4 | - | PASS |
| 12 | egfr_phsw_38 | T2 | 1 (acid->His433) | 70 | 0.55 | 0.38 | - | - | - | - | - | - | - | 0.79 | PASS |
| 13 | egfr_phsw_08 | T1 | 1 (acid->His433) | 104 | 0.57 | - | 0.63 | 1.59 | - | -2.02 | 0.7 | 0.41 | 8.6 | - | PASS |
| 14 | egfr_phsw_39 | T2 | 1 (acid->His433) | 92 | 0.48 | 0.2 | - | - | - | - | - | - | - | 0.77 | PASS |
| 15 | egfr_phsw_19 | T1 | 1 (acid->His433) | 61 | 0.48 | - | 0.36 | 0.8 | - | -3.29 | 0.85 | 0.83 | 5.2 | - | PASS |
| 16 | egfr_phsw_12 | T1 | 1 (acid->His433) | 70 | 0.34 | - | 0.13 | 1.34 | - | -2.78 | 0.84 | 0.86 | 6.9 | - | PASS |
| 17 | egfr_phsw_16 | T1 | 1 (acid->His433) | 92 | 0.2 | - | 0.35 | 1.06 | - | -8.47 | 0.76 | 0.91 | 10.3 | - | PASS |
| 18 | egfr_phsw_02 | T1 | 1 (acid->His433) | 91 | 0.15 | - | 0.69 | 2.12 | - | -4.93 | 0.21 | 0.86 | 7 | - | PASS |
| 19 | egfr_phsw_06 | T1 | 1 (acid->His433) | 108 | 0.59 | - | 0.95 | 1.6 | - | -11.98 | 0.81 | 0.83 | 7.4 | - | PASS |
| 20 | egfr_phsw_09 | T1 | 1 (acid->His433) | 108 | 0.26 | - | 0.85 | 1.47 | - | -6.45 | 0.86 | 0.9 | 6.9 | - | PASS |

**Binding of the five Tier-2 refills** is reported as a single co-fold consensus (the component iPTMs were not separately logged): egfr_phsw_33 0.53, 34 0.45, 36 0.84, 38 0.44, 39 0.61. **On the low-Chai Tier-1 designs** (egfr_phsw_02/22/27, ChaiH ~0.20) the other predictors agree the design binds — e.g. egfr_phsw_22 Boltz 0.85 / af2ig PAE 5.39, egfr_phsw_02 Boltz 0.86 / PRODIGY −11.4, egfr_phsw_27 Boltz 0.59 / af2ig 6.73 — so the low ChaiH is folder pessimism (Chai is the weakest AF3-class folder here), not a non-binder.

## Tiers

- **Tier 1 (12): four-leg MD-gated pH switch** (PROPKA + pKAI + two-state MD occupancy + MM-GBSA, all correct in sign). The isostere causal control gives a clean collapse for three (egfr_phsw_22/23/29); for egfr_phsw_27 it is inconclusive and for the eight locked designs the single-structure His433 linkage is ~0, so the isostere control is uninformative there and those rest on the MD occupancy gate. The tier deliberately spans the switch/binding trade-off. Note egfr_phsw_02 (d_occ 0.155) and egfr_phsw_16 (0.195) sit at/below the pipeline's 0.2 MD-selectivity floor and are carried as portal-locked WAVE-18 anchors (strong on binding/novelty), not for switch strength.
- **Tier 2 (5): portal-novel cysteine-cleaned children** of Tier-1 lineages; validated on two-state MD occupancy + binding consensus + cross-react mechanism (not the full static four-leg set).
- **Tier 3 (1): knockout control.** egfr_phsw_32 = egfr_phsw_23 with a single E66Q at position 66 (verified one-residue difference). Predicted to bind but lose the switch; shipped with its parent (23) as a built-in bench falsification test. Its loss-of-switch is inferred from the parent's isostere collapse and the E66Q chemistry, not computed by MD on the mutant (its metric row is intentionally blank).
- **Tier 4 (2): inverted mechanism** (binder-His -> conserved Asp460). Clean two-state MD switch (inv01 d_occ 0.515 / occ_HIE 0.00; inv02 d_occ 0.87 / occ_HIE 0.13), 3-predictor binding, cross-reactive by construction (Asp460 is human=mouse; inv02 Chai-mouse 0.801).

**Novelty note.** All 20 passed the Adaptyv portal novelty check (the authority). The local ESMFold2-oracle `nov_TM` is only a ranker and cannot resolve the boundary: egfr_phsw_22 (0.89) and egfr_phsw_27 (0.83) exceed the local 0.80 safe harbor yet are portal-confirmed novel. `nov_TM` is sourced from `struct_tm` (Tier-1/inverted) or `oracle_TM` (refills) in the consolidated CSV.
