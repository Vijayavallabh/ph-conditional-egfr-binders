# De novo pH-conditional EGFR domain III binders

Twenty de novo single-chain protein binders to the EGFR domain III (cetuximab/panitumumab) epitope,
engineered to bind **more tightly at tumor pH (~6.5) than at blood pH (7.4)** and to cross-react on **human
and mouse** EGFR. Anthropic × Adaptyv 2026 Protein Design Competition, Challenge 1 (Track 3), in-silico.
Objective priority: **pH-selectivity > human/mouse cross-reactivity > affinity**, with novelty as a gate.

The design switch is a single histidine titration near the tumor/blood boundary, realized by **two
mechanistically independent mechanisms**:

1. **Binder acid → antigen His433** (18 designs). A binder Asp/Glu salt-bridges EGFR His433 only when His433
   is protonated at acidic pH — a de novo single-chain realization of the Liu-2022 antibody mechanism
   (G532, PMID 36458200).
2. **Binder histidine → conserved antigen Asp460** (2 designs). The switch is inverted onto the binder, so
   the titrating residue sits where its pKa is set by design, keyed to the human/mouse-conserved Asp460.

## Results at a glance

- 2 independent switch mechanisms (neither has a published de novo single-chain precedent)
- 12/20 four-leg MD-gated pH switches (PROPKA + pKAI + two-state MD occupancy + two-state MM-GBSA, sign-consensus)
- 20/20 passed the Adaptyv portal structural-novelty check (≥3/4)
- 20/20 cross-reactive on a four-signal, co-folder-independent consensus
- 20/20 zero free cysteines, single chain, 56–108 aa
- Built-in controls: an FcRn:Fc positive control for the pH engine, and a one-residue E66Q knockout shipped
  beside its parent as a bench falsification test

## Approach

1. **Generation.** RFdiffusion partial diffusion of folding-validated binders (and an all-atom
   RFdiffusion3 two-hotspot motif for the inverted mechanism) → soluble ProteinMPNN (cysteine disallowed,
   switch residue locked) → AF2 initial-guess self-consistency against a full-atom domain III target.
2. **pH-selectivity validation** on the exact submitted sequence (the pH-blind folders cannot judge this):
   PROPKA Wyman proton-linkage, the pKAI machine-learned pKa shift, **explicit-protonation two-state MD
   occupancy** (the selectivity gate that separates a true switch from a constitutive always-on bridge), and
   two-state MM-GBSA — read by sign-consensus. Externally sanity-checked on the natural FcRn:Fc histidine
   switch (PDB 1FRT).
3. **Cross-reactivity.** Deterministic contact-conservation severity + matched mouse−human co-fold delta +
   SurfDiff surface comparison + Chai mouse/human iPTM.
4. **Design-space survey (method novelty).** A chemistry-fair GFN2-xTB + APBS + charge-scaling + novelty
   scoring stack showing that His-based salt bridges are the achievable pH-switch frontier on this epitope.
5. **Selection.** Four tiers spanning both mechanisms and the switch/binding trade-off, lineage-capped, with
   the knockout control.

See **`docs/METHODS.md`** for the full methodology and **`docs/DESIGN_TABLE.md`** for per-design metrics.

## Layout

```
submission/designs.csv    the 20 submitted designs (name, sequence, molecule_class), mechanism-forward order
submission/metrics.csv    per-design metrics on the exact submitted sequence
docs/METHODS.md           full methods (target, mechanisms, validation, frontier survey, selection)
docs/DESIGN_TABLE.md      per-design metric table
docs/SUBMISSION_METHODOLOGY.md   short methodology statement
scripts/                  the design + validation pipeline (see below)
data/epitope/epitope_map.json    conserved anchors / divergent positions (human P00533 vs mouse Q01279)
configs/target.yaml       target definition (EGFR domain III, PDB 6ARU)
```

## Reproduction

The pipeline orchestrates external structure-prediction, design, and physics tools and is **not turnkey**:
it requires AlphaFold2 initial-guess, RFdiffusion / RFdiffusion3, ProteinMPNN, Boltz-2, Chai-1, OpenMM
(Amber ff14SB / GBn2), PROPKA3, pKAI, GFN2-xTB (xtb), and APBS, with their own weights and environments,
which are not vendored here. `scripts/` contains the code for each step; `docs/METHODS.md` → "Reproduce"
maps each command to its script. `scripts/validate_submission.py` checks the submission CSV standalone:

```
python scripts/validate_submission.py submission/designs.csv --track 3
```

## Attribution

Designed by Vijayavallabh for the Anthropic × Adaptyv 2026 Protein Design Competition, Challenge 1.
Precedent and method citations are listed in `docs/METHODS.md`.
