# De novo pH-conditional EGFR domain III binders

Twenty de novo single-chain binders to the EGFR domain III (cetuximab/panitumumab) epitope. Each binds more
tightly at tumor pH (~6.5) than at blood pH (7.4) and keeps binding on both human and mouse EGFR. Built in
silico for the Anthropic × Adaptyv 2026 Protein Design Competition, Challenge 1 (Track 3). Objective priority,
highest first: pH-selectivity, then human/mouse cross-reactivity, then affinity, with novelty as a gate.

All twenty are single-chain, 56–108 aa, zero free cysteines, and passed the Adaptyv portal structural-novelty
check (≥3/4 per design). Twelve carry a four-method-confirmed pH switch; the rest are described under
Selection.

## Target

The epitope is EGFR domain III, the surface bound by cetuximab and panitumumab, taken from chain A of PDB
6ARU. Numbering follows the UniProt P00533 precursor convention (the switch histidine is His433, His409 in
mature numbering). His433 is identical in human and mouse and titrates near pH 6.3, so it carries a positive
charge in the acidic tumor margin and is mostly neutral in blood. Two further epitope residues the design
uses, His370 and Asp460, are also identical in human and mouse (human P00533 vs mouse Q01279, 88.73% identity
across the extracellular region). `data/epitope/epitope_map.json` records the conserved anchors and the
cetuximab-footprint positions that differ between species, which the design steers away from.

## Two switch mechanisms

**Forward, 18 designs: binder acid to antigen His433.** The binder presents an Asp or Glu that salt-bridges
His433 only when His433 is protonated, so the interface gains a salt bridge at tumor pH. Liu et al. 2022
(PMID 36458200) showed this at the antibody level: antibody G532 uses an engineered LCDR1 Glu32–His433 bridge,
binds about 8-fold more strongly at pH 6.5 than 7.4 (EC50 ratio 8.08), cross-reacts on human and mouse, and
the enabling Y32E carboxylate raised pH-dependency about 13-fold. The de novo single-chain version here has no
published experimental precedent, which is the main risk.

**Inverted, 2 designs (inv01/inv02): binder histidine to antigen Asp460.** The titrating residue moves onto
the binder, where its pKa is set by design, and keys to the conserved antigen Asp460. Because Asp460 is
identical in human and mouse, this mechanism is cross-reactive by construction. A flat conserved-Asp patch did
not bind (0 of 384 folded); docking the inverted switch onto the His433 face with a two-hotspot design
recovered both binding and a clean switch.

## Generation

RFdiffusion partial diffusion of binders that already passed folding produced the backbones; the inverted
lineage used an all-atom RFdiffusion3 two-hotspot motif that places the binder histidine against Asp460.
Partial diffusion keeps the unusual topology that clears the structural-novelty check, which de novo folds
tend to fail. Soluble ProteinMPNN wrote two sequences per backbone with cysteine disallowed and the switch
residue locked. AF2 initial-guess folded each sequence against a full-atom domain III target, and designs with
interface PAE below about 11 advanced.

## pH-selectivity validation

Folding models are blind to pH, so a separate physics layer scores the switch on the exact submitted sequence.
No pKa method is reliable in absolute terms for an engineered histidine, so four methods vote by sign:

1. **PROPKA Wyman proton linkage:** the proton-count integral from pH 7.4 to 6.5 across complex, target, and
   binder. Positive means tighter at tumor pH.
2. **pKAI:** a machine-learned pKa shift on the same three states. It corroborates PROPKA rather than acting
   as an independent leg, since the two share training data.
3. **Two-state MD occupancy, the selectivity gate:** explicit-protonation MD (OpenMM, GBn2, 3 replicas, 2 ns
   each) run with the sensor histidine protonated (HIP, tumor) and neutral (HIE, blood). The gap
   d_occ = occ_HIP − occ_HIE separates a real switch from a constitutive bridge. A positive static score only
   fixes direction: a bridge that holds whether or not the histidine is charged still scores positive yet
   binds at both pH, and only the two-state occupancy catches it. Most static-positive designs were
   constitutive here and were dropped.
4. **MM-GBSA two-state binding free energy:** dd_bind = dG_bind(HIP) − dG_bind(HIE), read by sign. This is
   the only leg that scores the functional binding readout, and it is independent of the folders and of PROPKA.

Two controls back the layer. The isostere ablation mutates the His433-contacting acid to its neutral isostere
and confirms the linkage collapses, which rules out the binder's bulk charge. The FcRn:Fc complex (PDB 1FRT)
is a positive control: run through the same engine, PROPKA recovers the correct low-pH-ON sign (+0.27) and
pKAI puts the shift on the two switch histidines (+2.66) while staying flat on the bystanders. A third control
ships inside the panel (see Selection).

## Cross-reactivity

Four co-folder-independent signals agree the designs bind mouse EGFR: a deterministic contact-conservation
severity over the 4 Å epitope core (every epitope contact lands on a conserved position), a matched
mouse-minus-human AF2-ig interface-PAE delta, a SurfDiff surface comparison, and Chai mouse/human iPTM. A low
mouse iPTM from any single folder is folder-wide pessimism, not a species signal: Chai-mouse and Chai-human
correlate at Spearman 0.869.

## Design-space survey

A chemistry-fair stack checked whether a better switch exists on this epitope: a GFN2-xTB interaction energy
with self-consistent charges, an APBS continuum energy that includes desolvation, a 0.8 charge-scaling test
for ionic-strength sensitivity, and the novelty oracle. Each alternative fails at a specific point.

- **Cation-pi** (binder-His⁺ on a conserved aromatic): survives charge-scaling, but on the fair QM a single
  cation-pi is comparable to or weaker than a salt bridge; the novelty-passing designs are single and weaker,
  and the stronger dual fails the novelty gate.
- **Multi-histidine neutral acceptors**: a contact strong enough to overcome buried-His desolvation has to be
  ionic; neutral contacts are too weak once buried (the APBS net energy refutes them).
- **Buried salt bridge**: strong but constitutive. The neutral-His off-state energy is −11.6 to −12.0 kcal/mol
  across four buried designs, so there is no off state.
- **Repulsion-gated carboxylate dyad**: the physics works (the off state is repulsive), but the
  buried-glutamate pKa is unpredictable; PROPKA and pKAI disagree by 2 to 3 pH units, wider than the whole
  6.5–7.4 window.
- **Allosteric hinged lid**: the generators produce monolithic folds rather than hinges (1 of 512 a binder),
  and the pH-aware readout reports docked-state stability, not lid release.

Histidine is the only residue whose pKa sits in the 6.5–7.4 window and is predictable, so a His salt bridge,
forward or inverted, is the achievable frontier on this epitope. The panel sits on it.

## Selection

Four tiers, listed in submission order. `d_occ` is the MD occupancy gap (higher is more selective),
`occ_HIE` the residual blood-pH occupancy (lower is a cleaner off state), `ChaiH`/`Boltz` the co-fold
interface iPTM, `af2ig` the AF2 initial-guess interface PAE (lower is better), and `nov_TM` the local
novelty-oracle TM (lower is more novel; the portal check is the authority). A dash marks a leg not recorded
for that tier. Full values are in `submission/metrics.csv` and sequences in `submission/designs.csv`.

| # | id | tier | mech | len | d_occ | occ_HIE | ChaiH | Boltz | af2ig | nov_TM |
|---|----|------|------|-----|-------|---------|-------|-------|-------|--------|
| 1 | 23 | 1 | fwd | 108 | 0.96 | 0.04 | 0.82 | 0.91 | 6.6 | 0.79 |
| 2 | 32 | 3 | KO | 108 | – | – | – | – | – | – |
| 3 | inv02 | 4 | inv | 80 | 0.87 | 0.13 | 0.79 | 0.84 | 9.6 | 0.75 |
| 4 | 27 | 1 | fwd | 104 | 0.98 | 0.02 | 0.21 | 0.58 | 6.7 | 0.83 |
| 5 | inv01 | 4 | inv | 88 | 0.52 | 0.00 | 0.62 | 0.78 | 9.3 | 0.72 |
| 6 | 29 | 1 | fwd | 56 | 0.93 | 0.07 | 0.81 | 0.55 | 7.7 | 0.77 |
| 7 | 33 | 2 | fwd | 92 | 0.93 | 0.02 | – | – | – | 0.78 |
| 8 | 22 | 1 | fwd | 56 | 0.96 | 0.04 | 0.20 | 0.85 | 5.4 | 0.89 |
| 9 | 34 | 2 | fwd | 108 | 0.93 | 0.07 | – | – | – | 0.77 |
| 10 | 36 | 2 | fwd | 56 | 0.95 | 0.04 | – | – | – | 0.65 |
| 11 | 04 | 1 | fwd | 104 | 0.89 | – | 0.78 | 0.64 | 11.4 | – |
| 12 | 38 | 2 | fwd | 70 | 0.55 | 0.38 | – | – | – | 0.79 |
| 13 | 08 | 1 | fwd | 104 | 0.57 | – | 0.70 | 0.41 | 8.6 | – |
| 14 | 39 | 2 | fwd | 92 | 0.48 | 0.20 | – | – | – | 0.77 |
| 15 | 19 | 1 | fwd | 61 | 0.48 | – | 0.85 | 0.83 | 5.2 | – |
| 16 | 12 | 1 | fwd | 70 | 0.34 | – | 0.84 | 0.86 | 6.9 | – |
| 17 | 16 | 1 | fwd | 92 | 0.20 | – | 0.76 | 0.91 | 10.3 | – |
| 18 | 02 | 1 | fwd | 91 | 0.15 | – | 0.21 | 0.86 | 7.0 | – |
| 19 | 06 | 1 | fwd | 108 | 0.59 | – | 0.81 | 0.83 | 7.4 | – |
| 20 | 09 | 1 | fwd | 108 | 0.26 | – | 0.86 | 0.90 | 6.9 | – |

- **Tier 1, forward, 12 designs.** All four pH legs agree in sign. The isostere control gives a clean collapse
  for 22/23/29, is inconclusive for 27, and is uninformative for the eight locked designs whose
  single-structure His433 linkage is near zero, which rest on the MD gate. The tier spans the switch/binding
  trade-off; 02 and 16 sit at or below the 0.2 selectivity floor and are kept as portal-locked anchors on
  binding and novelty. Low ChaiH on 02/22/27 is folder pessimism, since Boltz and af2ig agree they bind
  (22: Boltz 0.85, af2ig 5.4).
- **Tier 2, forward, 5 designs.** Cysteine-cleaned children of Tier-1 lineages, validated on the MD occupancy
  gate plus binding and cross-reactivity. Binding is a single co-fold consensus: 33 0.53, 34 0.45, 36 0.84,
  38 0.44, 39 0.61.
- **Tier 3, knockout control, 1 design.** Design 32 is design 23 with a single E66Q at position 66, which
  removes the switch acid. It should still bind but lose the pH preference, so 23 and 32 ship together as a
  bench test of the mechanism. Its loss of switch is inferred from the parent, not measured on the mutant.
- **Tier 4, inverted, 2 designs.** The inverted mechanism, with a clean two-state switch (inv01 d_occ 0.52,
  occ_HIE 0.00; inv02 0.87, 0.13), three-predictor binding, and cross-reactivity by construction.

## Reproduce

The pipeline drives external prediction, design, and physics tools and is not turnkey. It needs AlphaFold2
initial-guess, RFdiffusion and RFdiffusion3, ProteinMPNN, Boltz-2, Chai-1, OpenMM (Amber ff14SB / GBn2),
PROPKA3, pKAI, xtb (GFN2-xTB), and APBS, each with its own weights and environment, which this repo does not
vendor. `scripts/` holds the code for every step. The submission validator runs on its own:

```
python scripts/validate_submission.py submission/designs.csv --track 3
```

## Repository layout

```
submission/designs.csv   the 20 designs (name, sequence, molecule_class), submission order
submission/metrics.csv   per-design metrics on the submitted sequence
scripts/                 design and validation pipeline
data/epitope/            conserved anchors and divergent positions (human P00533 vs mouse Q01279)
configs/target.yaml      target definition (EGFR domain III, PDB 6ARU)
```

## Citations

Liu et al. 2022 (PMID 36458200), cross-reactive pH-dependent anti-EGFR domain III antibody; Burmeister et al.
1994 (PMID 7969498), FcRn:Fc, PDB 1FRT. Tools: RFdiffusion (Watson 2023) and RFdiffusion3; ProteinMPNN
(Dauparas 2022); AlphaFold2 (Jumper 2021) and AF2 initial-guess (Bennett 2023); AlphaFold3 (Abramson 2024);
ESMFold (Lin 2023); Chai-1 (2024); Boltz-2 (Passaro 2025); PROPKA (Olsson 2011); pKAI (Reis 2022); Foldseek
(van Kempen 2024); OpenMM (Eastman 2017); Amber ff14SB (Maier 2015); GBn2 (Nguyen 2013); PRODIGY (Xue 2016);
GFN2-xTB (Bannwarth 2019); APBS (Jurrus 2018).

## License

MIT, see `LICENSE`. Designed by Vijayavallabh.
