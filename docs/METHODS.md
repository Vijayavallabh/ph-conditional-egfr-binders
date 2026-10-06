# Methods: twenty de novo pH-conditional EGFR domain III binders

Twenty de novo single-chain binders to the clinically validated cetuximab/panitumumab epitope of EGFR domain
III, each engineered to tighten at tumor pH (~6.5) and release at blood pH (7.4) while cross-reacting on human
and mouse EGFR, through **two mechanistically independent pH switches**. Anthropic x Adaptyv 2026, Challenge 1
(Track 3), in-silico only. Objective priority: pH-selectivity first, then human/mouse cross-reactivity, then
affinity, with novelty as a gate. The 20 submitted sequences are in `submission/designs.csv` (ordered
mechanism-forward); per-design metrics in `submission/metrics.csv` and `docs/DESIGN_TABLE.md`. All 20 passed
the Adaptyv portal structural-novelty check (≥3/4 per design).

## Results at a glance

- **2 independent switch mechanisms** (a binder acid to target His433; a binder histidine to the conserved
  target Asp460) — neither has a published de novo single-chain precedent.
- **12 / 20 four-leg MD-gated pH switches** confirmed in sign by all of PROPKA, pKAI, two-state MD occupancy,
  and two-state MM-GBSA; the remaining 8 carry the MD and binding legs (see tiers).
- **20 / 20 passed the Adaptyv portal novelty check** (≥3/4 per design; user-run 2026-10-06).
- **20 / 20 cross-reactive** on a four-signal, co-folder-independent consensus; His433, His370 and Asp460 are
  byte-identical in human and mouse.
- **20 / 20 zero free cysteines, single chain, 56–108 aa.**
- **Two built-in controls:** an FcRn:Fc positive control proving the pH engine recovers a known low-pH-ON
  histidine switch, and a one-residue E66Q knockout (egfr_phsw_32) shipped beside its parent (egfr_phsw_23)
  as a bench falsification test.
- Switch selectivity up to d_occ 0.978 (egfr_phsw_27); OFF-state occupancy down to occ_HIE 0.00 (inv01).

## Method novelty

Three things here have no published single-chain precedent, and this track rewards method novelty:

1. **Two de-novo-novel pH-switch mechanisms** in one panel — the standard binder-acid-to-antigen-His433 salt
   bridge (proven at the antibody level by Liu 2022, realized here de novo) and an **inverted** switch that
   puts the titrating histidine on the binder, keyed to the conserved antigen Asp460.
2. **An MD-gated two-state selectivity test** that distinguishes a true switch from a constitutive always-on
   bridge — the step that a static pKa score cannot do and that most pH-design pipelines omit.
3. **A chemistry-fair frontier map** (GFN2-xTB + APBS + charge-scaling + novelty oracle) establishing that
   His salt bridges are the achievable pH-switch frontier on this epitope. Full account below.

## Target and mechanism

The target is the domain III epitope of human EGFR, the cetuximab and panitumumab surface. We use the EGFR
chain (chain A) of PDB 6ARU, a cetuximab-Fab bound EGFR ectodomain structure (3.20 Angstrom), and design
against its domain III. The goal is a binder that tightens at tumor pH (~6.5) relative to blood (7.4) and
holds on both human and mouse EGFR. Numbering is UniProt P00533 precursor throughout (the switch histidine
is His433, = His409 mature; Liu 2022 uses the same convention).

Both mechanisms exploit a histidine that titrates across the tumor/blood window (imidazole pKa ~6.3, mostly
protonated at tumor pH 6.5 and mostly neutral at blood pH 7.4), but they place the titrant differently.

**Mechanism 1 — binder acid to the antigen histidine His433 (18 designs: 17 active switches + 1 knockout
control).** Each of the 17 active designs places an aspartate or glutamate so it salt-bridges His433 only
when His433 is protonated; the interface gains a salt bridge at low pH. Liu et al. 2022 (PMID 36458200;
[DOI](https://doi.org/10.1016/j.omto.2022.11.001)) proved this at the antibody level: their anti-EGFR
antibody G532 carries an engineered LCDR1 Glu32 that salt-bridges EGFR His433, cross-reacts on human and
mouse, binds ~8-fold more strongly at pH 6.5 than 7.4 by ELISA (EC50 ratio 8.08) and is tumor-selective in
mice; the enabling Y32E carboxylate raised pH-dependency ~13-fold. We realize this de novo in a single chain,
which has no published experimental precedent and is the main risk.

**Mechanism 2 — binder histidine to the conserved antigen acid Asp460 (2 designs, egfr_phsw_inv01/inv02).**
The switch is inverted: the binder presents a histidine that salt-bridges the conserved antigen Asp460 only
when the binder histidine is protonated at tumor pH. This moves the titrating group onto the binder, where
its local environment, and therefore its pKa, is set by design rather than inherited from a buried antigen
histidine. Asp460 is identical in human and mouse, so the mechanism is cross-reactive by construction. A
first attempt against a flat conserved-Asp patch did not bind (0 of 384 folded); re-docking the inverted
switch onto the same proven-bindable His433 face with a two-hotspot design (His433 for the footprint, Asp460
for the switch) recovered both binding and a clean switch. This binder-side inverted pH switch also has no
published de novo precedent, and it is a second, mechanistically independent route to the same objective.

His433 and its neighbors are identical across the two species. We steered toward the six conserved domain III
anchors (precursor positions 408/432/433/436/464/489, identical in human and mouse) and away from the six
cetuximab-footprint positions that differ between species (377 R→K, 442 S→G, 467 K→R, 491 I→M, 492 S→N,
497 N→K). These come from our human-versus-mouse alignment (human P00533 vs mouse ortholog Q01279, 88.73%
whole-ECD identity over residues 25–645), recorded in `data/epitope/epitope_map.json`. The two
Mechanism-2 anchors, His370 (H/H) and Asp460 (D/D), are conserved by direct P00533-vs-Q01279 comparison
(no indel in the local alignment) and are recorded alongside His433 in that file. The His433 11-residue
context `GRTKQHGQFSL` and the His370 context `ISGDLHILPVA` are byte-identical between the species.

## Generation

Backbones came from RFdiffusion by partial diffusion of binders that had already passed folding (noise
partial_T=20 for Mechanism 1; partial_T 12/18/24 for the inverted lineage, whose backbones additionally came
from an all-atom RFdiffusion3 two-hotspot motif placing the binder histidine against Asp460). Partial
diffusion is the route that clears the structural novelty check: de novo runs land on common idealized folds,
while partial diffusion of an already-novel binder keeps the unusual topology. Sequences came from soluble
ProteinMPNN with cysteine disallowed and the switch residue locked (two sequences per backbone), each folded
on itself with AF2 initial-guess against a full-atom domain III target (the target must be full-atom or the
initial-guess template collapses). The generation gate kept interface PAE below ~11 for newly generated
designs; one carried-over anchor, egfr_phsw_04, sits just above at 11.4 and is kept on its strong orthogonal
binding signals (Chai 0.78, Boltz 0.64, PRODIGY −9.5) and its four-leg switch.

## Novelty

Adaptyv scores novelty by sequence identity (MMseqs2 vs SwissProt and PDB) and by structure (ESMFold →
Foldseek/TM-align vs PDB and AFDB, domain-segmented), requiring ≥3 of 4. The binding set is far below the
sequence thresholds, so the structure axis is the gate. **All 20 designs passed the portal novelty check,
user-run on 2026-10-06, including both inverted-mechanism designs** — the portal is the authority and is what
we report. We also rank novelty locally with a calibrated ESMFold2+TM-align oracle (safe-harbor TM < 0.80),
but it is only a ranker and cannot resolve the boundary where these designs sit: egfr_phsw_22 (local TM 0.89)
and egfr_phsw_27 (0.83) exceed the local safe harbor yet are portal-confirmed novel. The inverted mechanism
adds method-level novelty on top of structural novelty.

## pH-switch validation

Folding models (AF2, AF3, Boltz, Chai, ESMFold) are blind to pH, so the switch is a separate physics layer
validated on the exact submitted sequence by largely independent methods. We trust agreement in direction
(sign) across methods, not any single number, because no pKa method is reliable in absolute terms for an
engineered histidine.

- **PROPKA Wyman proton linkage:** ddG_pH = RT ln10 × the integral of the proton-count difference from pH 7.4
  to 6.5, from PROPKA pKa values of complex, target, and binder. Positive = more protons at the epitope at
  6.5, i.e. tighter at tumor pH.
- **pKAI:** an independent machine-learned pKa predictor on complex, target, and binder; we score the His pKa
  shift between bound and free (positive favors protonation on binding). pKAI and PROPKA share training data,
  so we treat them as corroborating, not fully independent, pKa legs.
- **Two-state MD occupancy — the selectivity gate:** implicit-solvent (GBn2) explicit-protonation MD (OpenMM,
  3 replicas, 200 ps equilibration, 2 ns production each). We run the complex with the sensor histidine
  protonated (HIP, tumor pH) and neutral (HIE, blood pH) and measure the fraction of the trajectory the bridge
  is formed: occ_HIP and occ_HIE. A positive static score only fixes direction; it does not prove selectivity,
  because a constitutive bridge that sits on the sensor whether or not it is charged scores positive and even
  passes the isostere control yet binds at both pH. A buried His–carboxylate bridge can push the histidine pKa
  well above 7 (in T4 lysozyme the buried His31–Asp70 bridge puts His31 at pKa 9.1), leaving no OFF state, and
  a single snapshot cannot tell that from a true switch. Only the two-state occupancy can. The gap
  d_occ = occ_HIP − occ_HIE separates a real switch (high occ_HIP, low occ_HIE) from a constitutive one. Most
  designs that cleared the static legs were constitutive on MD and were dropped; this is the step that
  separates the selected designs from the rest. Occupancies are estimates on limited sampling (2 ns/replica),
  carried with their replica spread and read as a ranked selectivity signal, not a converged free energy.
- **MM-GBSA two-state binding free energy:** single-trajectory MM-GBSA (Amber ff14SB, GBn2) on the same two
  states, dd_bind = dG_bind(HIP) − dG_bind(HIE). Negative = tighter with the sensor protonated. Read by sign,
  not magnitude: most systematic error cancels between two protonation states of the same complex. This is the
  only leg that scores the functional binding readout directly, independent of the pH-blind co-folders and of
  PROPKA.
- **Isostere ablation:** mutate the sensor-contacting acid to its neutral isostere (Asp→Asn, Glu→Gln), hold
  the rest fixed, and confirm the linkage collapses, ruling out the binder's bulk charge. Of the twelve Tier-1
  designs, three (egfr_phsw_22/23/29) show a clean collapse; egfr_phsw_27 is inconclusive; and the eight
  locked designs have a single-structure His433 linkage near zero, which makes the isostere control
  uninformative for them (not a refutation) — they rest on the MD occupancy gate.

The twelve Tier-1 designs are confirmed in sign by all four of PROPKA, pKAI, MD occupancy, and MM-GBSA.
Designs where PROPKA and pKAI disagreed on the His433 shift direction were dropped before selection.

### Engine positive control

Before trusting the static engine on an engineered histidine, we checked that it recovers a natural
low-pH-ON histidine switch: the rat FcRn:Fc complex (PDB 1FRT). Through the exact production legs, PROPKA
gives the correct low-pH-ON sign (ddG_pH +0.27 kcal/mol) and pKAI puts the shift on the two switch
histidines His310/His435 (mean dpKa +2.66) while staying flat on the four bystander histidines (max 0.19), a
negative control inside the positive control. Both legs pass, so the engine is not silently non-titrating.
1FRT is 4.5 Angstrom, so this tests recovery of the correct sign, not fine numeric accuracy.

### Knockout control

egfr_phsw_32 is egfr_phsw_23 with a single point mutation — the His433-salt-bridging Glu66 changed to Gln
(E66Q, position 66; the two submitted sequences differ at exactly that one residue, verified byte-for-byte).
The isostere analysis of the parent predicts this removes the switch while leaving the fold and the rest of
the interface intact, so egfr_phsw_32 should still bind but lose the pH preference. Shipping it beside its
parent is a built-in bench falsification test: if the parent is pH-selective and the mutant is not, the
designed bridge is the switch. Its predicted loss of switch is inferred from the parent's verified isostere
collapse and the E66Q chemistry, not computed by MD on the mutant. It remains de novo (a one-residue variant
of a design that passed the novelty gate).

## Cross-reactivity

For each binder we build a mouse domain III by substituting the human-versus-mouse differences, repack, and
refold with AF2 initial-guess alongside a human control with the identical repack; the mouse-minus-human
interface-PAE delta isolates the species substitutions from folding noise (cross-reactive within delta < 3.0,
preferred < 1.5). This is backed by three co-folder-independent signals:
a deterministic contact-conservation severity over the divergent contacts in the 4 Angstrom epitope core
(core severity 0 for the panel — all epitope contacts land on conserved positions), a SurfDiff structure
comparison of human-vs-mouse domain III, and matched Chai mouse/human iPTM. A single co-folder's mouse iPTM
is folder-wide pessimism, not a species signal: Spearman(Chai-mouse, Chai-human) is 0.869 with mean absolute
difference 0.053, so we never gate on it alone. The Mechanism-1 designs sit on the conserved His433 face; the
inverted designs are cross-reactive by construction (Asp460 is conserved; inv02 additionally has Chai-mouse
0.801, inv01 rests on the Asp460 conservation argument).

## Affinity and developability

Binding is read per predictor — AF2 initial-guess interface PAE, Chai and Boltz co-fold iPTM, and PRODIGY dG
where available — rather than collapsed into one number, because the folders disagree (Chai is the weakest
AF3-class folder here and reads pessimistically). On the three low-Chai Tier-1 designs the other predictors
agree the design binds (egfr_phsw_22 Boltz 0.85 / af2ig PAE 5.39; egfr_phsw_02 Boltz 0.86 / PRODIGY −11.4;
egfr_phsw_27 Boltz 0.59 / af2ig 6.73). Affinity is deliberately mid-range: higher affinity tends to come from
pH-independent contacts that wash out the selectivity, which is the top objective. Every design has zero free
cysteines (ProteinMPNN with Cys disallowed, a hard filter), is single chain, and is 56–108 residues.

## Design-space survey: why this mechanism

Before freezing on His433 salt bridges we surveyed whether a far-better switch exists on this epitope, using a
chemistry-fair scoring stack rather than the fixed-charge MD occupancy alone (which is biased: it under-counts
cation-pi and over-counts salt bridges). The stack is a GFN2-xTB interaction energy with self-consistent
charges (fair across chemistries), an APBS continuum net energy that accounts for desolvation, an
electronic-continuum-correction charge-scaling test (0.8 scaling) for ionic-strength robustness, and the
novelty oracle. Each alternative we could motivate was tested and rejected at a concrete wall:

- **Cation-pi** (binder-His⁺ on a conserved neutral aromatic): charge-scaling-robust, but on the fair QM a
  single cation-pi is comparable-to-weaker than a single salt bridge; the only novelty-passing designs are
  single (weaker), and the stronger dual configuration fails the novelty gate.
- **Multi-histidine neutral-acceptor collective:** a contact strong enough to overcome the buried-histidine
  desolvation penalty must be ionic; robust neutral contacts are too weak once buried (APBS net refutes).
- **Buried salt bridge:** strong but generically constitutive — the neutral-His off-state interaction energy
  is −11.6 to −12.0 kcal/mol across four buried designs (a neutral histidine still H-bonds a desolvated
  carboxylate), so there is no OFF state. A chemistry property, not a design flaw.
- **Repulsion-gated carboxylate dyad:** sound physics (the off-state is actively repulsive), but the
  make-or-break buried-glutamate pKa is unpredictable — PROPKA and pKAI disagree by 2–3 pH units per
  candidate, wider than the whole 6.5–7.4 window.
- **Allosteric hinged-lid swing:** primarily not buildable with the available generators (RFdiffusion +
  ProteinMPNN produced monolithic folds, 1 of 512 a binder, no hinge), and only partially verifiable (the
  pH-aware MD handle reports differential docked-state stability, not direct lid release).

The result is a three-way trade-off between desolvation, strength, and charge-scaling robustness: the only
residue whose pKa sits naturally in the 6.5–7.4 window and is predictable is histidine, so a His salt bridge
— forward to His433 or inverted from a binder histidine to Asp460 — is the achievable frontier on this
epitope. The panel sits on that frontier by design; we report it as a bounding result that justifies the
mechanism choice.

## Selection and panel composition

Twenty single-chain designs in four tiers (per-design numbers in `DESIGN_TABLE.md`):

- **Tier 1, twelve four-leg-confirmed Mechanism-1 designs** (02/04/06/08/09/12/16/19/22/23/27/29): all four pH
  legs correct in sign, binding by multi-predictor read, 0-Cys, portal-novel. Ranked by MD selectivity
  (d_occ 0.155–0.978) then binding, deliberately spanning the switch/binding trade-off. egfr_phsw_02
  (d_occ 0.155) and egfr_phsw_16 (0.195) sit at/below the pipeline's 0.2 MD-selectivity floor and are carried
  as portal-locked anchors (strong on binding/novelty), not for switch strength.
- **Tier 2, five portal-novel refills** (33/34/36/38/39): cysteine-cleaned children of Tier-1 lineages,
  validated on two-state MD occupancy + binding consensus + cross-react mechanism (not the full static
  four-leg set). Added to fill the panel to 18 after seven earlier candidates failed the portal novelty check.
- **Tier 3, one knockout control** (32): the E66Q matched pair with egfr_phsw_23 described above.
- **Tier 4, two inverted-mechanism designs** (inv01/inv02): clean two-state MD switch (d_occ 0.515 and 0.87,
  occ_HIE 0.00 and 0.13), 3-predictor binding, cross-reactive by construction.

The strongest risk hedge is the **two mechanistically independent switches**: a failure mode specific to the
buried-His433 route (e.g. a mispredicted antigen-His pKa) leaves the inverted binder-side switch untouched,
and vice versa. Lineage diversity is secondary and uneven — the per-lineage cap of two applies within Tier 1;
panel-wide, lineage 06 contributes four members (06, 23, refill 34, and the knockout 32, which is 23's E66Q),
lineage 16 three (16, 33, 39), and six lineages are singletons.

## Caveats

- Folding models do not see pH; all pH behavior rests on the pKa and MD layer, and MD occupancies are
  estimates on limited sampling, read as a ranked signal.
- No single pKa value is trustworthy for a buried histidine; the evidence is cross-method sign agreement plus
  the MD geometry and the isostere control. pKAI and PROPKA share training data (corroborating, not fully
  independent); the MD occupancy and MM-GBSA legs are the independent ones.
- The Tier-2 refills and Tier-4 inverted designs carry the MD and binding legs, not the full static four-leg
  set that the Tier-1 twelve carry (stated per design in the table).
- The knockout's loss of selectivity is predicted from the parent's isostere analysis and the E66Q chemistry,
  not measured by MD on the mutant.
- The inverted-mechanism salt bridge, like any surface salt bridge, weakens under the charge-scaling test; it
  ships as a mechanistically distinct hedge, not as a more-robust switch than Mechanism 1.
- The local novelty proxy is a ranker, not the gate; the portal's domain-segmented check is the authority and
  all 20 passed it.

## Reproduce

- Generation (Mechanism 1): `scripts/wave18_novel_wave.sh runs/pd_diverse_seeds.txt 400 20 pdnov "0 1 2 3 4 5 6 7"`
- Inverted mechanism (two-hotspot RFd3 His433+Asp460): `scripts/run_moonshot_rfd3.sh`, `scripts/moonshot_mpnn_af2ig.sh`, `scripts/run_pd_inv02.sh`
- MD selectivity gate (two protonation states): `scripts/run_md_campaign.py --pdb-list <list> --out <dir> --replicas 3` (inverted: `scripts/md_pH_switch.py --binder-his --switch-his-pos <pos>`)
- MM-GBSA two-state dd_bind: `scripts/md_pH_bind.py --pdb <design>__best.pdb --out <dir>`
- pKAI / PROPKA Wyman / isostere / cross-react / FcRn: `scripts/obj1_pkai.py`, `scripts/obj1_switch.py`, `scripts/wave16k_isostere_ablation.py`, `scripts/run_xreact.sh`, `scripts/obj1_fcrn_control.py`
- Chemistry-fair frontier stack: `scripts/qm_switch_fair.py --solvent {water,chcl3,vacuum}` (GFN2-xTB) + APBS net linkage + novelty oracle
- Metrics + validate: `submission/metrics.csv`; `uv run python scripts/validate_submission.py submission/designs.csv --track 3`

## Citations

- Liu et al. 2022, PMID 36458200, [DOI](https://doi.org/10.1016/j.omto.2022.11.001): cross-reactive,
  pH-dependent anti-EGFR domain III antibody via the His433 salt bridge (the Mechanism-1 precedent; His370 and
  His433 both in domain III and conserved human/mouse).
- Human-vs-mouse EGFR alignment (P00533 vs Q01279, 88.73% whole-ECD identity), `data/epitope/epitope_map.json`.
- Burmeister, Huber, Bjorkman 1994, PMID 7969498: rat FcRn:Fc (PDB 1FRT), the engine positive control.
- RFdiffusion (Watson et al. 2023) and RFdiffusion3; ProteinMPNN (Dauparas et al. 2022); AlphaFold2 (Jumper et
  al. 2021) and AF2 initial-guess (Bennett et al. 2023); AlphaFold3 (Abramson et al. 2024); ESMFold (Lin et
  al. 2023); Chai-1 (Chai Discovery 2024); Boltz-2 (Passaro et al. 2025); PROPKA (Olsson et al. 2011); pKAI
  (Reis et al. 2022); Foldseek (van Kempen et al. 2024); OpenMM (Eastman et al. 2017); Amber ff14SB (Maier et
  al. 2015) with GBn2 (Nguyen et al. 2013); PRODIGY (Xue et al. 2016); GFN2-xTB (Bannwarth et al. 2019); APBS
  (Jurrus et al. 2018).
