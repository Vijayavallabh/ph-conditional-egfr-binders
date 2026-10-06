# Submission methodology

Twenty de novo single-chain binders to the human EGFR domain III epitope (the clinically validated
cetuximab/panitumumab surface), designed to bind preferentially at tumor pH (~6.5 vs blood 7.4) and to
cross-react on human and mouse EGFR. Objective priority: pH-selectivity, then cross-reactivity, then affinity;
novelty is a gate. All 20 passed the Adaptyv portal novelty check (≥3/4, user-run 2026-10-06), are 0-Cys
single chains of 56–108 aa.

**Design.** Backbones by RFdiffusion partial diffusion of folding-validated binders (and an all-atom
RFdiffusion3 two-hotspot motif for the inverted mechanism); sequences by soluble ProteinMPNN (cysteine
disallowed, switch residue locked); self-consistency by AF2 initial-guess against a full-atom domain III
target. Partial diffusion is what clears the structural novelty check.

**Two switch mechanisms.** (1) A binder Asp/Glu that salt-bridges EGFR His433 when His433 protonates at acidic
pH — the Liu-2022 G532 antibody mechanism (PMID 36458200), realized de novo single-chain (18 designs: 17
active switches plus one E66Q knockout control). (2) An inverted switch in which a binder histidine
salt-bridges the human/mouse-conserved Asp460, putting the titrating residue on the binder where its pKa is
designable (2 designs). His433, His370 and Asp460 are identical in human and mouse.

**pH-selectivity validation (on the exact submitted sequence).** For the 12 Tier-1 designs, a four-leg
sign-consensus: PROPKA Wyman proton-linkage, the pKAI machine-learned pKa shift (corroborating PROPKA; shared
training data), explicit-protonation two-state MD occupancy (occ_HIP vs occ_HIE; the gap d_occ is the
selectivity gate), and two-state MM-GBSA binding free energy (the independent functional leg). The five
refills and two inverted designs carry the two-state MD occupancy and binding legs. The two-state MD occupancy
is decisive: a positive static pKa score only fixes the switch direction, while a constitutive always-on
bridge scores positive yet binds at both pH; only MD occupancy separates them. An FcRn:Fc (PDB 1FRT) positive
control confirms the engine recovers a known low-pH-ON histidine switch. We read agreement in sign across
methods, not any single pKa value.

**Cross-reactivity.** A four-signal consensus independent of any single co-folder: deterministic
contact-conservation severity over the 4 Å epitope core, matched mouse−human AF2-ig interface-PAE delta,
SurfDiff structure comparison, and Chai mouse/human iPTM. All 20 are cross-reactive; the inverted designs are
cross-reactive by construction (Asp460 conserved).

**Controls and developability.** A matched knockout pair (egfr_phsw_23 and its single-residue E66Q mutant
egfr_phsw_32) is a built-in bench falsification test (the mutant's predicted loss of switch is inferred from
the parent's isostere collapse, not computed on the mutant). Binding is read per predictor (AF2-ig interface
PAE, Chai/Boltz iPTM, PRODIGY) rather than one blended number, since the folders disagree and Chai reads
pessimistically here.

**Method novelty.** The mechanism choice is backed by a survey of the switch design space with a
chemistry-fair scoring stack (GFN2-xTB interaction energies, APBS continuum net electrostatics,
electronic-continuum-correction charge-scaling, and the novelty oracle). Every alternative (cation-π,
multi-histidine neutral collective, buried salt bridge, repulsion-gated carboxylate dyad, allosteric hinged
lid) meets a concrete physical wall, bounding the problem: His-based salt bridges are the achievable
pH-switch frontier on this epitope.
