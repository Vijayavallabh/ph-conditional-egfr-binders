#!/usr/bin/env python
"""
WAVE 16b - POSITIVE CONTROL for the objective-1 pH-switch engine on a NATURAL His switch.

Validates the validator that selects/ranks the submitted 20: do our EXISTING obj1 legs correctly
recover a gold-standard, well-characterised bind-at-low-pH histidine switch (and stay flat on the
bystander residues that are NOT at the switch)?

System: rat FcRn:Fc complex (PDB 1FRT; Burmeister, Huber & Bjorkman 1994, Nature 372:379, PMID 7969498,
doi:10.1038/372379a0 -- the 1FRT deposition itself). The pH-dependent salt-bridge mechanism is resolved
at higher resolution in a DIFFERENT (engineered heterodimeric-Fc) structure, Martin et al. 2001
(Mol Cell 7:867, PMID 11336709) -- cited for mechanism only, NOT as 1FRT. Fc His310/His435
(chain C) protonate at acidic pH and salt-bridge FcRn acidic residues (chain A) => the complex binds
at pH<=6.5 and releases at 7.4 -- the SAME mechanism and SAME SIGN as our His433 design (the
His-on-one-partner / acid-on-the-other role swap is sign-neutral for the Wyman linkage).

Two legs, run in their own envs (reuses the obj1 ENGINES that score the 20; each leg notes exactly
what is and is NOT shared with production):
  --engine propka  (envs/pyrosetta; reuses obj1_switch's PROPKA run + _f Wyman kernel): proton-linkage
                   ddg_pH over pH 6.5->7.4, summed over ALL matched titratable pairs -- a TOTAL linkage,
                   NOT the His433/Glu496 sensor-restricted subset production ranks by (FcRn's switch His
                   differ), so this leg validates the linkage KERNEL + its SIGN, not the sensor-selection.
                   PASS = correct SIGN (>0 = binds tighter at low pH). PROPKA is a known-weak His predictor
                   (it credits carboxylate desolvation shifts), so only the SIGN is trusted here -- which
                   is exactly why the project also runs pKAI + MD (sign-consensus).
  --engine pkai    (envs/pkai; reuses obj1_pkai.run_pkai VERBATIM, the production pKA engine):
                   His-specific dpKa = pKa(complex)-pKa(Fc-alone)
                   for the switch His (310/435) vs the bystander His (268/285/429/433). PASS = switch His
                   dpKa>0 (FcRn stabilises protonated His = low-pH-ON) AND bystander His ~flat (built-in
                   negative control).

Implements the project's own "validation gates must fail loud / add a positive control" lesson and
guards against a silent non-titration bug (cf. the dead Rosetta pH_mode).

Run both:
  PATH=$PWD/envs/pyrosetta/bin:$PATH envs/pyrosetta/bin/python scripts/obj1_fcrn_control.py <1FRT.pdb> --engine propka
  envs/pkai/bin/python scripts/obj1_fcrn_control.py <1FRT.pdb> --engine pkai
"""
from __future__ import annotations
import os, sys, argparse, tempfile, shutil

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))

FCRN_CHAIN = "A"   # FcRn alpha-chain. beta2m (chain B) is the MHC-like light chain; the His310/His435
                   # pH-switch salt bridges are to the FcRn alpha2-domain acidics, so B is omitted from the
                   # 2-body cycle. Verified robust on 1FRT: including B leaves the switch dpKa essentially
                   # unchanged (His310 +4.71->+4.98, His435 +0.61->+0.54; bystanders flat both ways).
FC_CHAIN = "C"     # Fc (carries the His310/His435 switch)
SWITCH_HIS = (310, 435)
BYSTANDER_HIS = (268, 285, 429, 433)   # His433 sits near the CH2-CH3 FcRn footprint; empirically flat
                                       # here (|dpKa| 0.19) but it is the least-independent of the four.


def write_chains(pdb, chains, out):
    keep = set(chains)                 # hoisted once (was rebuilt per line)
    with open(out, "w") as o:
        for ln in open(pdb):
            if ln[:6] == "ATOM  " and ln[21] in keep:
                o.write(ln)
        o.write("END\n")


def run_propka(pdb):
    import obj1_switch as o1
    PKA_CLAMP, RT_LN10 = o1.PKA_CLAMP, o1.RT_LN10

    # KEEP IN SYNC with obj1_switch.ddg_ph_from_pdb: clamp() + the grid/trapz/contrib linkage below
    # mirror code nested (hence non-importable) in obj1_switch. PKA_CLAMP/RT_LN10 are imported (track
    # automatically); the integration kernel is copied, so a change there to the pH window (6.5..7.4),
    # grid step or temperature must be mirrored here or this control stops certifying the live scorer.
    def clamp(groups):
        out = {}
        for rn, n, c, p in groups:
            if rn in PKA_CLAMP:
                lo, hi = PKA_CLAMP[rn]; out[(rn, n, c)] = min(max(p, lo), hi)
        return out

    wd = tempfile.mkdtemp(prefix="fcrn_")
    try:
        write_chains(pdb, [FCRN_CHAIN, FC_CHAIN], os.path.join(wd, "cx.pdb"))
        write_chains(pdb, [FCRN_CHAIN], os.path.join(wd, "tg.pdb"))
        write_chains(pdb, [FC_CHAIN], os.path.join(wd, "bd.pdb"))
        try:
            gcx = clamp(o1._parse_pka(o1._run_propka("cx.pdb", wd)))
            gtg = clamp(o1._parse_pka(o1._run_propka("tg.pdb", wd)))
            gbd = clamp(o1._parse_pka(o1._run_propka("bd.pdb", wd)))
        except FileNotFoundError:
            sys.exit("[FATAL] propka3 not on PATH -- re-run with PATH=$PWD/envs/pyrosetta/bin:$PATH "
                     "(see the module docstring).")
        if not gcx or not gtg or not gbd:
            sys.exit("[FATAL] PROPKA produced no groups -> engine/propka broken.")
        free = {}; free.update(gtg); free.update(gbd)
        pairs = [(k, gcx[k], free[k]) for k in gcx if k in free]
        grid = [6.5 + 0.05 * k for k in range(19)]
        trapz = lambda v: sum((v[i] + v[i + 1]) / 2 * (grid[i + 1] - grid[i]) for i in range(len(v) - 1))
        contrib = lambda pc, pf: RT_LN10 * trapz([o1._f(pc, p) - o1._f(pf, p) for p in grid])
        rows, total = [], 0.0
        for k, pc, pf in pairs:
            gi = contrib(pc, pf); total += gi; rows.append((k, gi, pf, pc))
        rows.sort(key=lambda x: x[1], reverse=True)
        print(f"[PROPKA-Wyman] FcRn:Fc ddg_pH_total = {total:+.3f} kcal/mol "
              f"({'LOW-pH-ON (correct sign)' if total > 0 else 'WRONG sign'})")
        for (rn, n, c), gi, pf, pc in rows[:6]:
            if abs(gi) > 0.02:
                print(f"    driver {rn}{n}{c}  contrib {gi:+.3f}  pKa {pf:.1f}->{pc:.1f}")
        ok = total > 0.1
        print(f"PROPKA VERDICT: {'PASS (correct sign on a natural low-pH-ON His switch)' if ok else 'FAIL'}")
        print("  caveat: PROPKA credits carboxylate desolvation, not the His (documented His-weakness)")
        print("          -> only the SIGN is trusted; His attribution is the pKAI leg's job (--engine pkai).")
        return 0 if ok else 2
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def run_pkai(pdb):
    from obj1_pkai import run_pkai as pkai_pkas
    wd = tempfile.mkdtemp(prefix="fcrnpkai_")
    try:
        cx = os.path.join(wd, "cx.pdb"); fc = os.path.join(wd, "fc.pdb")
        write_chains(pdb, [FCRN_CHAIN, FC_CHAIN], cx)
        write_chains(pdb, [FC_CHAIN], fc)
        CX = pkai_pkas(cx, "pKAI"); FC = pkai_pkas(fc, "pKAI")
        print("[pKAI] Fc His pKa shift on FcRn binding (dpKa = pKa_complex - pKa_FcAlone; >0 = low-pH-ON):")
        sw, by = [], []
        for n in SWITCH_HIS + BYSTANDER_HIS:
            c, f = CX.get((FC_CHAIN, n)), FC.get((FC_CHAIN, n))
            if not (c and f):
                continue
            if c[0] != "HIS" or f[0] != "HIS":   # only bin genuine His (match obj1_pkai's guard) -- a
                print(f"    C{n}: not HIS (cx={c[0]} / free={f[0]}) -> skipped")  # non-His here = wrong chain
                continue
            d = c[1] - f[1]
            (sw if n in SWITCH_HIS else by).append(d)
            tag = "  <== switch His" if n in SWITCH_HIS else "  (bystander/null)"
            print(f"    C{n}: {c[0]} pKa_cx={c[1]:.2f} pKa_free={f[1]:.2f} dpKa={d:+.2f}{tag}")
        if not sw:   # fail loud: the switch His MUST be measured, else a wrong-chain/numbering PASS is vacuous
            sys.exit(f"[FATAL] no switch His (310/435) resolved as HIS in Fc chain '{FC_CHAIN}' "
                     "-> wrong chain or numbering; refusing a vacuous verdict.")
        if not by:   # the bystander NEGATIVE control must actually be exercised (no flat-by-default pass)
            sys.exit("[FATAL] no bystander His resolved -> the flat-bystander null was not tested.")
        sw_mean = sum(sw) / len(sw)
        by_absmax = max(abs(x) for x in by)
        print(f"  switch-His mean dpKa = {sw_mean:+.2f} ; bystander |dpKa| max = {by_absmax:.2f} "
              f"({len(sw)} switch / {len(by)} bystander His measured)")
        ok = sw_mean > 0.3 and by_absmax < 0.5
        print(f"pKAI VERDICT: {'PASS (His-specific leg recovers the switch; flat on bystanders)' if ok else 'FAIL'}")
        return 0 if ok else 2
    finally:
        shutil.rmtree(wd, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--engine", choices=["propka", "pkai"], default="propka")
    a = ap.parse_args()
    if not os.path.exists(a.pdb):
        sys.exit(f"[FATAL] FcRn:Fc PDB not found: {a.pdb}")
    print(f"FcRn:Fc obj1 positive control ({a.engine}) on {a.pdb} "
          f"(FcRn={FCRN_CHAIN}, Fc={FC_CHAIN}; switch His {SWITCH_HIS})")
    return run_propka(a.pdb) if a.engine == "propka" else run_pkai(a.pdb)


if __name__ == "__main__":
    sys.exit(main())
