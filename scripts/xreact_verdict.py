#!/usr/bin/env python3
"""WAVE-15 cross-reactivity verdict: fuse the DIRECT af2ig mouse-vs-human co-fold with the MECHANISM
severity signal (xreact_mechanism.py) into a 2-signal per-design verdict (objective 2).

Protocol (run_xreact.sh): mutate target chain B to mouse DIII at all divergent positions + repack
(mouse complex); repack-only control (human); af2ig-fold the binder against BOTH. The CONTROLLED
cross-reactivity measure is delta = pae_interaction(mouse) - pae_interaction(human_repack): both get
the identical repack, so delta isolates the human->mouse substitutions (absolute PAE is inflated ~1
unit by the repack and is a separate BINDING-STRENGTH axis, not cross-reactivity).

    delta < 1.5  & severity_core == 0 -> cross-reactive (2-signal)
    delta < 1.5                       -> cross-reactive (fold; mech-watch)
    1.5 <= delta < 3.0                -> cross-reactive (mild degrade)
    delta >= 3.0                      -> species-specific RISK
    (both species' repacked interface-PAE >= 9.5 adds a "[weak binder]" affinity note -- symmetric,
     orthogonal to cross-reactivity, not a gate)

Usage: uv run python scripts/xreact_verdict.py
       [--mech results/wave15_xreact_mechanism.csv]
       [--mouse-dir runs/af2ig_xreact_fin20_mouse] [--human-dir runs/af2ig_xreact_fin20_human]
       [--out results/wave15_xreact_verdict.csv]
"""
from __future__ import annotations
import argparse, csv, glob, re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def parse_sc(run_dir: str):
    """pae_interaction per design. Reads the merged top-level out.sc FIRST, then per-gpu shards, keeping
    the LAST-SEEN value (shards win) so a stale leftover top-level out.sc cannot override fresh shards."""
    out = {}
    paths = glob.glob(str(Path(run_dir) / 'out.sc')) + sorted(glob.glob(str(Path(run_dir) / 'gpu*' / 'out.sc')))
    for path in paths:
        hdr = None
        with open(path) as fh:
            for l in fh:
                p = l.split()
                if not p or p[0] != 'SCORE:':
                    continue
                if 'description' in p:
                    hdr = p[1:]
                    continue
                if hdr is None or len(p[1:]) != len(hdr):   # skip malformed / mis-tokenised rows
                    continue
                d = dict(zip(hdr, p[1:]))
                m = re.search(r'(egfr_phsw_\d+)', d.get('description', ''))
                if not m:
                    continue
                try:
                    pae = float(d['pae_interaction'])
                except (KeyError, ValueError):
                    continue
                out[m.group(1)] = pae
    return out


def _sev(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 99          # missing/unparseable severity -> treat as NOT mechanism-confirmed (conservative)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mech', default=str(REPO / 'results/wave15_xreact_mechanism.csv'))
    ap.add_argument('--mouse-dir', default=str(REPO / 'runs/af2ig_xreact_fin20_mouse'))
    ap.add_argument('--human-dir', default=str(REPO / 'runs/af2ig_xreact_fin20_human'))
    ap.add_argument('--out', default=str(REPO / 'results/wave15_xreact_verdict.csv'))
    a = ap.parse_args()

    mouse, human = parse_sc(a.mouse_dir), parse_sc(a.human_dir)
    with open(a.mech) as fh:
        mech = {r['name']: r for r in csv.DictReader(fh)}

    rows = []
    for name in sorted(mech):
        mp, hp = mouse.get(name), human.get(name)
        sev_core, sev_sum = _sev(mech[name].get('xreact_severity_core4A')), _sev(mech[name].get('xreact_severity'))
        mech_ok = sev_core == 0
        delta = (mp - hp) if (mp is not None and hp is not None) else None   # raw, for thresholds
        weak = (mp is not None and hp is not None and min(mp, hp) >= 9.5)
        if delta is None:
            verdict = 'fold pending'
        elif delta < 1.5 and mech_ok:
            verdict = 'cross-reactive (2-signal)'
        elif delta < 1.5:
            verdict = 'cross-reactive (fold; mech-watch)'
        elif delta < 3.0:
            verdict = 'cross-reactive (mild degrade)'
        else:
            verdict = 'species-specific RISK'
        if weak:
            verdict += ' [weak binder]'
        rows.append(dict(name=name,
                         mouse_pae=(round(mp, 2) if mp is not None else None),
                         human_pae=(round(hp, 2) if hp is not None else None),
                         delta_pae=(round(delta, 2) if delta is not None else None),
                         sev_core=(sev_core if sev_core != 99 else ''),
                         sev_sum=(sev_sum if sev_sum != 99 else ''),
                         chai_mouse=mech[name].get('chai_mouse', ''), chai_human=mech[name].get('chai_human', ''),
                         verdict=verdict))

    def fx(v, f='%6.2f'):
        return (f % v) if isinstance(v, float) else ('%6s' % ('-' if v in (None, '') else v))

    print('%-13s %8s %8s %8s %5s %5s %8s %8s   %s' % (
        'name', 'mousePAE', 'humanPAE', 'deltaPAE', 'sevC', 'sevS', 'chaiMo', 'chaiHu', 'verdict'))
    for r in sorted(rows, key=lambda r: (r['mouse_pae'] if r['mouse_pae'] is not None else 99)):
        print('%-13s %8s %8s %8s %5s %5s %8s %8s   %s' % (
            r['name'], fx(r['mouse_pae']), fx(r['human_pae']), fx(r['delta_pae']),
            str(r['sev_core']), str(r['sev_sum']), r['chai_mouse'], r['chai_human'], r['verdict']))
    ok = sum(1 for r in rows if 'cross-reactive' in r['verdict'])
    print(f"\naf2ig parsed mouse {len(mouse)} / human {len(human)} ; cross-reactive: {ok}/{len(rows)}")

    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['name', 'mouse_pae', 'human_pae', 'delta_pae',
            'sev_core', 'sev_sum', 'chai_mouse', 'chai_human', 'verdict'])
        w.writeheader(); w.writerows(rows)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
