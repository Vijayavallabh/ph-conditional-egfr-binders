#!/usr/bin/env python3
"""WAVE-15 mechanism-based cross-reactivity signal (objective 2: human/mouse cross-reactivity).

Co-folders (Chai/Boltz/AF3) are noisy and species-blind in practice; this adds an orthogonal,
deterministic, leakage-free signal. For each exact-sequence binder:target complex it finds the
target (chain B) residues the binder (chain A) contacts (<5 A any-atom), maps them to crop ordinal
(normalised = chainB resnum - (min chainB resnum - 1); == resnum for the 6ARU DIII crop which is
numbered 1..191), and scores each contacted position that DIFFERS between human (P00533) and mouse
(Q01279) EGFR DIII by substitution severity:

    charge reversal (+/-)            = 3   (salt bridge broken in mouse)
    charge gain/loss (charged<->0)   = 2
    class change (hyd/aro/Pro)       = 1
    conservative (R<->K, I<->M, ...) = 0   (tolerated -> mouse binding preserved)

A low severity within the tight 4 A core => the binder's epitope substitutions are tolerated in mouse
=> cross-reactive by mechanism, independent of any co-folder. Divergent set + residues come from
data/targets/mouse_domainIII_map.json (make_mouse_target.py).

Usage: uv run python scripts/xreact_mechanism.py
       [--indir runs/md_bind_gplus/in] [--map data/targets/mouse_domainIII_map.json]
       [--manifest results/submission_manifest_v2_wave14.csv] [--out results/wave15_xreact_mechanism.csv]
"""
from __future__ import annotations
import argparse, csv, glob, json
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
CUT, CORE = 5.0, 4.0
CHG = {'D': -1, 'E': -1, 'K': 1, 'R': 1}          # His neutral at pH 7.4 for flip scoring
HYD, ARO = set('AVLIMFWC'), set('FWY')


def severity(h: str, mm: str) -> int:
    if h == mm:
        return 0
    ch, cm = CHG.get(h, 0), CHG.get(mm, 0)
    if (ch > 0 and cm < 0) or (ch < 0 and cm > 0):
        return 3
    if (ch != 0) != (cm != 0):
        return 2
    if 'P' in (h, mm) or (h in ARO) != (mm in ARO) or (h in HYD) != (mm in HYD):
        return 1
    return 0


def analyze(pdb, POS, DIV, CONS, his433):
    A, B = [], {}
    with open(pdb) as fh:
        for l in fh:
            if l[:4] != 'ATOM':
                continue
            if l[16] not in (' ', 'A'):          # keep only the primary altLoc
                continue
            try:
                x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54]); num = int(l[22:26])
            except ValueError:
                continue
            if l[21] == 'A':
                A.append((x, y, z))
            elif l[21] == 'B':
                B.setdefault(num, []).append((x, y, z))
    if not A or not B:
        return None
    A = np.array(A)
    offset = min(B) - 1                            # normalise chain-B resnum -> crop ordinal (1..191)
    dist = {num - offset: float(np.sqrt(((np.array(xs)[:, None, :] - A[None, :, :]) ** 2).sum(-1)).min())
            for num, xs in B.items()}
    con = {ordv for ordv, d in dist.items() if d < CUT}
    div_c = sorted(con & DIV)
    sev_sum = sum(severity(POS[o]['human'], POS[o]['mouse']) for o in div_c)
    sev_core = sum(severity(POS[o]['human'], POS[o]['mouse']) for o in div_c if dist[o] < CORE)
    detail = ';'.join(f"{o}:{POS[o]['human']}>{POS[o]['mouse']}({severity(POS[o]['human'], POS[o]['mouse'])},{dist[o]:.1f}A)"
                      for o in div_c)
    return {
        'n_contacts': len(con), 'n_divergent_contacts': len(div_c),
        'n_conserved_anchor_contacts': len(con & CONS), 'his433_contact': his433 in con,
        'xreact_severity': sev_sum, 'xreact_severity_core4A': sev_core,
        'min_dist_div_in_epitope_A': round(min((dist[o] for o in div_c), default=99.0), 2),
        'divergent_detail': detail,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', default=str(REPO / 'runs/md_bind_gplus/in'))
    ap.add_argument('--map', default=str(REPO / 'data/targets/mouse_domainIII_map.json'))
    ap.add_argument('--manifest', default=str(REPO / 'results/submission_manifest_v2_wave14.csv'))
    ap.add_argument('--out', default=str(REPO / 'results/wave15_xreact_mechanism.csv'))
    a = ap.parse_args()

    with open(a.map) as fh:
        m = json.load(fh)
    POS = {int(p['ordinal']): p for p in m['positions']}
    DIV = {o for o, p in POS.items() if not p.get('same', True)}
    CONS = {int(k) for k in m['conserved_anchor_check']}
    # His433 == author 6ARU 409 (and residue 'H'); fail loud if the map is renumbered (no silent default).
    his433 = next((o for o, p in POS.items() if p['human'] == 'H' and int(p.get('author_6aru', 0)) == 409), None)
    if his433 is None:
        raise SystemExit("His433 (human 'H', author_6aru 409) not found in mouse map — check the map/numbering")

    chai = {}
    if Path(a.manifest).exists():
        with open(a.manifest) as fh:
            for r in csv.DictReader(fh):
                n = r.get('name')
                if n:
                    chai[n] = (r.get('chai_mouse_iptm', ''), r.get('chai_iptm_human', ''))

    rows = []
    for pdb in sorted(glob.glob(str(Path(a.indir) / '*egfr_phsw_*.pdb'))):
        r = analyze(pdb, POS, DIV, CONS, his433)
        if r is None:
            continue
        r['name'] = Path(pdb).stem
        r['chai_mouse'], r['chai_human'] = chai.get(r['name'], ('', ''))
        rows.append(r)

    rows.sort(key=lambda r: (r['xreact_severity_core4A'], r['xreact_severity'], r['n_divergent_contacts']))
    print('%-13s %4s %4s %7s %9s %8s %8s   %s' % (
        'name', 'ndiv', 'nanc', 'sevSum', 'sevCore4A', 'chaiMo', 'chaiHu', 'divergent_detail'))
    for r in rows:
        print('%-13s %4d %4d %7d %9d %8s %8s   %s' % (
            r['name'], r['n_divergent_contacts'], r['n_conserved_anchor_contacts'],
            r['xreact_severity'], r['xreact_severity_core4A'], r['chai_mouse'], r['chai_human'],
            r['divergent_detail']))
    clean = [r['name'] for r in rows if r['xreact_severity_core4A'] == 0]
    print(f"\ncore-4A severity == 0 (all tight divergent contacts tolerated): {len(clean)}/{len(rows)}")

    cols = ['name', 'n_contacts', 'n_divergent_contacts', 'n_conserved_anchor_contacts', 'his433_contact',
            'xreact_severity', 'xreact_severity_core4A', 'min_dist_div_in_epitope_A',
            'chai_mouse', 'chai_human', 'divergent_detail']
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda r: r['name']):
            w.writerow({k: r[k] for k in cols})
    print('wrote', a.out)


if __name__ == '__main__':
    main()
