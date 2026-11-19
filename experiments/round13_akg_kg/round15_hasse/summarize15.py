"""Collect Round 15 results into round15.json and the two tables of ROUND15.md (summary15.md)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import _paths  # noqa: E402  (run from anywhere)
import json, glob
NG = json.load(open('../round14b_euler/noghosts.json'))['A']
key14 = dict(ZETA='ZETA_f45.0', EPS='EPS_f45.0', DH='DH_f125.0', BEUR2='BEUR2_f45.0', FTH='FTH0.2_f45.0')
out = dict(blocks={}, grams={}, ff=json.load(open('ff_hasse.json')),
           checks=dict(B0=json.load(open('check_core.json')), shift=json.load(open('check_shift.json')), quad=json.load(open('check_quad.json'))))
for f in sorted(glob.glob('blocks_*_l*.json')):
    d = json.load(open(f)); out['blocks'][f"{d['model']}_l{d['l']}"] = {k: v for k, v in d.items() if k != 'rows'}
for f in sorted(glob.glob('gram_*_l*.json')):
    d = json.load(open(f)); out['grams'][f"{d['model']}_{d['family']}"] = d
json.dump(out, open('round15.json', 'w'), indent=1)
MS = ['ZETA', 'FTAU', 'EPS', 'DH', 'BEUR2', 'FTH']
T1 = ['| model | λ_min(W_l) l=0.6 | l=0.6: violating shifts / scanned | first violating s (first prime p) | min μ | max ρ | l=0.3: violations, first s | l=0.3 max ρ |', '|---|---|---|---|---|---|---|---|']
for m in MS:
    b = out['blocks'][f'{m}_l0.6']; c = out['blocks'][f'{m}_l0.3']; pv = b['prime_shifts_violating']
    fv = f"{b['first_violation_s']:.2f} (p = {pv[0]})" if b['first_violation_s'] is not None else 'none'
    fc = f"{c['n_violations']}, {c['first_violation_s']:.2f}" if c['first_violation_s'] is not None else 'none'
    T1.append(f"| {m} | {b['W_lam'][0]:.3g} | {b['n_violations']}/{b['n_shifts']} | {fv} | {b['min_mu']:.2e} | {b['max_rho']:.4f} | {fc} | {c['max_rho']:.4f} |")
def first_ghost(d):
    r = next((r for r in d['rows'] if r['n_ghosts'] > 0), None); return f"{r['width']:.2f}" if r else '—'
T2 = ['| model | P2: width 8.2 λ_min / ghosts / onset width | P3: 8.3 | P5: 7.0 | INT (K=40, width 4.29) | PP (K=128, width 5.45) | 14B ghost onset L |', '|---|---|---|---|---|---|---|']
for m in MS:
    cells = []
    for fam in ['P2', 'P3', 'P5', 'INT', 'PP']:
        d = out['grams'][f'{m}_{fam}']; r = d['rows'][-1]
        cells.append(f"{r['lam_min']:.1e} / {r['n_ghosts']} / {first_ghost(d)}")
    on = next((r['L'] for r in NG.get(key14.get(m, ''), []) if r['n_ghosts'] > 0), None)
    T2.append(f"| {m} | " + ' | '.join(cells) + f" | {on if on is not None else ('none to 8' if m in key14 else '— (14F: 0 ghosts L=2,4)')} |")
s = '\n'.join(T1) + '\n\n' + '\n'.join(T2) + '\n'
open('summary15.md', 'w').write(s); print(s)
