"""Round 14E2: collect mp80/gap_*.json into round14e2.json and print the markdown tables used in ROUND14E2.md."""
import json, glob, os
HERE = os.path.dirname(os.path.abspath(__file__)); rows = []
for f in sorted(glob.glob(os.path.join(HERE, 'mp80', 'gap_*.json'))):
    for r in json.load(open(f)): r['control'] = f.endswith('_control.json'); rows.append(r)
json.dump(rows, open(os.path.join(HERE, 'round14e2.json'), 'w'), indent=1)
def g(x, k, d='–'): v = x.get(k); return d if v is None else v
print('| model | L | fmax (M) | TAU | eps | dim S | bilinear ||WS||/||W|| | lam_1(W) | restriction lam_min(Q^T W Q) | quotient lam_min (Schur W/W_SS) | interlace cap lam_{k+1}(W) | sin angle(S, bottom-k eigvecs) | ghosts W -> restr |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in sorted(rows, key=lambda r: (r['control'], r['model'] != 'ZETA', r['L'], r['F'], r['tau'] < 1e-20, -r['eps'])):
    if r['dimS'] == 0: continue
    print(f"| {r['model']}{' (control F-0.1·2^-s)' if r['control'] else ''} | {r['L']} | {r['F']} ({r['M']}) | {r['tau']:.0e} | {r['eps']:.0e} | {r['dimS']} | {g(r,'bilinear_resid_rel')} | {r['low_W'][0]} | {r['lam_min_restr']} | {g(r,'lam_min_quot')} | {r['interlace_upper']} | {g(r,'sin_max_angle_S_vs_bottom_k_eigvecs')} | {r['n_neg_W']} -> {g(r,'n_neg_restr')} |")
print()
for r in rows:
    if r['dimS'] == 0: print('empty S:', r['model'], r['L'], r['F'], r['tau'], r['eps'], 'rho top 1-rho:', r['rho_top'][:2])
