"""Round 34 figures + control summary."""
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = json.load(open('r34_results.json')); C = json.load(open('r34_controls.json'))
summ = {}
for nm, m in C['models'].items():
    E = np.array(m['Ek']); z = np.linspace(-1, 1, 201); Ez = np.array([np.polyval(E[::-1], t) for t in z])
    summ[nm] = dict(E0=E[0], E1=E[1], E2=E[2], R=float(E[1:].sum()), G=float(E.sum()), E_minus1=float(np.polyval(E[::-1], -1)),
                    maxabs=float(np.abs(E).max()), kmax=int(np.argmax(np.abs(E))), sum_abs=float(np.abs(E).sum()),
                    cancellation=float(E.sum()/np.abs(E).sum()), min_Ez=float(Ez.min()), argmin_z=float(z[np.argmin(Ez)]),
                    G_over_H=float(E.sum()/E[0]), sign_check=m.get('sign_check'))
json.dump(summ, open('r34_controls_summary.json', 'w'), indent=1)
for k, v in summ.items(): print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != 'sign_check'})
fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
a = ax[0, 0]
for X, col in ((10000, 'tab:blue'), (100000, 'tab:green'), (1000000, 'tab:orange'), (10000000, 'tab:red')):
    r = [q for q in R['rows'] if q['X'] == X][0]; E = np.array(r['Ek']); k = np.arange(len(E)); s = np.abs(E) > 1e-6
    a.semilogy(k[s], np.abs(E[s]), '-', color=col, lw=1, label=f'X = {X:.0e}')
    a.scatter(k[s][E[s] > 0], np.abs(E[s][E[s] > 0]), marker='^', color=col, s=18); a.scatter(k[s][E[s] < 0], np.abs(E[s][E[s] < 0]), marker='v', color=col, s=18)
    a.scatter([1], [abs(E[1])], s=90, facecolors='none', edgecolors='k')
a.set_xlabel('layer k = ω(u)+ω(v) of the reduced ratio v/u'); a.set_ylabel('|E_k|  (▲ > 0, ▼ < 0)')
a.set_title('(a) ζ: ratio-class layers of the Fejér square; circled = prime class (Codex)', fontsize=10); a.legend(fontsize=8)
a = ax[0, 1]
for X, col in ((100000, 'tab:green'), (1000000, 'tab:orange'), (10000000, 'tab:red')):
    cv = R['curve'][str(X)]; z = np.array(cv['z']); a.semilogy(z, cv['E'], color=col, label=f'E_X(z), X = {X:.0e}')
    m = np.array(cv['SD_main']); s = m > 0; a.semilogy(z[s], m[s], '--', color=col, lw=.8)
a.set_xlabel('z'); a.set_ylabel('E_X(z) = mean of ‖f_c‖² over random signs (mean z)')
a.set_title('(b) layer polynomial E_X(z) ≥ 0; dashed: Selberg–Delange main ∝ 1/Γ(−z)² (zero at z=0,1)', fontsize=9)
a.legend(fontsize=8)
a = ax[1, 0]
X = np.array([q['X'] for q in R['rows']])
for key, lab, col in (('sum_abs', 'Σ_k |E_k|', 'k'), ('maxabs', 'max_k |E_k|', 'gray'), ('E1_prime', '|E_1| prime class', 'tab:blue'),
                      ('E0_H', 'E_0 = H (Haar)', 'tab:green'), ('G_E1', 'E(1) = G (real μ)', 'tab:red')):
    a.loglog(X, [abs(q[key]) for q in R['rows']], 'o-', color=col, label=lab, ms=4)
a.loglog(X, [q['G_E1'] for q in R['twist_rows']], 's--', color='tab:purple', label='G, Codex sparse twist', ms=4)
a.loglog(X, X, ':', color='gray', lw=.7); a.loglog(X, np.sqrt(X)/np.log(X)*50, ':', color='tab:blue', lw=.7)
a.set_xlabel('X'); a.legend(fontsize=8); a.set_title('(c) the payback: layers ≍ X, prime class ≍ √X/log X, total G = O(1) in range', fontsize=10)
a = ax[1, 1]
for nm, col in (('ZETA', 'k'), ('ZETA2', 'tab:purple'), ('FTH0.1', 'tab:orange'), ('EPS', 'tab:green'), ('DH', 'tab:blue'), ('BEUR', 'tab:brown')):
    E = np.array(C['models'][nm]['Ek']); k = np.arange(len(E)); s = np.abs(E) > 1e-9
    a.semilogy(k[s], np.abs(E[s])/np.abs(E).sum(), 'o-', ms=3, color=col, label=f"{nm}: G/Σ|E_k| = {summ[nm]['cancellation']:.1e}")
a.set_xlabel('layer k (Ω(u)+Ω(v))'); a.set_ylabel('|E_k| / Σ|E_k|'); a.legend(fontsize=7.5)
a.set_title('(d) controls at X = 3·10⁴, each with its own μ: the payback is generic to Dirichlet inverses', fontsize=9.5)
fig.suptitle('Round 34 — ratio classes of the R33 Fejér square (Codex grouping, layer generating polynomial); no zeros as inputs', fontsize=11)
fig.tight_layout(); fig.savefig('fig1_ratio_classes.png', dpi=130)
