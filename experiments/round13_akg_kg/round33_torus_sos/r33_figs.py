"""Round 33 figures."""
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
z = np.load('r33_series.npz'); R = json.load(open('r33_results.json')); V = json.load(open('r33_validate_zeros.json'))
R0 = 0.04620998623083794
col = {'ZETA': 'k', 'ZETA2': 'tab:purple', 'FTH0.1': 'tab:orange', 'FTH0.25': 'tab:red', 'EPS': 'tab:green', 'DH': 'tab:blue', 'BEUR': 'tab:brown'}
lab = {'ZETA': r'$\zeta$', 'ZETA2': r'$\zeta^2$', 'FTH0.1': r'$F_{0.1}$', 'FTH0.25': r'$F_{0.25}$', 'EPS': 'EPS', 'DH': 'DH', 'BEUR': 'Beurling'}
fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
a = ax[0, 0]
for n in col:
    L = z[f'{n}_L']; s = L >= 3; a.semilogy(L[s], z[f'{n}_Qc'][s]/L[s], color=col[n], label=lab[n], lw=1.2)
a.axhline(V['rows'][-1]['kappa_mu'], color='k', ls=':', lw=1); a.text(13, 0.034, r'$\kappa_\mu$ (300 zeros, test data only)', fontsize=8)
a.set_xlabel('L = log X'); a.set_ylabel(r'$Q_c(L)/L$'); a.legend(fontsize=8, ncol=4, loc='upper left')
a.set_title(r'(a) Kronecker-line mass $Q_c=\int_0^L|M_{1/2}(e^u)-1/b(\frac{1}{2})|^2du$, per unit L', fontsize=10)
a = ax[0, 1]
for n in col:
    L = z[f'{n}_L']; s = L >= 3; a.plot(L[s], z[f'{n}_Rc'][s]/z[f'{n}_Hc'][s], color=col[n], lw=1)
a.axhline(0, color='gray', lw=.5); a.set_xlabel('L'); a.set_ylabel(r'$R_c/H_c$')
a.set_title(r'(b) off-diagonal remainder / Haar diagonal: $\|f_c\|^2_{\nu_L}=H_c+R_c$', fontsize=10)
a = ax[1, 0]
W = 4.0
for key, cc, nm, ref, rl in (('ZETA_PSI', 'k', r'$\psi$ side $(\Lambda-1)$', R0, r'$R_0=(\log\xi)^{\prime\prime}(\frac{1}{2})$ (integers + $\Gamma$, Arb)'),
                             ('ZETA', 'tab:blue', r'Möbius side', V['rows'][-1]['kappa_mu'], r'$\kappa_\mu$ (300 zeros, test data only)')):
    L = z[f'{key}_L']; Q = z[f'{key}_Qc']; k = int(round(W/(L[1]-L[0]))); x = L[k:]; sl = (Q[k:] - Q[:-k])/W
    s_ = x >= 6; a.plot(x[s_], sl[s_], color=cc, lw=1, label=nm + r': $[Q(L)-Q(L-4)]/4$')
    a.axhline(ref, color=cc, ls='--', lw=1, label=rl)
a.set_ylim(0, 0.08); a.set_xlabel('L'); a.legend(fontsize=7.5, loc='lower right')
a.set_title(r'(c) $\zeta$: windowed slopes of the Kronecker-line mass (no zeros used except the dotted $\kappa_\mu$)', fontsize=9)
a = ax[1, 1]
L = z['ZETA_L']; s = L >= 3
a.plot(L[s], np.abs(z['ZETA_mc'][s]), 'k', label=r'$|m(L)-cL|$ (R32 floor, shifted)')
a.plot(L[s], np.sqrt(L[s]*z['ZETA_Qc'][s]), 'tab:blue', label=r'SOS bound $\sqrt{L\,Q_c}$')
a.plot(L[s], np.sqrt(L[s]*z['ZETA_Gc'][s]), 'tab:cyan', lw=.7, label=r'$\sqrt{L\,G_c}$')
a.plot(L[s], 2*np.exp(L[s]/2), 'tab:red', ls='--', label=r'termwise transfer cost $2\sqrt{X}$ (Prop 32.3)')
a.set_yscale('log'); a.set_xlabel('L'); a.legend(fontsize=8); a.set_title(r'(d) $\zeta$: SOS bound vs the $\sqrt{X}$ route', fontsize=10)
fig.suptitle('Round 33 — torus sum-of-squares on the Fejér (Kronecker-line) measure; each model with its own μ; no zeros as inputs', fontsize=11)
fig.tight_layout(); fig.savefig('fig1_torus_sos.png', dpi=130)
