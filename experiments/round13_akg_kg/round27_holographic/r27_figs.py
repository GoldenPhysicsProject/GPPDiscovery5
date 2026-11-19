"""Round 27 figures (from saved JSON + cheap recomputation). Writes fig*.png in this folder."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from r22_models import lam_table, TH
# fig1: four lifts
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, (b, g, title) in zip(axs, [(0.5, 14.13, 'zeta: on the seam, the quartet collapses (chi=-1 lifts = chi=+1 lifts)'),
                                   (0.75, 14.13, 'F_theta: split quartet, separate chi=-1 class')]):
    pts = {'++  s': (b, g), '--  s\u0304': (b, -g), '+-  1-s': (1 - b, -g), '-+  1-s\u0304': (1 - b, g)}
    ax.axvline(0.5, color='k', lw=1.2); ax.axvspan(0, 1, color='0.95')
    for k, (x, y) in pts.items():
        c = 'tab:blue' if k[:2] in ('++', '--') else 'tab:red'
        ax.plot(x, y, 'o', ms=13 if k[:2] in ('++', '--') else 8, color=c, alpha=0.8)
        ax.annotate(k, (x, y), xytext=((-70 if k[:2] in ('++','--') else 10) if b == 0.5 else 8, 6 if y > 0 else -14), textcoords='offset points', fontsize=10, color=c)
    ax.annotate('', xy=(1 - b, g), xytext=(b, g), arrowprops=dict(arrowstyle='<->', color='tab:red')) if b != 0.5 else None
    ax.annotate('', xy=(b, -g), xytext=(b, g), arrowprops=dict(arrowstyle='<->', color='tab:blue', ls='--'))
    ax.text(0.52, 0, 'D: s -> s\u0304 (sides)', fontsize=9, color='tab:blue')
    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-20, 20); ax.set_xlabel('Re s'); ax.set_ylabel('Im s'); ax.set_title(title, fontsize=9.5)
fig.suptitle('Four lifts (q,t): ++ = s, -- = s\u0304, +- = 1-s, -+ = 1-s\u0304;  blue chi=+1, red chi=-1;  seam Re s = 1/2 = fixed set of R_t: s -> 1-s\u0304')
fig.tight_layout(); fig.savefig('fig1_four_lifts.png', dpi=130); plt.close(fig)
# fig2: sector scan
sc = json.load(open('r27_scan_ZETA_FTH_EPS_DH_1.0-5.0.json')); Ls = [1.0, 2.0, 3.0, 4.0, 5.0]
fig, ax = plt.subplots(figsize=(8, 5))
cols = dict(ZETA='k', FTH='tab:red', EPS='tab:green', DH='tab:purple')
for m in cols:
    for p, ls in (('even', '-'), ('odd', '--')):
        ax.plot(Ls, [sc[f'{m}_L{L}'][p]['lam_min_rel'] for L in Ls], ls, marker='o', color=cols[m], label=f'{m} {p}')
ax.set_yscale('symlog', linthresh=1e-9); ax.axhspan(-1e-9, 1e-9, color='0.85', label='double-precision floor')
ax.axhline(0, color='k', lw=0.6); ax.set_xlabel('window L (primes n < e^L)'); ax.set_ylabel('lambda_min / max|W| (symlog)')
ax.set_title('Weil form on the two J-sectors of the Theta-real form (low-frequency Galerkin)'); ax.legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig('fig2_sector_scan.png', dpi=130); plt.close(fig)
# fig3: F_theta attribution
q = json.load(open('r27_quartet.json')); keys = list(q)
fig, ax = plt.subplots(figsize=(8, 4.2)); x = np.arange(len(keys))
ax.bar(x - 0.2, [q[k]['plus_part'] for k in keys], 0.4, label='chi=+1 squares  v^T(4aa^T)v', color='tab:blue')
ax.bar(x + 0.2, [-q[k]['minus_part'] for k in keys], 0.4, label='chi=-1 squares  -v^T(4bb^T)v', color='tab:red')
ax.plot(x, [q[k]['lam_min'] for k in keys], 'k*', ms=12, label='lambda_min (total)')
ax.set_xticks(x); ax.set_xticklabels([k.replace('FTH_', '') for k in keys]); ax.axhline(0, color='k', lw=0.6)
ax.set_title('F_theta: the negativity is the separated chi=-1 class of the split quartets'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig('fig3_fth_chi_attribution.png', dpi=130); plt.close(fig)
# fig4: charge per prime power + phase non-convergence
N = 600; lz = lam_table('ZETA', N); pp = np.nonzero(lz > 0)[0]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
for m, c in (('ZETA', 'k'), ('FTH', 'tab:red'), ('EPS', 'tab:green'), ('DH', 'tab:purple')):
    lm = lam_table(m, N); n = np.arange(2, N + 1)
    q = np.where(lz[2:] > 0, lm[2:] / np.where(lz[2:] > 0, lz[2:], 1), np.nan)
    axs[0].plot(n, q, '.', color=c, ms=4, label=f'{m}: Lambda_M/Lambda at prime powers')
    off = (lz[2:] == 0) & (np.abs(lm[2:]) > 1e-9)
    axs[0].plot(n[off], lm[2:][off], 'x', color=c, ms=3, alpha=0.5)
axs[0].axhline(1, color='k', lw=0.6); axs[0].set_ylim(-3, 9); axs[0].set_xlabel('n')
axs[0].set_title("charge per prime power ('counted once' = 1); x = mass off prime powers", fontsize=9); axs[0].legend(fontsize=7)
Nn = 2_000_000; sv = np.ones(Nn + 1, bool); sv[:2] = False
for p in range(2, int(Nn**0.5) + 1):
    if sv[p]: sv[p*p::p] = False
P = np.nonzero(sv)[0].astype(float)
for t, c in ((10.0, 'tab:blue'), (14.134725, 'tab:orange')):
    cum = np.cumsum(np.angle((1 - P**(-0.5 - 1j*t)) / (1 - P**(-0.5 + 1j*t))))
    idx = np.unique(np.logspace(0, math.log10(len(P) - 1), 400).astype(int))
    axs[1].semilogx(P[idx], cum[idx], color=c, label=f't = {t:g}')
axs[1].set_xlabel('X'); axs[1].set_ylabel('sum_{p<=X} arg gamma_p(1/2+it)')
axs[1].set_title('each local factor is a pure phase, but the product does not converge on the line', fontsize=9); axs[1].legend()
fig.tight_layout(); fig.savefig('fig4_charge_and_phase.png', dpi=130); plt.close(fig)
print('ok')
# fig5: no budget (eps-rigidity) and truncation
b = json.load(open('r27_budget.json'))
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
for k, c in (('L1.0_even', 'tab:blue'), ('L1.0_odd', 'tab:cyan'), ('L1.5_even', 'tab:red'), ('L1.5_odd', 'tab:orange')):
    d = b[k]['lam_min_vs_eps']; e = np.array([float(x) for x in d]); v = np.array(list(d.values())); o = np.argsort(e)
    axs[0].plot(e[o], v[o], 'o-', color=c, label=k)
axs[0].set_xscale('symlog', linthresh=1e-6); axs[0].set_yscale('symlog', linthresh=1e-10); axs[0].axhline(0, color='k', lw=0.6)
axs[0].axvline(0, color='k', lw=0.6, ls=':'); axs[0].set_xlabel('eps  (Lambda -> (1+eps) Lambda)'); axs[0].set_ylabel('lambda_min / max|W|')
axs[0].set_title('zeta: scale all primes by (1+eps)\nPSD window [-6%, +0.02%] at L=1 shrinks to |eps|<1e-6 at L=1.5', fontsize=8.5); axs[0].legend(fontsize=8)
for k, c in (('L1.5_even', 'tab:red'), ('L1.5_odd', 'tab:orange')):
    t = b[k]['truncation_lam_min_rel']; P = [int(x) for x in t]; v = list(t.values())
    axs[1].plot(P, v, 'o-', color=c, label=k)
axs[1].set_yscale('symlog', linthresh=1e-10); axs[1].axhline(0, color='k', lw=0.6)
axs[1].set_xlabel('keep primes p <= P (window has p in {2, 3})'); axs[1].set_ylabel('lambda_min / max|W|')
axs[1].set_title('truncate primes (L=1.5): charge kept,\npositivity lost until every prime in the window is present', fontsize=8.5); axs[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig('fig5_no_budget_and_truncation.png', dpi=130); plt.close(fig)
print('fig5 ok')
