"""Round 30-b figures: fig1 Z_M(t) on the principal series (window mechanism); fig2 Pontryagin index kappa(L)."""
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
Z = np.load('r30b_window_Z.npz'); ts = Z['ts']; m = ts <= 70
fig, ax = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
for a, nm in zip(ax, ['ZETA', 'FTH', 'EPS', 'DH']):
    a.plot(ts[m], Z[nm][m], lw=0.9); a.axhline(0, color='k', lw=0.5); a.set_ylabel(nm)
ax[1].set_title('F_theta: Z = |Lambda_zeta(1/4+it)|^2/|G| >= 0 -- no sign change ever (window functional never fires)', fontsize=9)
ax[0].set_title('Z_M(t) = Lambda_M(1/2+it)/|Q^s G(s)| (real by shadow symmetry + reality); no zeros used', fontsize=9)
ax[-1].set_xlabel('t'); plt.tight_layout(); plt.savefig('fig1_window_Z.png', dpi=120)
C = json.load(open('r30b_casimir.json'))
plt.figure(figsize=(7, 4.5))
for nm in ['ZETA', 'FTH', 'EPS', 'DH']:
    k = C[nm]['0.05']['kappa_by_L']; plt.plot([r['L'] for r in k], [r['kappa'] for r in k], 'o-', label=nm)
plt.xlabel('L (dilation range, log scale)'); plt.ylabel('kappa = # negative eigenvalues (grid h=0.05)')
plt.title('Unitarity Gram of the scaling flow: zeta kappa = 0; every fake kappa > 0', fontsize=10)
plt.legend(); plt.tight_layout(); plt.savefig('fig2_pontryagin_kappa.png', dpi=120)
