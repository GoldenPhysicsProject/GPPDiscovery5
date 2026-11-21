"""Round 32 figures (from r32_series.npz, r32_results.json, r32_selberg.json)."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
z = np.load('r32_series.npz'); R = json.load(open('r32_results.json')); S = json.load(open('r32_selberg.json'))
col = dict(ZETA='k', ZETA2='tab:purple', **{'FTH0.1': 'tab:orange', 'FTH0.25': 'tab:red'}, EPS='tab:blue', DH='tab:green', BEUR='tab:brown')
lab = dict(ZETA='ζ', ZETA2='ζ² (F₀: on-line, double zeros)', **{'FTH0.1': 'F_θ, θ=0.1', 'FTH0.25': 'F_θ, θ=0.25'},
           EPS='Epstein x²+5y²', DH='Davenport–Heilbronn', BEUR='Beurling (DMV-type)')

# fig1: the Moebius floor
fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
for k in ('FTH0.25', 'EPS', 'FTH0.1', 'ZETA2', 'DH', 'ZETA'):
    L = z[f'{k}_Lg']; m = z[f'{k}_m1t']; s = L >= 3
    ax[0].plot(L[s], m[s], color=col[k], lw=0.6 if k != 'ZETA' else 1.0, label=lab[k])
ax[0].set_yscale('symlog', linthresh=0.02); ax[0].axhline(0, color='gray', lw=0.5)
ax[0].set_xlabel('L = log X'); ax[0].set_ylabel('m̃₁(L) = Σ_{d≤X} μ_M(d) d^{-1/2}(L − log d) − L/b(½) + b′(½)/b(½)²')
ax[0].set_title('Möbius floor m̃₁, own μ_M per model (X ≤ 10⁹; EPS, DH 10⁷)', fontsize=10)
ax[0].legend(fontsize=8, loc='upper left')
for k in ('FTH0.25', 'EPS', 'FTH0.1', 'ZETA2', 'DH', 'ZETA', 'BEUR'):
    g = R[k]['growth_m1t']; ax[1].semilogy(g['L'], g['W'], color=col[k], marker='.', ms=3, lw=1, label=f"{lab[k]}")
ax[1].set_xlabel('L (window [L, 1.5L])'); ax[1].set_ylabel('max |m̃₁| on the window')
ax[1].set_title('Window envelope: ζ flat ≤ 0.022; ζ² ∝ L; F_θ, EPS grow;\nDH flat (zeros too high); Beurling = left-of-line transient', fontsize=10)
ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig1_mobius_floor.png', dpi=130); plt.close()

# fig2: zeta decomposition D = T1 + T2 + T3 - E_mu + R0
L = z['ZETA_Lg']; s = L >= 4
fig, ax = plt.subplots(figsize=(12, 5))
for nm, v, c in (('4√X(1 − m⁽¹⁾)  [Möbius, σ=1]', z['ZETA_T1'], 'tab:blue'), ('−ζ′(½)·m̃₁  [Möbius floor]', z['ZETA_T3'], 'tab:orange'),
                 ('−2ζ(½)·m̃₂', z['ZETA_T2'], 'tab:green'), ('−E_μ (residual)', -z['ZETA_Emu'], 'tab:gray')):
    ax.plot(L[s], v[s], lw=0.5, color=c, label=nm)
ax.plot(L[s], z['ZETA_D'][s], 'k', lw=1.4, label='D(L) (R31 floor) = sum + R₀')
ax.axhline(0, color='gray', lw=0.5); ax.set_xlabel('L'); ax.legend(fontsize=8, ncol=2)
ax.set_title('ζ: the floor written in Möbius coordinates (identity A3). Every piece is two-sided; only the sum is one-sided (≥ +0.0216)')
plt.tight_layout(); plt.savefig('fig2_zeta_mobius_decomposition.png', dpi=130); plt.close()

# fig3: fixed-level Selberg mollifier
fig, ax = plt.subplots(1, 2, figsize=(15, 5))
for j, nm in enumerate(('ZETA', 'FTH')):
    c = S['curves'][nm]; Lc = np.array(c['L'])
    ax[j].plot(Lc, c['D'], 'k', lw=1.3, label='D (true floor)')
    for y, cc in ((1000, 'tab:blue'), (100000, 'tab:red')):
        ax[j].plot(Lc, c[f'Dy_{y}'], color=cc, lw=0.8, label=f'D_y, y = {y:g}')
        ax[j].axvline(math.log(y), color=cc, ls=':', lw=0.8)
    ax[j].set_yscale('symlog', linthresh=0.05); ax[j].set_xlabel('L'); ax[j].legend(fontsize=8)
    ax[j].set_title(('ζ' if nm == 'ZETA' else 'F_θ, θ = 0.25') + ': Selberg-mollified D_y = D for X ≤ y,\nthen an explicit (zero-free) drift', fontsize=10)
plt.tight_layout(); plt.savefig('fig3_selberg_fixed_level.png', dpi=130); plt.close()

# fig4: |mu_M| on dyadic blocks
fig, ax = plt.subplots(figsize=(9, 5))
for k in ('ZETA', 'ZETA2', 'FTH0.1', 'FTH0.25', 'EPS', 'DH'):
    v = np.array(R[k]['dyadic_max_mu'], float); j = np.arange(len(v)); ok = v > 0
    ax.semilogy(j[ok], v[ok], color=col[k], marker='.', label=lab[k])
ax.set_xlabel('j  (block 2^j ≤ n < 2^{j+1})'); ax.set_ylabel('max |μ_M(n)|'); ax.legend(fontsize=8)
ax.set_title('Counted-once on the inverse side: |μ| ≤ 1 only for ζ (and Beurling, by construction)')
plt.tight_layout(); plt.savefig('fig4_mu_growth.png', dpi=130); plt.close()
print('D mean (ZETA, L in [4, Lmax])', float(np.mean(z['ZETA_D'][s])))
