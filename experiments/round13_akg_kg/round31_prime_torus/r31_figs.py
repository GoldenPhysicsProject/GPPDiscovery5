"""Round 31 figures. fig1: the prime torus (winding on T^2, T^3, integers as the lattice face). fig2: the two faces and
turn-vs-stretch. fig3: the integer floor D(4 log k), four models. fig4: Bagchi limit law (time vs torus)."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
TH = 0.25
fl = json.load(open('r31_floor.json')); tor = json.load(open('r31_torus.json'))
l2, l3, l5 = math.log(2), math.log(3), math.log(5)

# ---------------- fig1 ----------------
fig = plt.figure(figsize=(17, 6.4))
ax = fig.add_subplot(1, 3, 1)
t = np.linspace(0, 300, 300000); x = np.mod(t*l2/(2*np.pi), 1); y = np.mod(t*l3/(2*np.pi), 1)
br = np.where((np.abs(np.diff(x)) > .5) | (np.abs(np.diff(y)) > .5))[0] + 1
for seg in np.split(np.arange(len(t)), br):
    ax.plot(x[seg], y[seg], lw=.7, color=plt.cm.viridis(t[seg[0]]/300))
ax.plot(0, 0, 'r*', ms=14); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')
ax.set_xlabel('turn of the 2-circle   (t·log2 / 2π mod 1)'); ax.set_ylabel('turn of the 3-circle   (t·log3 / 2π mod 1)')
ax.set_title('Prime torus, two circles: the line t ↦ (2^{-it}, 3^{-it})\nnever closes (log3/log2 irrational), fills the square\n(colour = time t ∈ [0, 300]; ★ = the base point t = 0)', fontsize=9)
ax = fig.add_subplot(1, 3, 2, projection='3d')
t = np.linspace(0, 120, 40000); P = np.mod(np.outer(t, [l2, l3, l5])/(2*np.pi), 1)
br = np.where(np.abs(np.diff(P, axis=0)).max(1) > .5)[0] + 1
for seg in np.split(np.arange(len(t)), br):
    ax.plot(P[seg, 0], P[seg, 1], P[seg, 2], lw=.5, color=plt.cm.plasma(t[seg[0]]/120))
ax.scatter([0], [0], [0], color='r', marker='*', s=120)
ax.set_xlabel('2-circle'); ax.set_ylabel('3-circle'); ax.set_zlabel('5-circle')
ax.set_title('Three circles (2, 3, 5): same winding in a cube\nwhose opposite faces are glued (t ∈ [0, 120])', fontsize=9)
ax = fig.add_subplot(1, 3, 3)
u = np.array([l2, l3]); u /= np.linalg.norm(u)
pts = [(a, b) for a in range(0, 7) for b in range(0, 5)]
for a, b in pts:
    n = 2**a*3**b; proj = (a*l2 + b*l3)/np.hypot(l2, l3)
    ax.plot(a, b, 'o', color='k', ms=4)
    foot = proj*u
    ax.plot([a, foot[0]], [b, foot[1]], ':', color='gray', lw=.6)
    ax.plot(*foot, '|', color='C3', ms=8, mew=1.5)
    if n <= 72: ax.annotate(str(n), (a, b), textcoords='offset points', xytext=(3, 3), fontsize=7)
L = np.array([0, 9]); ax.plot(L*u[0], L*u[1], color='C3', lw=1.2)
ax.set_xlim(-.5, 7); ax.set_ylim(-.5, 5.5); ax.set_aspect('equal')
ax.set_xlabel('a  (power of 2)'); ax.set_ylabel('b  (power of 3)')
ax.set_title('The other face: integers n = 2^a 3^b are lattice points;\nprojecting onto the direction (log2, log3) gives log n.\nIrrational slope ⇒ no two land together: each integer counted once', fontsize=9)
plt.tight_layout(rect=(0, 0, 1, 0.97)); plt.savefig('fig1_prime_torus.png', dpi=140); plt.close()

# ---------------- fig2 ----------------
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
ax = axs[0]
lx = np.linspace(-3, 3, 600); xx = np.exp(lx)
def theta(x):
    n = np.arange(1, 60)[:, None]; return 1 + 2*np.exp(-np.pi*n**2*x[None, :]).sum(0)
ax.plot(lx, xx**0.25*theta(xx), label='x^{1/4} θ(x)   (points face)')
ax.plot(lx, (1/xx)**0.25*theta(1/xx), '--', label='x^{-1/4} θ(1/x)   (waves face, mirrored)')
ax.axvline(0, color='k', lw=.8); ax.text(0.05, ax.get_ylim()[1]*0.9 if False else 2.5, 'mirror x = 1', fontsize=8)
ax.set_xlabel('log x   (left: 0 < x < 1,  right: 1 < x < ∞)'); ax.set_ylim(0.8, 3)
ax.set_title('Two faces in space: with the half-density weight the theta curve\nis exactly symmetric about x = 1 (Poisson/Jacobi; the curves coincide)', fontsize=9); ax.legend(fontsize=8)
ax = axs[1]
ax.axvspan(1, 1.6, color='C2', alpha=.15); ax.axvspan(.5, 1, color='C0', alpha=.12)
ax.axvline(.5, color='k', lw=2, label='mirror line Re s = ½ (FE of both)')
ax.axvline(1, color='C2', lw=1, ls='--', label='ζ: L∞ edge of torus model (σ = 1)')
ax.axvline(.5 + .012, color='C0', lw=1, ls=':', label='ζ: L² edge of torus model (σ = ½): touches the mirror')
ax.axvline(1 + TH, color='C3', lw=1, ls='--', label='F_θ: L∞ edge (1+θ)')
ax.axvline(.5 + TH, color='C3', lw=2, ls=':', label='F_θ: L² edge (½+θ)  ← its zeros sit here')
ax.axvline(.5 - TH, color='C3', lw=2, ls=':', alpha=.5, label='F_θ mirror zeros (½−θ)')
ax.set_xlim(0, 1.6); ax.set_ylim(0, 1); ax.set_yticks([]); ax.set_xlabel('Re s')
ax.set_title('Two faces in frequency: where the torus model lives\n(counted once a_p=1 puts the L² edge ON the mirror; F_θ\'s a_p=p^θ+p^{-θ} pushes it to ½+θ)', fontsize=9)
ax.legend(fontsize=7, loc='upper right')
ax = axs[2]
for key, c, lab in (('ZETA_0.6', 'C0', 'ζ at σ=0.6 (inside torus strip): flat = pure turning'),
                    ('FTH_0.6', 'C3', 'F_θ at σ=0.6 (mirror < σ < ½+θ): drifts = stretch'),
                    ('FTH_0.85', 'C1', 'F_θ at σ=0.85 (inside its own torus strip): flat'),
                    ('DH_0.6', 'C4', 'Davenport–Heilbronn at σ=0.6')):
    if key not in tor['T3']: continue
    r = tor['T3'][key]; H = [q['H'] for q in r]; m = [q['mean'] for q in r]
    ax.semilogx(H, m, 'o-', color=c, label=lab)
    if key == 'FTH_0.6': ax.semilogx(H, [q['pred_drift'] for q in r], 'x--', color=c, lw=.8, label='prediction (½+θ−σ)·log(1.05H/2π)  [|χ| ≠ 1]')
ax.set_xlabel('height H (samples t ∈ [H, 1.1H])'); ax.set_ylabel('mean of log|Z(σ+it)|')
ax.set_title('"Turns, never stretches", checked along the flow:\nmean log|Z| stationary ⇔ the flow is a rotation of the torus model', fontsize=9)
ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig('fig2_faces_turn_vs_stretch.png', dpi=140); plt.close()

# ---------------- fig3 ----------------
fig, axs = plt.subplots(1, 2, figsize=(15, 5))
cols = dict(ZETA='C0', FTH='C3', EPS='C1', DH='C4')
for name, m in fl['models'].items():
    axs[0].plot(m['Lgrid'], m['Dgrid'], color=cols[name], lw=.8, label=f'{name}: D(L), real L')
    axs[0].plot(4*np.log(np.array(m['k'])), m['D'], '.', color=cols[name], ms=4)
axs[0].axhline(0, color='k', lw=.6); axs[0].set_ylim(-3, 1.2); axs[0].set_xlabel('L = log X  (dots: L = 4 log k, the F13S samples)'); axs[0].set_ylabel('D(L)')
axs[0].set_title('Integer floor D(4 log k), from the integers alone (explicit formula, no zeros)\nζ: X ≤ 10^8 (k ≤ 100);  EPS, DH: X ≤ 10^7', fontsize=9); axs[0].legend(fontsize=8)
for name, m in fl['models'].items():
    W = m['windows']; axs[1].plot([w['K'] for w in W], [w['min'] for w in W], 'o-', color=cols[name], label=f'{name}: min over [K, K^{{3/2}}]')
    axs[1].plot([w['K'] for w in W], [w['frac_ge0'] for w in W], 's:', color=cols[name], ms=3, alpha=.6)
axs[1].axhline(0, color='k', lw=.6); axs[1].set_xlabel('K'); axs[1].set_ylim(-5, 1.1)
axs[1].set_title('Window minima (solid) and the fraction of k in the window with D ≥ 0 (dotted).\nF_θ keeps ≈ ½ of its samples ≥ 0: a positive-density floor cannot separate it', fontsize=9)
axs[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig3_integer_floor.png', dpi=140); plt.close()
print('figs ok')
