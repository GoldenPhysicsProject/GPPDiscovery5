"""Round 22, approach 3: de Branges / Hermite-Biehler frontier. E_h(z) = Xi_M(z + i h), Xi_M(z) = xi_M(1/2 + i z).
HB (|E_h(conj z)| < |E_h(z)| for Im z > 0) <=> |xi_M(1/2 + a + ix)| < |xi_M(1/2 + b + ix)| for 0 <= |a| < b with
a = y - h, b = y + h (FE + real coefficients). We tabulate f(sigma, x) = log|xi_M(1/2 + sigma + i x)| on a grid
(no zeros used) and report, for each h, whether some (x, y) violates HB, i.e. f(|y-h|, x) >= f(y+h, x)."""
import sys, json, numpy as np, mpmath as mp
mp.mp.dps = 20
def L(s, chi, q): return q**(-s)*mp.fsum(chi[r % q]*mp.zeta(s, mp.mpf(r)/q) for r in range(1, q + 1) if chi[r % q] != 0)
c4 = {0: 0, 1: 1, 2: 0, 3: -1}; c5 = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
c20 = {r: c4[r % 4]*c5[r % 5] for r in range(20)}
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1); dh = {0: 0, 1: 1, 2: k, 3: -k, 4: -1}
def xi_zeta(s): return s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def xi(name, s):
    if name == 'ZETA': return xi_zeta(s)
    if name == 'FTH': return xi_zeta(s - 0.25)*xi_zeta(s + 0.25)
    if name == 'DH': return (5/mp.pi)**((s + 1)/2)*mp.gamma((s + 1)/2)*L(s, dh, 5)
    if name == 'EPS':
        Z = mp.zeta(s)*L(s, c20, 20) + L(s, c4, 4)*L(s, c5, 5)
        return s*(s - 1)*(mp.sqrt(20)/(2*mp.pi))**s*mp.gamma(s)*Z
name = sys.argv[1]; XMAX = float(sys.argv[2])
if name == 'EPS':   # check genus identity against the direct lattice sum at s = 3
    direct = mp.fsum((m*m + 5*n*n)**-3 for m in range(-60, 61) for n in range(-30, 31) if (m, n) != (0, 0))
    print('EPS genus check', direct, mp.zeta(3)*L(3, c20, 20) + L(3, c4, 4)*L(3, c5, 5), flush=True)
sig = np.round(np.arange(0, 1.301, 0.05), 3); xs = np.round(np.arange(0.1, XMAX + 1e-9, 0.1), 2)
f = np.zeros((len(sig), len(xs)))
for j, x in enumerate(xs):
    for i, s in enumerate(sig):
        v = abs(xi(name, mp.mpf(0.5) + mp.mpf(s) + 1j*mp.mpf(x)))
        f[i, j] = float(mp.log(v)) if v > 0 else -1e300
H = [0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]
idx = {round(float(s), 3): i for i, s in enumerate(sig)}
out = dict(x_range=[float(xs[0]), float(xs[-1])], sigma_grid=sig.tolist(), h={})
for h in H:
    worst = (-1e9, None, None)
    for y in np.round(np.arange(0.05, 0.701, 0.05), 3):
        a, b = round(abs(y - h), 3), round(y + h, 3)
        if a not in idx or b not in idx: continue
        d = f[idx[a]] - f[idx[b]]          # >= 0 means violation
        j = int(np.argmax(d))
        if d[j] > worst[0]: worst = (float(d[j]), float(xs[j]), float(y))
    out['h'][str(h)] = dict(HB_violated=worst[0] >= 0, max_log_ratio=worst[0], at_x=worst[1], at_y=worst[2])
    print(name, h, out['h'][str(h)], flush=True)
json.dump(out, open(f'r22_debranges_{name}.json', 'w'), indent=1)
