"""Grok S65 H65a: rigorous weighted Schur test for C_p = sum_{n<e^{2L}} Lambda(n) n^{-1/2} (S_{log n}+S_{-log n}) on L^2[-L,L]
with an OPTIMISED piecewise-constant weight (nonlinear Perron / Collatz-Wielandt iteration), vs Mori VIII's w = 1+(2/3)(x/L)^2.
Schur test (Mori VIII sec. 2, Liu App. B.4): ||C_p|| <= ess sup_x R_w(x), R_w(x) = sum c_n [1_I(x+l_n) w(x+l_n) + 1_I(x-l_n) w(x-l_n)] / w(x).
Rigorous upper bound for step w on N uniform cells: for x in cell j, w(x +- l) <= max of w over all cells meeting [a_j +- l, a_{j+1} +- l]
(index range widened by 1 cell on each side to absorb float error in l and cell ends); 1_I <= 1 on partially covered windows.
Final sums: nonneg float terms, c_n rounded up from Arb, result multiplied by (1+1e-12) (> n_terms * 2u, Higham). => CERTIFIED upper bound.
Also optional: signed coefficients (DH) via |c_n| (bounds || |C| || >= ||C||).
usage: perron_schur.py L N iters [dh]"""
import sys, json, math, time
import numpy as np
from flint import arb
Ls = sys.argv[1]; N = int(sys.argv[2]); IT = int(sys.argv[3]); MODE = sys.argv[4] if len(sys.argv) > 4 else 'zeta'; DH = MODE.startswith('dh')
L = float(eval(Ls, {'log': math.log}))
def mang(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else 0
if not DH:
    ns = [n for n in range(2, int(math.exp(2*L))+1) if mang(n) and math.log(n) < 2*L]
    cn = {n: float((arb(mang(n)).log()/arb(n).sqrt()).upper()) for n in ns}
else:
    s5 = math.sqrt(5); kap = (math.sqrt(10-2*s5)-2)/(s5-1)
    chi = {0:(0,0),1:(1,0),2:(0,1),3:(0,-1),4:(-1,0)}
    NN = int(math.exp(2*L)) + 1
    a = {n: chi[n%5][0] + kap*chi[n%5][1] for n in range(1, NN+1)}
    b = {1: 0.0}
    for n in range(2, NN+1): b[n] = a[n]*math.log(n) - sum(b[d]*a[n//d] for d in range(2, n) if n % d == 0)
    ns = [n for n in range(2, NN+1) if math.log(n) < 2*L and abs(b[n]) > 1e-15]
    if MODE == 'dhpos': ns = [n for n in ns if b[n] > 0]
    if MODE == 'dhneg': ns = [n for n in ns if b[n] < 0]
    # Arb recomputation of the DH coefficients (kappa algebraic), upper endpoints of |b(n)|/sqrt(n)
    s5a = arb(5).sqrt(); kapa = ((10-2*s5a).sqrt()-2)/(s5a-1)
    aa = {n: arb(chi[n%5][0]) + kapa*chi[n%5][1] for n in range(1, NN+1)}
    ba = {1: arb(0)}
    for n in range(2, NN+1): ba[n] = aa[n]*arb(n).log() - sum((ba[d]*aa[n//d] for d in range(2, n) if n % d == 0), arb(0))
    for n in ns: assert (ba[n] > 0) == (b[n] > 0) and (ba[n] < 0) == (b[n] < 0)
    for n in range(2, NN+1):
        if n not in ns and math.log(n) < 2*L: assert ba[n].contains(0) and abs(b[n]) < 1e-15
    cn = {n: float((abs(ba[n])/arb(n).sqrt()).upper())*(1+1e-15) for n in ns}
h = 2*L/N; edges = -L + h*np.arange(N+1)
def windows(l):
    lo = np.floor((edges[:-1] + l + L)/h).astype(np.int64) - 1
    hi = np.floor((edges[1:] + l + L)/h).astype(np.int64) + 1
    inside = (hi >= 0) & (lo <= N-1)
    return np.clip(lo, 0, N-1), np.clip(hi, 0, N-1), inside
W = []
for n in ns:
    l = math.log(n)
    for s in (l, -l): W.append((cn[n],) + windows(s))
def winmax(w, lo, hi):
    # range max via sparse table
    K = int(math.log2(N)) + 1; tab = [w]
    for k in range(1, K):
        prev = tab[-1]; sh = 1 << (k-1)
        tab.append(np.maximum(prev, np.concatenate([prev[sh:], np.full(sh, -np.inf)])))
    ln = hi - lo + 1; k = np.floor(np.log2(ln)).astype(np.int64)
    out = np.empty(N)
    for kk in np.unique(k):
        m = k == kk; t = tab[kk]
        out[m] = np.maximum(t[lo[m]], t[hi[m] - (1 << kk) + 1])
    return out
def Chat(w):
    acc = np.zeros(N)
    for c, lo, hi, ins in W:
        mx = winmax(w, lo, hi); acc += np.where(ins, c*mx, 0.0)
    return acc
x = (edges[:-1] + edges[1:])/2
res = {"L": Ls, "N": N, "mode": MODE, "terms": ns}
w0 = 1 + (2/3)*(x/L)**2
# Mori's weight (step version is an upper bound on the smooth-weight Schur sup, up to O(h))
R0 = Chat(w0)/w0; res["schur_mori_step"] = float(R0.max()*(1+1e-12))
w = np.ones(N); best = np.inf; t0 = time.time(); hist = []
for it in range(IT):
    Cw = Chat(w); R = Cw/w; sup = float(R.max()*(1+1e-12))
    if sup < best: best = sup; wbest = w.copy()
    hist.append(sup)
    w = Cw + 2.0*w; w = (w + w[::-1])/2; w /= w.max()    # shifted nonlinear power step, symmetrised
res["schur_perron_cert"] = best; res["hist"] = hist[::max(1, IT//10)]; res["wall"] = time.time()-t0
res["w_min_over_max"] = float(wbest.min()/wbest.max())
np.save(f"w_L{Ls.replace('/','_')}_N{N}_{MODE}.npy", wbest)
print(json.dumps(res), flush=True)
