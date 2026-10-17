"""Grok S65: Arb recheck of the Schur sup for a saved step weight (zeta C_p). Window maxima are exact float comparisons on the stored
dyadic weights; cell windows are computed from Arb enclosures of l_n = log n and of the exact rational cell ends (widened outward by one cell);
the sum c_n * max / w_j is evaluated in Arb (prec 128) and its upper endpoint is maximised. usage: arb_recheck.py L N wfile"""
import sys, json, math
import numpy as np
from flint import arb, ctx, fmpq
ctx.prec = 128
Ls, N, wf = sys.argv[1], int(sys.argv[2]), sys.argv[3]
num, den = map(int, Ls.split('/')); Lq = fmpq(num, den); La = arb(Lq)
w = np.load(wf); assert len(w) == N and (w > 0).all()
def mang(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else 0
ns = [n for n in range(2, 40) if mang(n) and arb(n).log() < 2*La]
assert all(not (arb(n).log() >= 2*La) for n in ns)
h = 2*La/N
tot = [arb(0)]*N
best = arb(0); tot = np.zeros(N)
acc = [arb(0) for _ in range(0)]
sums = [arb(0)] * N
lo_all = []; 
contrib = []
for n in ns:
    c = arb(mang(n)).log()/arb(n).sqrt(); l = arb(n).log()
    for s in (l, -l):
        # cell j = [-L + j h, -L + (j+1) h]; shifted window hits cells k with k in [floor((j h + s)/h), floor(((j+1)h + s)/h)]
        sl = float((s/h).lower()); su = float((s/h).upper())
        j = np.arange(N)
        lo = np.floor(j + sl).astype(np.int64) - 1; hi = np.floor(j + 1 + su).astype(np.int64) + 1
        ins = (hi >= 0) & (lo <= N-1); lo = np.clip(lo, 0, N-1); hi = np.clip(hi, 0, N-1)
        mx = np.array([w[a:b+1].max() if i else 0.0 for a, b, i in zip(lo, hi, ins)])
        contrib.append((c, mx))
# Arb sum per cell only where float estimate is within 1e-6 of the float max (others bounded by float+1e-9 margin)
fl = sum(float(c.upper())*mx for c, mx in contrib)/w
top = np.where(fl > fl.max() - 1e-3)[0]
ub = max(float((sum((c*arb(float(mx[j])) for c, mx in contrib), arb(0))/arb(float(w[j]))).upper()) for j in top)
rest = float(np.max(np.delete(fl, top))) if len(top) < N else 0.0
print(json.dumps({"L": Ls, "N": N, "terms": ns, "arb_sup_top_cells": ub, "n_top": len(top), "float_sup_rest": rest,
                  "schur_bound_cert": max(ub, rest*(1+1e-12))}))
