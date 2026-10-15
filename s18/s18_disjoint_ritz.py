"""S18 claim 2 (+6): sup_g |F_g(log2)|/||g||^2 on the Legendre*(1-u^2)^4 basis (Ritz, NUMERICAL) with Arb-certified Rayleigh
quotients of explicit rational g (CERTIFIED lower bounds for the sup). Theory (PROVED): <= 1/2 for h <= log2 (disjoint strips),
= cos(pi/(ceil(2h/t)+1)) numerical radius of the nilpotent compression P_h U_t P_h (Haagerup-de la Harpe) in general.
usage: python s18_disjoint_ritz.py M1 M2 ..."""
import sys
from s18_ritz_common import *
t = arb(2).log()
HS = [("0.36", arb("0.36")), ("0.45", arb("0.45")), ("h*=log3/2", HSTAR), ("0.65", arb("0.65")), ("0.69", arb("0.69")), ("0.75 (>log2)", arb("0.75")), ("1.0 (>log2)", arb(1))]
for M in map(int, sys.argv[1:]):
    for hn, h in HS:
        B = build(M, h); N, T2 = B["N"], B["T"][2]
        n = None; k = 1
        while not (k*t >= 2*h): k += 1
        rad = mp.cos(mp.pi/(k+1))
        for sn, I in SECT(M):
            ev, V, Li = geig(T2, N, I)
            j = max(range(len(ev)), key=lambda i: abs(ev[i])); v = Li.T*V[:, j]
            c, _ = rq(T2, N, I, v)
            print(f"M={M:2d} h={hn:13s} {sn:4s}  Ritz max|F(log2)|/N = {mp.nstr(abs(ev[j]),12):16s} (sign {'+' if ev[j]>0 else '-'})  "
                  f"Arb RQ = {c.str(12):30s} theory sup = cos(pi/{k+1}) = {mp.nstr(rad,10)}  <=1/2: {bool(abs(c) < arb(1)/2)}", flush=True)
