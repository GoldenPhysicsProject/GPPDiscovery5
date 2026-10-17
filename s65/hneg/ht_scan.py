"""Grok S64 H64c stage 1: zeros of the de Bruijn-Newman function H_t(x) = int_0^inf e^{t u^2} Phi(u) cos(x u) du
for t < 0 (backward heat flow), x in [0, XMAX]. H_0(x) = xi(1/2 + i x/2)/8, so zeta zeros gamma <-> x = 2 gamma.
Composite Gauss-Legendre (panel 0.01, 20 nodes) on [0, 1.3] in Arb (prec 400). NUMERICAL (quadrature not rigorous).
Real zeros: sign changes on a grid of step DX + bisection. Count vs zeta zeros below XMAX/2 => missing = 2 x #complex pairs.
usage: ht_scan.py t XMAX DX out.json"""
import sys, json, time, math
from flint import arb, acb, ctx
import mpmath as mp
import os; ctx.prec = int(os.environ.get("HPREC", "400")); mp.mp.dps = 200
t = arb(sys.argv[1]); XMAX = float(sys.argv[2]); DX = float(sys.argv[3]); out = sys.argv[4]
PI = arb.pi()
def Phi(u):
    s = arb(0); n = 1
    while True:
        e4 = (4*u).exp(); term = (2*PI**2*n**4*(9*u).exp() - 3*PI*n**2*(5*u).exp())*(-PI*n*n*e4).exp()
        s += term
        if n > 1 and abs(float(term.mid())) < 1e-200: break
        n += 1
    return s

# Gauss-Legendre nodes on [-1,1] via mpmath
import mpmath
gx = []; gw = []
n = 20
for k in range(1, n+1):
    x = mpmath.cos(mpmath.pi*(k-mpmath.mpf(1)/4)/(n+mpmath.mpf(1)/2))
    for _ in range(100):
        p0, p1 = mpmath.mpf(1), x
        for j in range(2, n+1): p0, p1 = p1, ((2*j-1)*x*p1-(j-1)*p0)/j
        dp = n*(x*p1-p0)/(x*x-1); dx = p1/dp; x -= dx
        if abs(dx) < mpmath.mpf(10)**(-190): break
    gx.append(x); gw.append(2/((1-x*x)*dp*dp))
H = mpmath.mpf(1)/100; U = []; C = []
for pi_ in range(130):
    a = pi_*H
    for x, w in zip(gx, gw):
        u = arb(mpmath.nstr(a + H*(x+1)/2, 190)); ww = arb(mpmath.nstr(H*w/2, 190))
        U.append(u); C.append(ww*Phi(u)*(t*u*u).exp())
def Ht(x):
    x = arb(x) if not isinstance(x, arb) else x
    return sum((c*(x*u).cos() for u, c in zip(U, C)), arb(0))
def Htc(z):
    return sum((c*(z*u).cos() for u, c in zip(U, C)), acb(0))
t0 = time.time()
grid = [i*DX for i in range(1, int(XMAX/DX)+1)]
vals = [Ht(repr(x)) for x in grid]
sg = [1 if v > 0 else (-1 if v < 0 else 0) for v in vals]
zeros = []
for i in range(len(grid)-1):
    if sg[i]*sg[i+1] < 0:
        a, b, fa = arb(repr(grid[i])), arb(repr(grid[i+1])), vals[i]
        for _ in range(45):
            m = (a+b)/2; fm = Ht(m)
            if (fm > 0) == (fa > 0): a, fa = m, fm
            else: b = m
        zeros.append(float(((a+b)/2).mid()))
# local minima of |H|/scale without sign change: candidates for complex pairs
lv = [abs(float(v.mid())) for v in vals]
cand = [grid[i] for i in range(1, len(grid)-1) if lv[i] < lv[i-1] and lv[i] < lv[i+1] and sg[i-1] == sg[i] == sg[i+1]]
json.dump({"t": sys.argv[1], "XMAX": XMAX, "DX": DX, "real_zeros_x": zeros, "nreal": len(zeros),
           "absmin_nosign_candidates": cand, "H_at_1": str(vals[0]), "wall_s": time.time()-t0}, open(out, "w"))
print("t", sys.argv[1], "nreal", len(zeros), "cands", len(cand), "wall", time.time()-t0)
