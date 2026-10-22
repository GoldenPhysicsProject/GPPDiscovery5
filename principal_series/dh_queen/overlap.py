"""DH lowest even Ritz vector (s75/dhritz.py cos basis, L = HALF-width, f on [-L,L]) vs Queen's
f_Q(x)=e^{-i g0 x} sinh(b0 x) on [-L,L] (team/queen dh_sinh.py: same half-width convention, G supported on [0,2L], primes n<e^{2L}).
Re f_Q = cos(g0 x) sinh(b0 x) is ODD; Im f_Q = -sin(g0 x) sinh(b0 x) is EVEN. NUMERICAL (float64 grid + mpmath quad)."""
import json, glob, math, sys, numpy as np, mpmath as mp
from fractions import Fraction
g0, b0 = 85.699348, 0.308517
def load(fn):
    d = json.load(open(fn)); Ls = d['L']
    L = float(Fraction(Ls.replace('_','/'))) if '/' in Ls else float(Ls)
    v = np.array([float(Fraction(s)) for s in d['vec']]); return L, v, d
def evalf(L, v, x):
    a = (2*np.arange(len(v))+1)*np.pi/(2*L); return np.cos(np.outer(x, a)) @ v
def nodes(y):
    s = np.sign(y); s = s[s != 0]; return int(np.sum(s[1:] != s[:-1]))
files = sys.argv[1:] or sorted(glob.glob('/workspace/sm/GPPDiscovery5/s75/cosvec_DH_even_*.json'))
for fn in files:
    L, v, d = load(fn); K = len(v)
    x = np.linspace(-L, L, 400001); f = evalf(L, v, x)
    gE = -np.sin(g0*x)*np.sinh(b0*x); gO = np.cos(g0*x)*np.sinh(b0*x)
    dx = x[1]-x[0]; ip = lambda p, q: np.trapezoid(p*q, dx=dx)
    ovE = ip(f, gE)/math.sqrt(ip(f, f)*ip(gE, gE)); ovO = ip(f, gO)/math.sqrt(ip(f, f)*ip(gO, gO))
    ovC = math.sqrt(ip(f, gE)**2 + ip(f, gO)**2)/math.sqrt(ip(f, f)*(ip(gE, gE)+ip(gO, gO)))
    # projection of gE onto the K-mode even basis (Gram = L I): captured fraction
    a = (2*np.arange(K)+1)*np.pi/(2*L); c = np.array([ip(gE, np.cos(ai*x)) for ai in a])/L
    capt = L*np.sum(c**2)/ip(gE, gE)
    # tail: fraction of |v|^2 in modes with a_m > g0 + 10
    hi = a > g0 + 10; tail = float(np.sum(v[hi]**2)/np.sum(v**2))
    pk = int(np.argmax(np.abs(v))); 
    print(json.dumps(dict(file=fn.split('/')[-1], L=L, K=K, amax=round(a[-1],1), nodes_DHvec=nodes(f), nodes_ImQ=nodes(gE), nodes_ReQ=nodes(gO),
        pred_2Lg0_over_pi=round(2*L*g0/math.pi, 2), overlap_ImQ_even="%.3e"%ovE, overlap_ReQ_odd="%.3e"%ovO, overlap_complexQ="%.3e"%ovC,
        ImQ_captured_by_basis=round(capt, 8), peak_mode=pk, peak_freq=round(a[pk], 2), DHvec_tail_mass_a_gt_g0p10=tail, RQ=d["RQ"], mass_by_mode=[round(float(t),4) for t in (v[:6]**2/np.sum(v**2))])), flush=True)
