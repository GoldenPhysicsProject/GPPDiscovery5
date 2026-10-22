"""S84: L2 overlap |<f,g>|/(|f||g|) of two even cosine-basis vectors on the same window (basis orthogonal, Gram L*I -> coefficient
overlap on the common modes), plus pole-channel overlap with cosh(u/2) and first-zero channel. NUMERICAL (mpmath)."""
import json, sys, mpmath as mp
from fractions import Fraction
mp.mp.dps = 220
def vec(fn):
    d = json.load(open(fn)); return [mp.mpf(Fraction(q).numerator)/Fraction(q).denominator for q in d['vec']], d['L'].replace('_','/')
def ov(a, b):
    n = min(len(a), len(b)); s = mp.fsum(a[i]*b[i] for i in range(n))
    return abs(s)/mp.sqrt(mp.fsum(x**2 for x in a)*mp.fsum(x**2 for x in b))
def cosh_ov(v, L):
    Lm = mp.mpf(Fraction(L).numerator)/Fraction(L).denominator; K = len(v)
    a = [(2*m+1)*mp.pi/(2*Lm) for m in range(K)]
    c = [2*(-1)**m*a[m]*mp.cosh(Lm/2)/(a[m]**2+mp.mpf(1)/4) for m in range(K)]   # <cos(a_m u), cosh(u/2)> on [-L,L]
    nc = mp.sqrt(Lm + mp.sinh(Lm))                                                    # ||cosh(u/2)||
    return abs(mp.fsum(v[m]*c[m] for m in range(K)))/(mp.sqrt(Lm*mp.fsum(x**2 for x in v))*nc)
if __name__ == '__main__':
    fs = sys.argv[1:]; V = [vec(f) for f in fs]
    for i, f in enumerate(fs): print(f.split('/')[-1], "cosh(u/2) overlap", mp.nstr(cosh_ov(*V[i]), 8))
    for i in range(len(fs)):
        for j in range(i+1, len(fs)):
            if V[i][1] == V[j][1]: print(fs[i].split('/')[-1], fs[j].split('/')[-1], "overlap", mp.nstr(ov(V[i][0], V[j][0]), 10))
