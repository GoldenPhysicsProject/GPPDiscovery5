"""DH lowest EVEN eigenvector along the L ladder (L = half-width; s75/dhritz.py cos basis a_m=(2m+1)pi/(2L), mode 0 = cos(pi x/2L)).
Reports: Ritz lambda_1, interior sign changes on [-L,L], mode-0 / mode-1 weights, normalized overlap with the even part of Queen's
q(x)=e^{-i g0 x} sinh(b0 x) (Im q = -sin(g0 x)sinh(b0 x); Re q is odd => overlap 0 by parity), and Golden Eye's cross term
X = fhat(g0 - i b0) * conj(fhat(g0 + i b0)) / ||f||^2, fhat(z) = int_{-L}^{L} f(x) e^{i z x} dx (= int f cos(zx) for even f).
mpmath dps 60, closed-form integrals. NUMERICAL."""
import json, sys, glob, mpmath as mp, numpy as np
from fractions import Fraction
mp.mp.dps = 60
g0, b0 = mp.mpf('85.699348'), mp.mpf('0.308517')
def fr(s): q = Fraction(s); return mp.mpf(q.numerator)/q.denominator
def cosint(a, z, L):   # int_{-L}^{L} cos(a x) cos(z x) dx, z complex
    return mp.sin((a-z)*L)/(a-z) + mp.sin((a+z)*L)/(a+z)
def sinsinh(a, L):     # int_{-L}^{L} cos(a x) * (-sin(g0 x) sinh(b0 x)) dx via exponentials: -Im-part trick
    # -sin(g x) sinh(b x) = Re[ i/2 (e^{(b+ig)x}... ] ; use: sin(gx)sinh(bx) = Im[cos((g - i b)x)]*(-1)?  do directly:
    # cos((g-ib)x) = cos gx cosh bx + i sin gx sinh bx  => sin gx sinh bx = Im cos((g-ib)x)
    return -mp.im(cosint(a, g0 - 1j*b0, L))
rows = []
for fn in sys.argv[1:]:
    d = json.load(open(fn)); Ls = d['L']; L = fr(Ls.replace('_', '/')) if '/' in Ls else mp.mpf(Ls)
    v = [fr(s) for s in d['vec']]; K = len(v); a = [(2*m+1)*mp.pi/(2*L) for m in range(K)]
    n2 = L*mp.fsum(t*t for t in v)
    fm = mp.fsum(vi*cosint(ai, g0 - 1j*b0, L) for vi, ai in zip(v, a))
    fp = mp.fsum(vi*cosint(ai, g0 + 1j*b0, L) for vi, ai in zip(v, a))
    X = fm*mp.conj(fp)/n2
    qn2 = mp.quad(lambda x: (mp.sin(g0*x)*mp.sinh(b0*x))**2, mp.linspace(-L, L, 400))
    ov = mp.fsum(vi*sinsinh(ai, L) for vi, ai in zip(v, a))/mp.sqrt(n2*qn2)
    xs = np.linspace(-float(L), float(L), 200001); f = np.cos(np.outer(xs, [float(t) for t in a])) @ np.array([float(t) for t in v])
    core = np.abs(xs) < float(L) - 1e-4; s = np.sign(f[core]); nodes = int(np.sum(s[1:] != s[:-1]))
    w = [float(t*t/(n2/L)) for t in v]
    lam = mp.mpf(d['RQ'].strip('[').split(' ')[0])
    r = dict(L=Ls, K=K, lam1=mp.nstr(lam, 6), nodes=nodes, w0=round(w[0], 4), w1=round(w[1], 4), overlap_evenq=mp.nstr(ov, 4),
             X_re=mp.nstr(mp.re(X), 4), X_abs=mp.nstr(abs(X), 4), lam_over_absX=mp.nstr(lam/abs(X), 4))
    rows.append(r); print(json.dumps(r), flush=True)
