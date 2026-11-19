"""Round 15 core: bilinear Weil form between a base window and its log-shifts (the Z-analogue of Delta . Gamma_{F^n}).
Basis phi_i = Dirichlet sines on [-l/2, l/2] (14B weil_galerkin).  For a shift s (any real) define the M x M matrix
   B_s[i,j] = W-bilinear( phi_i , tau_s phi_j ),   tau_s phi(x) = phi(x - s)   (x = log scale, u = log x)
via the Round-13/14B linear functional on even autocorrelations
   Q(A) = C A(0) + 2 int_0^inf w(t)[A(0)-A(t)] dt + int_R A(t) k(t) dt - 2 sum_n c(n) n^-1/2 A(log n),
applied to Asym(t) = (X_ij(t-s) + X_ij(-t-s))/2,  X_ij(u) = int phi_i(x) phi_j(x+u) dx  (closed form, 14B Aclosed).
k(t) = sum over poles a of 2 cosh(a t) (complex pairs combined, real).  B_0 = 14B W exactly (checked in run).
Translation invariance:  W(tau_a f, tau_b h) = B_{b-a};  B_{-s} = B_s^T.  Integer/own-gamma data only; no zeros."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import _paths  # noqa: E402  (run from anywhere)
import sys, os, math, numpy as np
HERE = _paths.HERE
import weil_galerkin as wg, models14b as md
from scipy.special import digamma
XG, WG = np.polynomial.legendre.leggauss(48)

def X(l, M, u):
    """X_ij(u) = int phi_i(x) phi_j(x+u) dx, any real u."""
    if abs(u) >= l: return np.zeros((M, M))
    return wg.Aclosed(l, M, u) if u >= 0 else wg.Aclosed(l, M, -u).T

def ftau(tau=3.0):
    Z = md.zeta(); C = 2*digamma(0.25+0.5j*tau).real - 2*math.log(math.pi)
    return dict(name=f'FTAU{tau}', n=Z['n'], c=Z['c']*2*np.cos(tau*np.log(Z['n'])), C=C,
                w=lambda t: 2*np.cos(tau*t)*md.wz(t), kpole=lambda t: 4*np.cosh(t/2)*np.cos(tau*t))

def kpole(model):
    if 'kpole' in model: return model['kpole']
    P = model.get('poles', [])
    return lambda t: sum(wt*2*np.cosh(a*t) for a, wt in P) if P else 0.0*t

def panels(a, b, n):
    e = np.linspace(a, b, n+1); out = []
    for lo, hi in zip(e[:-1], e[1:]):
        out += [(0.5*(hi-lo)*x+0.5*(hi+lo), 0.5*(hi-lo)*w) for x, w in zip(XG, WG)]
    return out

def Bs(model, l, M, s, npan=24):
    """B_s for shift s (s may be negative)."""
    if s < 0: return Bs(model, l, M, -s, npan).T
    w = model['w']; kp = kpole(model)
    X0 = X(l, M, -s)                                   # = Asym(0)
    # arch: C*X0 + 2 int_0^inf w(t)[X0 - (X(t-s)+X(-t-s))/2] dt ; integrand support: t < s+l (beyond: w(t) X0)
    T = s + l
    arch = model['C']*X0
    from scipy.integrate import quad
    if np.any(X0): arch = arch + 2*X0*quad(lambda t: float(w(t)), T, np.inf, limit=400)[0]
    # nodes: geometric near 0 (w ~ 1/2t), dense near t = s (kinks of X at t-s = 0, +-l)
    br = sorted(set([0.0, T] + [x for x in (s-l, s, s+l, l, max(0, l-s)) if 0 < x < T]))
    for a, b in zip(br[:-1], br[1:]):
        if a == 0.0:
            ed = np.concatenate([[0], b*np.geomspace(1e-8, 1, 40)])
            nodes = sum((panels(lo, hi, 1) for lo, hi in zip(ed[:-1], ed[1:])), [])
        else: nodes = panels(a, b, npan)
        for t, wt in nodes:
            arch = arch + 2*wt*float(w(t))*(X0 - 0.5*(X(l, M, t-s) + X(l, M, -t-s)))
    # pole: int_R X(t-s) k(t) dt = int_{-l}^{l} X(u) k(u+s) du
    pole = np.zeros((M, M))
    for a, b in ((-l, 0.0), (0.0, l)):
        for u, wt in panels(a, b, npan): pole += wt*kp(u+s)*X(l, M, u)
    # prime: -2 sum c(n) n^-1/2 Asym(log n) = - sum c(n) n^-1/2 (X(log n - s) + X(-log n - s))
    prime = np.zeros((M, M)); ln = np.log(model['n'])
    sel = np.nonzero(np.abs(ln - s) < l)[0]
    for k in sel: prime += model['c'][k]/math.sqrt(model['n'][k])*X(l, M, ln[k]-s)
    if s < l:
        for k in np.nonzero(ln < l - s)[0]: prime += model['c'][k]/math.sqrt(model['n'][k])*X(l, M, -ln[k]-s)
    return dict(arch=arch, pole=pole, prime=prime, B=arch + pole - prime, N=X(l, M, -s))

MODELS = dict(ZETA=md.zeta, EPS=md.eps, DH=md.dh, BEUR2=md.beur2, FTH=lambda: md.fth(0.2), FTAU=lambda: ftau(3.0))
def Mof(l, fmax): return int(math.ceil(fmax*l/math.pi)) + 8
