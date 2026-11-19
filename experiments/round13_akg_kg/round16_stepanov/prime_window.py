"""Round 16 Part 3a: the prime window of the one-shift Hasse block, out to the largest feasible shift.
g = first Dirichlet sine on [-l/2, l/2], l = 0.6, ||g|| = 1 (Round 15 basis vector phi_1); X = autocorrelation of g (closed form).
For a shift s > l (Round 15 hasse_core.Bs, entry [0,0]):
   P(s)   = sum_n c(n) n^{-1/2} X(log n - s)                 prime window (the 'N_k' of Hasse)
   Pi(s)  = int X(u) k(u + s) du                           pole main term (the 'q^k + 1')
   A(s)   = - int w(t) X(t - s) dt                          archimedean tail (O(e^{-s/2}))
   b_s    = A + Pi - P = W(g, tau_s g)                     (checked against hasse_core.Bs to 1e-8 for s <= 8)
   d      = W(g) (diagonal, prime-free 'genus')
Hasse block <=> |b_s| <= d.  RH-strength target: |P(s) - Pi(s)| = O(1)  (<=> unnormalised error O(e^{s/2})).
The literal ratio P(s)/e^{s/2} is reported as requested; it is bounded for zeta by Chebyshev (PROVED), see ROUND16 s3.
Also records c(n) >= 0 per model (Principal Series ROUND16PS, e648fd3). No zeros used."""
import sys, os, math, json, numpy as np
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import models16 as M16
sys.path.insert(0, os.path.join(HERE, '..', 'round15_hasse'))
import hasse_core as hc, models14b as md
L = 0.6; TH = 0.2; TAU = 3.0
def X(u):
    a = np.abs(u); r = ((L - a)*np.cos(np.pi*a/L) + (L/np.pi)*np.sin(np.pi*a/L))/L
    return np.where(a < L, r, 0.0)
xg, wg = np.polynomial.legendre.leggauss(80)
def gl(f, a, b): x = 0.5*(b-a)*xg + 0.5*(b+a); return 0.5*(b-a)*np.sum(wg*f(x))
wz = md.wz
KP = dict(ZETA=lambda t: 2*np.cosh(t/2), EPS=lambda t: 2*np.cosh(t/2), BEUR2=lambda t: 2*np.cosh(t/2), DH=lambda t: 0*t,
          FTH=lambda t: 2*np.cosh((0.5+TH)*t) + 2*np.cosh((0.5-TH)*t), FTAU=lambda t: 4*np.cosh(t/2)*np.cos(TAU*t))
EPSm = None
def W_of(name):
    if name in ('ZETA', 'BEUR2'): return wz
    if name == 'FTAU': return lambda t: 2*np.cos(TAU*t)*wz(t)
    if name == 'FTH': return lambda t: (np.exp(-(0.5-TH)*t) + np.exp(-(0.5+TH)*t))/(-np.expm1(-2*t))
    if name == 'EPS': return lambda t: np.exp(-t/2)/(-np.expm1(-t))
    if name == 'DH': return lambda t: np.exp(-1.5*t)/(-np.expm1(-2*t))

def data():
    ln, lam = M16.zeta_lam(); F = M16.fakes()
    return dict(ZETA=(ln, lam, None), FTH=(ln, lam, lambda u: 2*np.cosh(TH*u)), FTAU=(ln, lam, lambda u: 2*np.cos(TAU*u)),
                EPS=(*F['EPS'], None), DH=(*F['DH'], None), BEUR2=(*F['BEUR2'], None))

def sgrid(smax):
    return np.unique(np.round(np.concatenate([np.arange(1.0, min(smax, 10)+1e-9, 0.01), np.arange(10, min(smax, 14)+1e-9, 0.025),
                                               np.arange(14, smax+1e-9, 0.05)]), 6))

def main():
    D = data(); out = {}
    diag = {nm: float(hc.Bs(hc.MODELS[nm](), L, 1, 0.0)['B'][0, 0]) for nm in D}
    for nm, (ln, c, wfun) in D.items():
        smax = round(ln[-1] - L - 0.01, 2); S = sgrid(smax); kp = KP[nm]; w = W_of(nm)
        pos = dict(c_nonneg=bool(np.all((c if wfun is None else c*wfun(ln)) >= -1e-12)))
        cc = c if wfun is None else c*wfun(ln); neg = np.nonzero(cc < -1e-12)[0]
        pos.update(n_negative=int(len(neg)), first_negative=[[float(np.exp(ln[i])), float(cc[i])] for i in neg[:4]])
        rows = []
        for s in S:
            i0, i1 = np.searchsorted(ln, s - L), np.searchsorted(ln, s + L)
            u = ln[i0:i1]; cw = c[i0:i1] if wfun is None else c[i0:i1]*wfun(u)
            P = float(np.sum(cw*np.exp(-u/2)*X(u - s)))
            Pi = gl(lambda v: X(v)*kp(v + s), -L, 0) + gl(lambda v: X(v)*kp(v + s), 0, L)
            A = -(gl(lambda t: w(t)*X(t - s), s - L, s) + gl(lambda t: w(t)*X(t - s), s, s + L))
            rows.append((s, P, Pi, A, A + Pi - P))
        R = np.array(rows); s, P, Pi, A, b = R.T
        e = np.exp(s/2); sel = s >= 2.0
        # growth exponent of the Hasse off-diagonal: slope of log(running max |b|) over the last half of the range
        rm = np.maximum.accumulate(np.abs(b)); half = s >= s[0] + 0.5*(s[-1]-s[0])
        slope = float(np.polyfit(s[half], np.log(rm[half]), 1)[0])
        out[nm] = dict(smax=float(smax), d=diag[nm], **pos,
                       sup_P_over_es2=float(np.max(np.abs(P[sel])/e[sel])), P_over_es2_at_smax=float(P[-1]/e[-1]),
                       sup_abs_b=float(np.max(np.abs(b[sel]))), s_at_sup_b=float(s[sel][np.argmax(np.abs(b[sel]))]),
                       sup_abs_b_over_d=float(np.max(np.abs(b[sel]))/diag[nm]),
                       first_s_violation=(float(s[sel][np.argmax(np.abs(b[sel]) > diag[nm])]) if np.any(np.abs(b[sel]) > diag[nm]) else None),
                       growth_slope_log_runmax_b=slope,
                       table=[dict(s=float(x), P=float(p_), Pi=float(pi_), arch=float(a_), b=float(bb), P_over_es2=float(p_/math.exp(x/2)))
                              for x, p_, pi_, a_, bb in R[::max(1, len(R)//60)]])
        o = out[nm]; print(nm, {k: o[k] for k in o if k != 'table'})
    # cross-check vs Round 15 hasse_core at s <= 8 (model data at NX = e^8.6 there)
    chk = []
    for nm in ('ZETA', 'EPS', 'DH', 'FTH', 'FTAU', 'BEUR2'):
        mdl = hc.MODELS[nm]()
        for s0 in (1.5, 4.0, 7.5):
            ref = float(hc.Bs(mdl, L, 1, s0)['B'][0, 0])
            ln, c, wfun = D[nm]; i0, i1 = np.searchsorted(ln, s0 - L), np.searchsorted(ln, s0 + L); u = ln[i0:i1]
            cw = c[i0:i1] if wfun is None else c[i0:i1]*wfun(u)
            P = float(np.sum(cw*np.exp(-u/2)*X(u - s0))); kp = KP[nm]; w = W_of(nm)
            Pi = gl(lambda v: X(v)*kp(v + s0), -L, 0) + gl(lambda v: X(v)*kp(v + s0), 0, L)
            A = -(gl(lambda t: w(t)*X(t - s0), s0 - L, s0) + gl(lambda t: w(t)*X(t - s0), s0, s0 + L))
            chk.append(dict(model=nm, s=s0, b_round16=A + Pi - P, b_round15=ref, absdiff=abs(A + Pi - P - ref)))
    print('max |b16 - b15| =', max(r['absdiff'] for r in chk))
    out['crosscheck_vs_round15'] = chk
    json.dump(out, open('prime_window.json', 'w'), indent=1)

if __name__ == '__main__': main()
