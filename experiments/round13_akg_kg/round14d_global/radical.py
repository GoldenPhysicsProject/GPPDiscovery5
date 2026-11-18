"""Round 14D face 3 (a): the Connes-Consani global map E_F(H)(x) = sum_n a_F(n) H(n x)  (u-form g(u) = e^{u/2} E_F H(e^u)),
built from each model's own integer data a_F(n) (zeta: 1; EPS: r_Q(n) = lattice count of x^2+5y^2; DH: periodic a(n);
F_theta: sum_{d|n} d^theta (n/d)^-theta; spliced Beurling: Beurling integers). Test: is E_F(H) in the radical of the
model's own Weil form, W_F(E_F H) = 0, and what do CROSS terms W_G(E_F H) give?  W evaluated from the prime side
on a grid (no zeros). H(x) = b(log x) p(log x), b a C^inf bump on [0.5, 2.5], p chosen so that the Mellin transform of H
vanishes at the model's poles (zeta, EPS, BEUR2: int H dx = 0; F_theta: at s = 1 +- theta; DH: none)."""
import sys, os, math, json, numpy as np
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'round14b_euler')); sys.path.insert(0, os.path.join(HERE, '..'))
import models14b as md
U0, U1, DU = -6.0, 2.5, 2e-3
u = np.arange(U0, U1+DU/2, DU); NMAXE = int(math.exp(U1 - U0)) + 2
def bump(v, a=0.5, b=2.5):
    y = (2*v - a - b)/(b - a); out = np.zeros_like(y, dtype=float); m = np.abs(y) < 1
    out[m] = np.exp(-1/(1 - y[m]**2)); return out
def make_H(exps):
    vv = np.linspace(0.5, 2.5, 20001); bb = bump(vv); k = len(exps)
    if k:
        A = np.array([[np.trapezoid(bb*vv**j*np.exp(s*vv), vv) for j in range(1, k+1)] for s in exps])
        rhs = -np.array([np.trapezoid(bb*np.exp(s*vv), vv) for s in exps]); c = np.linalg.solve(A, rhs)
    else: c = np.array([])
    return (lambda x: bump(np.log(x))*(1 + sum(cj*np.log(x)**(j+1) for j, cj in enumerate(c)))), c
def coeffs(name):
    n = np.arange(1, NMAXE+1).astype(float)
    if name == 'ZETA': return n, np.ones(len(n))
    if name == 'EPS':
        R = np.zeros(NMAXE+1); ymax = int(math.sqrt(NMAXE/5)) + 1
        for y in range(-ymax, ymax+1):
            xm = int(math.sqrt(max(NMAXE - 5*y*y, 0))); xs = np.arange(-xm, xm+1); q = xs**2 + 5*y*y
            q = q[(q > 0) & (q <= NMAXE)]; np.add.at(R, q, 1)
        return n, R[1:]
    if name == 'DH': cf = md.dh_coeff(); return n, np.array([cf[int(k) % 5] for k in n])
    if name == 'FTH':
        th = 0.2; b = np.zeros(NMAXE+1)
        for d in range(1, NMAXE+1): b[d::d] += d**th*(np.arange(d, NMAXE+1, d)/d)**(-th)
        return n, b[1:]
    if name == 'BEUR2':
        B = md.beur2(); pr = sorted(set(float(x) for x, c in zip(B['n'], B['c']) if abs(math.log(x) - c) < 1e-9))
        vals = [1.0]; stack = [(1.0, 0)]
        while stack:
            v, i = stack.pop()
            for j in range(i, len(pr)):
                w = v*pr[j]
                if w > NMAXE: break
                vals.append(w); stack.append((w, j))
        vals = np.array(sorted(vals)); return vals, np.ones(len(vals))
def E_u(H, name):
    n, a = coeffs(name); x = np.exp(u); G = np.zeros_like(u)
    for k, xx in enumerate(x):
        lo, hi = np.searchsorted(n, math.exp(0.5)/xx), np.searchsorted(n, math.exp(2.5)/xx, side='right')
        if hi > lo: G[k] = np.sum(a[lo:hi]*H(n[lo:hi]*xx))
    return G*np.exp(u/2)
def W_grid(model, g):
    N = len(g); A = np.correlate(g, g, 'full')[N-1:]*DU; t = np.arange(N)*DU
    parts = {}
    if model.get('C') is not None:
        w = model['w']; tt = np.maximum(t, 1e-12)
        integ = np.where(t > 0, w(tt)*(A[0] - A), 0.0)
        parts['arch'] = model['C']*A[0] + 2*np.trapezoid(integ, t) + 2*A[0]*quad(lambda s: float(w(s)), t[-1], np.inf)[0]
    else: parts['arch'] = 0.0
    parts['pole'] = sum(wt_*4*np.trapezoid(A*np.cosh(a*t), t) for a, wt_ in model.get('poles', []))
    m = np.log(model['n']) < t[-1]
    parts['prime'] = 2*sum(c/math.sqrt(nn)*np.interp(math.log(nn), t, A) for nn, c in zip(model['n'][m], model['c'][m]))
    parts['W'] = parts['arch'] + parts['pole'] - parts['prime']; parts['norm2'] = float(A[0])
    parts['rel'] = parts['W']/max(abs(parts['arch']), abs(parts['pole']), abs(parts['prime']))
    return {k: float(v) for k, v in parts.items()}
if __name__ == '__main__':
    models = {'ZETA': md.zeta(), 'EPS': md.eps(), 'DH': md.dh(), 'FTH': md.fth(0.2), 'BEUR2': md.beur2()}
    poles = {'ZETA': [1.0], 'EPS': [1.0], 'DH': [], 'FTH': [1.2, 0.8], 'BEUR2': [1.0]}
    out = {}
    for name in ['ZETA', 'EPS', 'DH', 'FTH', 'BEUR2']:
        H, c = make_H(poles[name]); g = E_u(H, name)
        out[name] = dict(c=[float(x) for x in c], g_left_edge=float(np.abs(g[:200]).max()), g_max=float(np.abs(g).max()))
        for wname in (['ZETA', name] if name != 'ZETA' else ['ZETA', 'EPS', 'BEUR2']):
            out[name][f'W_{wname}(E_{name}H)'] = W_grid(models[wname], g)
        gc = np.exp(u/2)*H(np.exp(u)); out[name]['control_W_own(H itself)'] = W_grid(models[name], gc)
        print(name, json.dumps(out[name]), flush=True)
    json.dump(out, open(os.path.join(HERE, 'radical.json'), 'w'), indent=1)
