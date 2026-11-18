"""Round 14E2: 80-digit Weil-form Galerkin matrix W (same Dirichlet-sine basis and formula as round14b weil_galerkin.operators).
Port of Golden Eye's scratch mpw.py (AUDIT14 6.3/6.8), generalised to the F_theta control.  Arch integral: Gauss-Legendre,
`panels` x 3*2^(deg-1) nodes (deg 7 = 192); prime and pole terms closed form.  Integer / own-gamma-factor data only; no zeros.
Usage: python3 mpW.py MODEL L M [dps] [deg]   MODEL in {ZETA, FTH}  -> mp80/W_<MODEL>_<L>_<M>.pkl (mpf matrix as strings)."""
import sys, os, pickle, time, mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
HERE = os.path.dirname(os.path.abspath(__file__))
def vonm(N):
    out = []
    for n in range(2, N+1):
        m = n; p = None
        for q in range(2, int(n**0.5)+1):
            if m % q == 0: p = q; break
        if p is None: out.append((n, mp.log(n))); continue
        while m % p == 0: m //= p
        if m == 1: out.append((n, mp.log(p)))
    return out
def Amat(L, M, t):
    if t >= L: return mp.zeros(M, M)
    a = [(j+1)*mp.pi/L for j in range(M)]; Y = L - t; R = mp.zeros(M, M)
    for i in range(M):
        for j in range(M):
            ai, aj = a[i], a[j]; d = ai-aj; s = ai+aj
            T1 = (mp.sin(d*Y - aj*t) - mp.sin(-aj*t))/d if i != j else Y*mp.cos(aj*t)
            T2 = (mp.sin(s*Y + aj*t) - mp.sin(aj*t))/s
            R[i, j] = (T1 - T2)/L
    return R
def model(name, theta=mp.mpf('0.2')):
    if name == 'ZETA':
        return dict(C=mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi), w=lambda t: mp.exp(-t/2)/(1-mp.exp(-2*t)),
                    poles=[(mp.mpf(1)/2, 1)], c=lambda n, lam: lam)
    if name == 'FTH':   # F_theta = Lambda(s-theta)Lambda(s+theta): weights 2cosh(theta log n)Lambda(n), Gamma_R(s-th)Gamma_R(s+th)
        return dict(C=mp.digamma(mp.mpf(1)/4+theta/2) + mp.digamma(mp.mpf(1)/4-theta/2) - 2*mp.log(mp.pi),
                    w=lambda t: (mp.exp(-(mp.mpf(1)/2-theta)*t) + mp.exp(-(mp.mpf(1)/2+theta)*t))/(1-mp.exp(-2*t)),
                    poles=[(mp.mpf(1)/2+theta, 1), (mp.mpf(1)/2-theta, 1)], c=lambda n, lam: lam*2*mp.cosh(theta*mp.log(n)))
def W_model(name, L, M, deg=7, panels=4):
    mo = model(name); L = mp.mpf(L); I = mp.eye(M); w = mo['w']
    arch = (mo['C'] + 2*mp.quad(w, [L, mp.inf]))*I
    nodes = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
    for k in range(panels):
        a, b = L*k/panels, L*(k+1)/panels
        for x, wt in nodes:
            t = (b-a)/2*x + (b+a)/2; ww = (b-a)/2*wt
            A = Amat(L, M, t); arch += 2*ww*w(t)*(I - (A + A.T)/2)
    k_ = [(j+1)*mp.pi/L for j in range(M)]; nrm = mp.sqrt(2/L); pole = mp.zeros(M, M)
    for aa, wtp in mo['poles']:
        def Ib(b): return [(mp.exp(b*L)*(b*mp.sin(k*L) - k*mp.cos(k*L)) + k)/(b*b+k*k) for k in k_]
        ep = [mp.exp(-aa*L/2)*v for v in Ib(aa)]; em = [mp.exp(aa*L/2)*v for v in Ib(-aa)]
        c = [nrm*(x+y)/2 for x, y in zip(ep, em)]; s = [nrm*(x-y)/2 for x, y in zip(ep, em)]
        for i in range(M):
            for j in range(M): pole[i, j] += wtp*2*(c[i]*c[j] - s[i]*s[j])
    prime = mp.zeros(M, M)
    for n, lam in vonm(int(mp.exp(L))+1):
        if mp.log(n) < L:
            A = Amat(L, M, mp.log(n)); prime += mo['c'](n, lam)/mp.sqrt(n)*(A + A.T)
    W = arch + pole - prime
    return (W + W.T)/2
def path(name, L, M, deg=7): return os.path.join(HERE, 'mp80', f'W_{name}_{float(L)}_{M}' + ('' if deg == 7 else f'_deg{deg}') + '.pkl')
def load(name, L, M, deg=7):
    d = pickle.load(open(path(name, L, M, deg), 'rb')); mp.mp.dps = max(mp.mp.dps, d['dps'])
    return mp.matrix([[mp.mpf(x) for x in row] for row in d['W']])
if __name__ == '__main__':
    name, L, M = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    mp.mp.dps = int(sys.argv[4]) if len(sys.argv) > 4 else 80; deg = int(sys.argv[5]) if len(sys.argv) > 5 else 7
    t0 = time.time(); W = W_model(name, L, M, deg=deg)
    ev = sorted(mp.eigsy(W, eigvals_only=True))
    pickle.dump(dict(model=name, L=L, M=M, dps=mp.mp.dps, deg=deg, W=[[mp.nstr(W[i, j], mp.mp.dps+5) for j in range(M)] for i in range(M)]),
                open(path(name, L, M, deg), 'wb'))
    print(name, L, M, 'deg', deg, 'low', [mp.nstr(x, 8) for x in ev[:8]], 'max', mp.nstr(ev[-1], 6), 'sec', round(time.time()-t0), flush=True)
