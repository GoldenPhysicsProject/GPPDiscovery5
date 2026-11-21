"""Round 31, task 2: the integer floor D(4 log k) >= -B on [K, ceil(K^{3/2})], four models, no zeros as inputs.
D_M(L) = sum_rho (1 - cos gamma L)/gamma^2 computed from the integers by the explicit formula (r22_fejer.rho at gamma0=0):
   D_M(L) = (pi L/2) * [pole + archimedean part of rho at X=e^L]  -  sum_{n<=e^L} Lam_M(n) n^{-1/2} (L - log n).
For ZETA we also evaluate GPPVerify's closed form (GppHeartLaplace.Dfun, main c81c031):
   Dfun(L) = 4(e^{L/2} - 2 + e^{-L/2}) + c1 L + C2 - sum_k e^{-b_k L}/b_k^2 - sum_n Lam(n) n^{-1/2} max(L - log n, 0),
   b_k = 2k + 1/2, c1 = (Re digamma(1/4) - log pi)/2, C2 = sum_k 1/b_k^2   -> CHECK that the two agree.
Prime-side sums use cumulative arrays over prime powers, so every sample k costs O(log) after one sieve."""
import json, math, sys, numpy as np, mpmath
from r22_fejer import rho
from r22_models import MODELS, lam_table, TH

def prime_powers(N):
    s = np.ones(N + 1, np.uint8); s[:2] = 0
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = 0
    P = np.nonzero(s)[0].astype(np.int64); del s
    ns, lg = [P], [np.log(P.astype(float))]
    q = P[P <= int(N**0.5)]
    pw = q*q; base = q.copy()
    while len(pw):
        keep = pw <= N; pw, base = pw[keep], base[keep]
        if not len(pw): break
        ns.append(pw.copy()); lg.append(np.log(base.astype(float))); pw = pw*base
    n = np.concatenate(ns); l = np.concatenate(lg); o = np.argsort(n)
    return n[o], l[o]            # n, Lambda(n)

class PrimeSide:
    """cumulative A(X)=sum_{n<=X} c_n, B(X)=sum c_n log n with c_n = Lam_M(n) n^{-1/2}."""
    def __init__(self, n, lam):
        self.n = n.astype(float); c = lam*self.n**-0.5
        self.A = np.cumsum(c); self.B = np.cumsum(c*np.log(self.n))
    def S(self, L):
        i = np.searchsorted(self.n, math.exp(L)*(1 + 1e-12), side='right')
        return 0.0 if i == 0 else L*self.A[i-1] - self.B[i-1]

from r22_fejer import K as _K, arch_I
def nonprime_rho(name, X):
    """rho_X(0) of r22_fejer without the prime sum: log Q + archimedean + pole terms."""
    M = MODELS[name]; L = math.log(X); g0 = np.array([0.0])
    out = 2*(1/np.pi)*M['logQ']*np.ones_like(g0)
    for (lmb, a, b) in M['arch']:
        out += 2*lmb*arch_I(a, b, L, g0)
    for z in M['poles']:
        r = -1j*(z - 0.5)
        out += 2*(_K(r - g0, L) + _K(r + g0, L)).real
    return float(out[0])
def D_model(name, L, ps):
    return math.pi*L/2*nonprime_rho(name, math.exp(L)) - ps.S(L)

# closed form for ZETA (Lean)
c1 = float((mpmath.digamma(0.25) - mpmath.log(mpmath.pi))/2)
C2 = float(mpmath.nsum(lambda k: 1/(2*k + 0.5)**2, [0, mpmath.inf]))
def Tfun(L):
    k = np.arange(0, 2000); b = 2*k + 0.5
    return float(np.sum(np.exp(-b*L)/b**2))
def Dfun_lean(L, ps):
    return 4*(math.exp(L/2) - 2 + math.exp(-L/2)) + c1*L + C2 - Tfun(L) - ps.S(L)

if __name__ == '__main__':
    NZ = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**8
    NO = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**7
    out = dict(N_zeta_fth=NZ, N_eps_dh=NO, c1=c1, C2=C2)
    n, lam = prime_powers(NZ)
    PS = {'ZETA': PrimeSide(n, lam), 'FTH': PrimeSide(n, lam*(n.astype(float)**TH + n.astype(float)**-TH))}
    for name in ('EPS', 'DH'):
        lt = lam_table(name, NO + 1); idx = np.nonzero(lt[2:])[0] + 2
        PS[name] = PrimeSide(idx, lt[idx]); del lt
        print('lam', name, 'done', flush=True)
    # CHECK: Lean closed form vs explicit-formula route (ZETA), several L
    chk = []
    for L in (1.0, 3.0, 6.0, 9.0, 12.0, 15.0):
        a = Dfun_lean(L, PS['ZETA']); b = D_model('ZETA', L, PS['ZETA']); chk.append((L, a, b, abs(a - b)))
        print('check', L, a, b, flush=True)
    out['check_lean_vs_explicit'] = chk
    # integer samples
    res = {}
    for name in ('ZETA', 'FTH', 'EPS', 'DH'):
        N = NZ if name in ('ZETA', 'FTH') else NO
        kmax = int(math.floor(N**0.25))
        ks = np.arange(2, kmax + 1)
        Dk = np.array([Dfun_lean(4*math.log(k), PS[name]) if name == 'ZETA' else D_model(name, 4*math.log(k), PS[name]) for k in ks])
        # also a fine real-L grid to compare sampling vs continuum
        Lg = np.arange(0.05, math.log(N), 0.01)
        Dg = np.array([D_model(name, L, PS[name]) for L in Lg])
        win = []
        for K in range(2, kmax + 1):
            K2 = math.ceil(K**1.5)
            if K2 > kmax: break
            sel = (ks >= K) & (ks <= K2)
            d = Dk[sel]
            win.append(dict(K=K, K32=K2, min=float(d.min()), kmin=int(ks[sel][np.argmin(d)]),
                            frac_ge0=float(np.mean(d >= 0)), frac_ge_m01=float(np.mean(d >= -0.1))))
        res[name] = dict(k=ks.tolist(), D=Dk.tolist(), windows=win, Lgrid=Lg.tolist(), Dgrid=Dg.tolist(),
                         min_grid=float(Dg.min()), L_at_min=float(Lg[np.argmin(Dg)]))
        print(name, 'kmax', kmax, 'min D(4logk)', float(Dk.min()), 'at k', int(ks[np.argmin(Dk)]),
              'grid min', float(Dg.min()), 'at L', float(Lg[np.argmin(Dg)]), flush=True)
        for w in win[::max(1, len(win)//6)]: print('  ', w, flush=True)
    out['models'] = res
    json.dump(out, open('r31_floor.json', 'w'))
