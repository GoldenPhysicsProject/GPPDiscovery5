"""Round 32: the fixed-level Selberg mollified floor (P(x) = x, i.e. Selberg's Lambda_y):
   Lambda_y = b_M * (mu_M log(y/.) 1_{<=y}),  S_y(L) = sum_{1<n<=X} Lambda_y(n) n^{-1/2}(L - log n)
            = sum_{d<=min(y,X)} mu_M(d) log(y/d) d^{-1/2} A_M(X/d) - L log y,
   D_y(L) = P_M(L) - S_y(L)   versus the true floor D(L) = P_M(L) - S(L).
PROVED: D_y = D for X <= y (Lambda_y = Lambda on 1 < n <= y).  For X > y, generating function b_M(s) N_y(s) has NO poles
from zeros (N_y is a finite Dirichlet polynomial), so D_y is an explicit elementary function of finitely many mu_M(d), d<=y.
N = 1e7; models ZETA, F_0.25, DH, EPS; y = 1e2, 1e3, 1e4, 1e5."""
import json, math, numpy as np
from r22_models import KAP
from r32_common import P_zeta, P_model
N = 10**7; YS = (100, 1000, 10000, 100000)
pz = json.load(open('r32_pass_1e09.json')); pc = json.load(open('r32_controls_1e07.json'))
n = np.arange(N + 1); nf = np.maximum(n, 1).astype(float)

def b_of(name):
    if name == 'ZETA': b = np.ones(N + 1)
    elif name == 'FTH':
        b = np.zeros(N + 1)
        for a in range(1, N + 1):
            b[a::a] += (n[a::a]/(a*a))**0.25
    elif name == 'DH': b = np.array([[0, 1, KAP, -KAP, -1][k % 5] for k in range(N + 1)], float)
    else:
        b = np.zeros(N + 1); M = int(math.isqrt(N)) + 1
        for x in range(-M, M + 1):
            y = np.arange(-M, M + 1); v = x*x + 5*y*y; v = v[(v > 0) & (v <= N)]; np.add.at(b, v, 1)
        b /= 2
    b[0] = 0; return b
def inverse(b, M):
    binv = np.zeros(M + 1); binv[1] = 1.0; acc = np.zeros(M + 1)
    for d in range(1, M + 1):
        if d > 1: binv[d] = -acc[d]
        if binv[d] != 0 and 2*d <= M: acc[2*d::d] += binv[d]*b[2:M//d + 1]
    return binv

Lg = np.arange(2.0, math.log(N) - 1e-9, 0.02)
out = {}; curves = {}
for name, src, pre in (('ZETA', pz, 'lam'), ('FTH', pz, 'lf0.25'), ('DH', pc, 'DH'), ('EPS', pc, 'EPS')):
    b = b_of(name); mu = inverse(b, max(YS))
    C0 = np.cumsum(b*nf**-0.5); C1 = np.cumsum(b*nf**-0.5*np.log(nf))
    T = np.array(src['T']); P0 = np.array(src['cum'][f'{pre}_P0']); P1 = np.array(src['cum'][f'{pre}_P1'])
    ig = np.searchsorted(T, np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64))
    # thresholds of the L-grid step 0.02 are not all in T (step 0.005 grid): use the nearest lower threshold <= X
    ig = np.searchsorted(T, np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64), side='right') - 1
    S = Lg*P0[ig] - P1[ig]        # exact up to primes in (T[ig], X]: none if T[ig] == floor(X)
    exact = T[ig] == np.floor(np.exp(Lg)*(1 + 1e-13)).astype(np.int64)
    Pp = np.array([P_zeta(L) if name == 'ZETA' else P_model(name, L, 0.25 if name == 'FTH' else None) for L in Lg])
    D = Pp - S
    out[name] = {}; curves[name] = dict(L=Lg.tolist(), D=D.tolist(), exact=exact.tolist())
    for y in YS:
        lam_d = mu[1:y + 1]*np.log(y/np.arange(1, y + 1)); dd = np.arange(1, y + 1).astype(float)
        Sy = np.empty_like(Lg)
        for i, L in enumerate(Lg):
            X = math.exp(L); m = min(y, int(X)); Y = X/dd[:m]; Yi = np.floor(Y).astype(np.int64)
            Sy[i] = float(np.sum(lam_d[:m]*dd[:m]**-0.5*(np.log(Y)*C0[Yi] - C1[Yi]))) - L*math.log(y)
        Dy = Pp - Sy
        Ny1 = float(np.sum(lam_d/dd))                 # N_y(1): the residue weight of b*N_y at the main pole (ZETA)
        below = (Lg <= math.log(y) - 1e-9) & exact
        out[name][str(y)] = dict(N_y1=Ny1, maxdiff_below_logy=float(np.max(np.abs(Dy - D)[below])) if below.any() else None,
                                 Dy_end=float(Dy[-1]), D_end=float(D[-1]), Dy_scaled_end=float(Dy[-1]*math.exp(-Lg[-1]/2)),
                                 pred_scaled_ZETA=4*(1 - Ny1), max_abs_mu_le_y=float(np.max(np.abs(mu[1:y + 1]))),
                                 sum_abs_lambda_over_sqrt=float(np.sum(np.abs(lam_d)/np.sqrt(dd))))
        curves[name][f'Dy_{y}'] = Dy.tolist()
        print(name, y, out[name][str(y)], flush=True)
    del b
json.dump(dict(results=out, curves=curves), open('r32_selberg.json', 'w'))
