"""Queen S1: DH Riesz mean at the off-line zero rho0 = 0.808517 + 85.699348 i vs the Perron residue prediction
R_L(t) = Re sum_{n<x} c(n) n^{-1/2-it} (1 - log n/log x), x = e^{2L}; residue of -DH'/DH at rho0 is -1, so
R_L(gamma0) ~ Re[-x^{w}/(w^2 log x)], w = rho0 - 1/2 - i gamma0 = beta' (beta' = 0.308517). Also zeta sup over t in [10,200]. NUMERICAL."""
import numpy as np, math, json
N = 400000
xi = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
per = {0: 0.0, 1: 1.0, 2: xi, 3: -xi, 4: -1.0}
a = np.array([per[n % 5] for n in range(N+1)]); a[0] = 0
S = np.zeros(N+1); c = np.zeros(N+1)
for d in range(2, N+1):
    c[d] = a[d]*math.log(d) - S[d]
    if c[d] != 0.0 and 2*d <= N:
        m = N//d; S[2*d::d] += c[d]*a[2:m+1]
lam = np.zeros(N+1); sv = np.ones(N+1, bool); sv[:2] = False
for p in range(2, N+1):
    if sv[p]:
        sv[2*p::p] = False; q = p
        while q <= N: lam[q] = math.log(p); q *= p
g0, b = 85.699348, 0.808517 - 0.5
res = []
for L in [3.0, 4.0, 5.0, 5.5, 6.0, 6.4]:
    X = min(int(math.exp(2*L)), N); k = np.arange(2, X+1); lw = np.clip(1 - np.log(k)/(2*L), 0, None)
    v = c[2:X+1]*k**-0.5*lw; lk = np.log(k)
    tt = np.linspace(g0 - 3, g0 + 3, 1201)
    R = np.array([np.sum(v*np.cos(t*lk)) for t in tt])
    Rg = float(np.sum(v*np.cos(g0*lk)))
    pred = -math.exp(2*L*b)/(b*b*2*L)
    vz = lam[2:X+1]*k**-0.5*lw
    tz = np.linspace(10, 200, 4001); Rz = np.array([np.sum(vz*np.cos(t*lk)) for t in tz])
    row = dict(L=L, R_gamma0=Rg, pred_residue=pred, ratio=Rg/pred, min_near=float(R.min()), t_min=float(tt[R.argmin()]),
               max_near=float(R.max()), zeta_sup_10_200=float(Rz.max()), zeta_inf_10_200=float(Rz.min()))
    res.append(row); print(row, flush=True)
json.dump(res, open("dh_riesz2.json", "w"), indent=1)
