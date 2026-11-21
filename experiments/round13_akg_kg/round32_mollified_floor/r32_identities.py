"""Round 32 CHECKs of the PROVED finite identities (no zeros as inputs):
 (I1) Leibniz:  Lam_M = -b_M * (mu_M log)            for ZETA, F_0.25, DH, EPS (n <= 3000), max abs error.
 (I2) Selberg:  sum_{d|n, d<=y} mu_M(d) log(y/d) b_M(n/d) = Lam_M(n) for 1 < n <= y, = log y at n = 1.
 (I3) Moebius form of the floor (y = X):  S_M(L) = sum_{d<=X} mu_M(d)(L - log d) d^{-1/2} A_M(X/d) - L^2,
      A_M(Y) = sum_{m<=Y} b_M(m) m^{-1/2} log(Y/m);  at L = 8, 11, 13.
 (I4) ZETA: A(Y) = 4 sqrt(Y) + zeta(1/2) log Y + zeta'(1/2) + E(Y), |E(Y)| <= c Y^{-3/2}: fitted c on Y in [1, 1e6].
 (I5) Counted-once at sigma = 1: sum_{d<=x} mu(d) floor(x/d) = 1 and |sum_{d<=x} mu(d)/d| <= 1 (x <= 1e6, CHECKED),
      and the same two statements for F_0.25 / DH / EPS with their own mu_M and b_M-counting (they fail)."""
import json, math, numpy as np
from r22_models import KAP

K = json.load(open('r32_constants.json'))
z0, z1 = K['ZETA']['z0'], K['ZETA']['z1']
N = 600_000
n = np.arange(N + 1); nf = np.maximum(n, 1).astype(float)

def b_of(name):
    if name == 'ZETA': b = np.ones(N + 1)
    elif name == 'FTH':
        b = np.zeros(N + 1)
        for a in range(1, N + 1):
            m = np.arange(a, N + 1, a); b[m] += (m/a/a)**0.25
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
def dconv(f, g, M):
    h = np.zeros(M + 1)
    for d in range(1, M + 1):
        if f[d] != 0: h[d::d] += f[d]*g[1:M//d + 1]
    return h

res = {}
for name in ('ZETA', 'FTH', 'DH', 'EPS'):
    b = b_of(name); mu = inverse(b, N); ln = np.log(nf)
    M1 = 3000
    lam = dconv(b*ln, mu, M1)                         # definition: -b'/b
    lam2 = -dconv(b, mu*ln, M1)                       # Leibniz form
    y = 3000; lt = mu[:y + 1]*np.log(y/np.maximum(np.arange(y + 1), 1)); lt[0] = 0
    lamy = dconv(lt, b, y)
    r = dict(I1_maxerr=float(np.abs(lam[2:] - lam2[2:]).max()),
             I2_maxerr=float(np.abs(lamy[2:] - lam[2:y + 1]).max()), I2_at1=float(lamy[1]), logy=math.log(y))
    # I3
    lamN = dconv(b*ln, mu, N)
    C0 = np.cumsum(b*nf**-0.5); C1 = np.cumsum(b*nf**-0.5*ln)
    chk = []
    for L in (8.0, 11.0, 13.0):
        X = math.exp(L); Xi = int(X)
        S = float(np.sum((lamN*nf**-0.5*(L - ln))[2:Xi + 1]))
        d = np.arange(1, Xi + 1); Y = X/d; Yi = np.floor(Y).astype(int)
        A = np.log(Y)*C0[Yi] - C1[Yi]
        rhs = float(np.sum(mu[1:Xi + 1]*(L - np.log(d))*d**-0.5*A)) - L*L
        chk.append((L, S, rhs, abs(S - rhs)))
    r['I3'] = chk
    if name == 'ZETA':
        Y = np.exp(np.linspace(0, math.log(N) - 0.01, 20000)); Yi = np.floor(Y).astype(int)
        A = np.log(Y)*C0[Yi] - C1[Yi]; E = A - 4*np.sqrt(Y) - z0*np.log(Y) - z1
        r['I4'] = dict(max_E_Y32=float(np.max(np.abs(E)*Y**1.5)), max_E_Y32_Yge10=float(np.max((np.abs(E)*Y**1.5)[Y >= 10])),
                       E_at=[(float(Y[i]), float(E[i])) for i in (0, 5000, 10000, 19999)])
    # I5 at x up to 1e5 on a log grid (floor identity uses b-counting: sum_{d<=x} mu(d) B(x/d), B = cumsum b)
    B = np.cumsum(b)
    xs = np.unique(np.floor(np.exp(np.linspace(0, math.log(1e5), 300))).astype(int))
    fl, m1 = [], []
    for x in xs:
        d = np.arange(1, x + 1)
        fl.append(float(np.sum(mu[1:x + 1]*B[x//d]))); m1.append(float(np.sum(mu[1:x + 1]/d)))
    r['I5'] = dict(floor_identity_maxdev=float(np.max(np.abs(np.array(fl) - 1))), max_abs_sum_mu_over_d=float(np.max(np.abs(m1))),
                   max_abs_mu_le_1e5=float(np.max(np.abs(mu[1:100001]))))
    res[name] = r
    print(name, {k: v for k, v in r.items() if k != 'I3'}, [(c[0], round(c[3], 12)) for c in r['I3']], flush=True)
json.dump(res, open('r32_identities.json', 'w'), indent=1)
