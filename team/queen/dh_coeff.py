"""Queen S1 sanity checks (NUMERICAL, float64 / mpmath). RH OPEN.
(1) Levy-Khintchine identity: Re psi(1/4+it/2) - psi(1/4) = int_0^inf (1-cos tu) e^{u/2}/sinh(u) du  (Gauss digamma integral).
(2) Davenport-Heilbronn: coefficients c(n) of -DH'/DH; support off prime powers, signs, growth of sum |c(n)| n^{-1/2} vs zeta.
(3) Riesz means R_L(t) = sum_{n<e^{2L}} c(n) n^{-1/2} cos(t log n) (1 - log n/(2L)), sup over t in [1,200]
    (= (1/2) <e_t, C_L e_t> for the window prime operator C_L and the test vector e_t = e^{itx}/sqrt(2L) on [-L,L]).
"""
import numpy as np, mpmath as mp, json, math
out = {}
mp.mp.dps = 30
nu = lambda u: mp.e**(u/2)/mp.sinh(u)
lk = []
for t in [0.5, 1, 5, 14.134725, 50]:
    lhs = mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2)) - mp.digamma(mp.mpf(1)/4)
    rhs = mp.quad(lambda u: (1 - mp.cos(t*u))*nu(u), mp.linspace(0, 60, 241) + [mp.inf])
    lk.append([t, float(lhs), float(rhs), float(abs(lhs-rhs))])
out["levy_khintchine"] = lk
out["Omega0"] = float(mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi))
print(lk, out["Omega0"], flush=True)
N = 400000
xi = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
per = {0: 0.0, 1: 1.0, 2: xi, 3: -xi, 4: -1.0}
a = np.array([per[n % 5] for n in range(N+1)]); a[0] = 0
S = np.zeros(N+1); c = np.zeros(N+1)
for d in range(2, N+1):
    c[d] = a[d]*math.log(d) - S[d]
    if c[d] != 0.0 and 2*d <= N:
        m = N//d
        S[2*d::d] += c[d]*a[2:m+1]
lam = np.zeros(N+1); sieve = np.ones(N+1, bool); sieve[:2] = False
for p in range(2, N+1):
    if sieve[p]:
        sieve[2*p::p] = False
        q = p
        while q <= N: lam[q] = math.log(p); q *= p
n = np.arange(N+1)
nz = np.nonzero(np.abs(c[1:]) > 1e-12)[0] + 1
out["DH_c_first"] = [[int(k), round(float(c[k]), 6)] for k in nz[:20]]
out["DH_frac_support_non_prime_power"] = float(np.mean(lam[nz] == 0))
out["DH_sign_neg_frac"] = float(np.mean(c[nz] < 0))
w = np.zeros(N+1); w[1:] = n[1:]**-0.5
cumDH = np.cumsum(np.abs(c)*w); cumZ = np.cumsum(lam*w)
out["mass_x_zeta_DH_2sqrtx"] = [[int(x), float(cumZ[x]), float(cumDH[x]), float(2*math.sqrt(x))] for x in [10, 100, 1000, 10**4, 10**5, N]]
print({k: out[k] for k in ("DH_c_first", "DH_frac_support_non_prime_power", "DH_sign_neg_frac", "mass_x_zeta_DH_2sqrtx")}, flush=True)
ts = np.linspace(1, 200, 8001)
rz = []
for L in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.4]:
    X = min(int(math.exp(2*L)), N)
    k = np.arange(2, X+1); lw = 1 - np.log(k)/(2*L); lw[lw < 0] = 0
    def sup(coef):
        v = coef[2:X+1]*k**-0.5*lw
        m = np.nonzero(v)[0]; vv = v[m]; lk_ = np.log(k[m])
        best, bt = -1e300, None
        for i in range(0, len(ts), 200):
            tt = ts[i:i+200]
            R = (vv[None, :]*np.cos(np.outer(tt, lk_))).sum(1)
            j = int(np.argmax(R))
            if R[j] > best: best, bt = float(R[j]), float(tt[j])
        return best, bt
    zs, zt = sup(lam); ds, dt = sup(c)
    rz.append([L, X, zs, zt, ds, dt, math.exp(L)])
    print(rz[-1], flush=True)
out["riesz_sup_[L,X,zeta_sup,t,DH_sup,t,e^L]"] = rz
json.dump(out, open("dh_coeff.json", "w"), indent=1)
