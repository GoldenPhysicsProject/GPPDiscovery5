"""Queen S1b: Kronecker/Bohr torus surrogate for DH zeros with Re s > 1 (none is explicitly computable: they sit at astronomically large heights).
For the window vector f = e^{-i t x} sinh(b x) on [-L, L], the prime part of Q_D(f) depends on t only through the torus point (p^{-it})_p, p < e^{2L}.
By Kronecker's theorem every torus point is a limit of actual heights t. We maximise
   P(phi) = 2 sum_{n<e^{2L}} c(n) n^{-1/2} K(log n) cos(Phi(n)),  Phi(n) = sum_p v_p(n) phi_p,
by coordinate ascent (64-point grid per prime, 6 sweeps, 4 random restarts). K(u) = <s, tau_u s>, s = sinh(b x) 1_[-L,L].
Reported: P*/||f||^2 (a lower bound for the sup over the torus) for DH (c = c_DH) and for zeta (c = Lambda), and the trivial mass bound
M = 2 sum |c(n)| n^{-1/2}|K(log n)|/||f||^2. Q_D(f)/||f||^2 = Omega_D(t) + O(1) - P, with Omega_D(t) ~ log t. NUMERICAL."""
import numpy as np, math, json, sys
N = 60000
xi = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
per = {0: 0.0, 1: 1.0, 2: xi, 3: -xi, 4: -1.0}
a = np.array([per[n % 5] for n in range(N+1)]); a[0] = 0
S = np.zeros(N+1); c = np.zeros(N+1)
for d in range(2, N+1):
    c[d] = a[d]*math.log(d) - S[d]
    if c[d] != 0.0 and 2*d <= N:
        m = N//d; S[2*d::d] += c[d]*a[2:m+1]
lam = np.zeros(N+1); sv = np.ones(N+1, bool); sv[:2] = False; primes = []
for p in range(2, N+1):
    if sv[p]:
        sv[2*p::p] = False; primes.append(p); q = p
        while q <= N: lam[q] = math.log(p); q *= p
# valuations
def vals(X, ps):
    V = {}
    for p in ps:
        ks = np.arange(p, X+1, p); v = np.zeros(len(ks), int); r = ks.copy()
        while True:
            m = r % p == 0
            if not m.any(): break
            v[m] += 1; r[m] //= p
        V[p] = (ks, v)
    return V
rng = np.random.default_rng(7)
grid = np.linspace(0, 2*np.pi, 64, endpoint=False)
def optimise(coef, L, b, X, V, ps):
    n = np.arange(X+1); lk = np.log(np.maximum(n, 1))
    K = np.where(n >= 2, np.sinh(b*(2*L-lk))/(2*b) - (2*L-lk)*np.cosh(b*lk)/2, 0.0)
    K[(n < 2) | (lk >= 2*L)] = 0.0
    nf = math.sinh(2*b*L)/(2*b) - L
    w = 2*coef[:X+1]*np.where(n >= 1, n, 1)**-0.5*K
    mass = float(np.sum(np.abs(w)))/nf
    best = -1e300
    for rs in range(4):
        phi = {p: (rng.uniform(0, 2*np.pi) if rs else 0.0) for p in ps}
        Phi = np.zeros(X+1)
        for p in ps:
            ks, v = V[p]; Phi[ks] += v*phi[p]
        for sweep in range(6):
            for p in ps:
                ks, v = V[p]
                base = Phi[ks] - v*phi[p]
                vals_ = (w[ks][None, :]*np.cos(base[None, :] + np.outer(grid, v))).sum(1)
                j = int(np.argmax(vals_)); dphi = grid[j]
                Phi[ks] = base + v*dphi; phi[p] = dphi
        val = float(np.sum(w*np.cos(Phi)))/nf
        best = max(best, val)
    return best, mass
out = []
for b in [0.6]:
    for L in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]:
        X = min(int(math.exp(2*L)), N); ps = [p for p in primes if p <= X]; V = vals(X, ps)
        dh, mdh = optimise(c, L, b, X, V, ps); ze, mze = optimise(lam, L, b, X, V, ps)
        row = dict(b=b, L=L, X=X, DH_Pstar=dh, DH_mass=mdh, ZETA_Pstar=ze, ZETA_mass=mze, eL=math.exp(L))
        out.append(row); print(row, flush=True)
json.dump(out, open("dh_torus.json", "w"), indent=1)
