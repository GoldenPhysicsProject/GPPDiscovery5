"""Grok S87: high-precision (mpmath, exact-polynomial) rebuild of the S86 parity window forms
A (arch+delta, no primes, no pole), P (prime powers), pole vectors, for ZETA / DH / F_theta.
Same normalisation as s86/h86_prime_path.py (S14). Basis: orthonormal shifted Legendre on [0,L].
Method: Hankel/Toeplitz convolution matrices are exact polynomials in sigma=s/L (monomial basis, exact
rational coefficients), contracted with kernel moments computed by mp Gauss-Legendre (analytic integrands),
then mapped to Legendre by the exact integer change of basis. Not interval arithmetic: HIGH-PRECISION NUMERICAL.
"""
import mpmath as mp
from math import comb, factorial
from mpmath.calculus.quadrature import GaussLegendre

def setup(dps):
    mp.mp.dps = dps
_GN = {}
def gnodes(deg=7):
    key = (deg, mp.mp.prec)
    if key not in _GN:
        nd = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
        _GN[key] = [((x + 1) / 2, w / 2) for x, w in nd]       # on [0,1]
    return _GN[key]

def legC(N):  # phi_i(u)=sqrt(2i+1) sum_k (-1)^{i+k} C(i,k)C(i+k,k) u^k
    C = [[mp.mpf(0)] * N for _ in range(N)]
    for i in range(N):
        s = mp.sqrt(2 * i + 1)
        for k in range(i + 1):
            C[i][k] = s * ((-1) ** (i + k) * comb(i, k) * comb(i + k, k))
    return C

def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[mp.fsum(A[i][k] * B[k][j] for k in range(m) if A[i][k] != 0) for j in range(p)] for i in range(n)]
def T(A): return [list(r) for r in zip(*A)]
def add(A, B, a=1, b=1): return [[a * x + b * y for x, y in zip(r, q)] for r, q in zip(A, B)]
def sym(A): return add(A, T(A), mp.mpf(1) / 2, mp.mpf(1) / 2)
def CMCt(C, M): return matmul(matmul(C, M), T(C))

def mang(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else mp.mpf(0)
    return mp.mpf(0)
def dh_weights(NN):
    s5 = mp.sqrt(5); kap = (mp.sqrt(10 - 2 * s5) - 2) / (s5 - 1)
    chi = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0)}
    aa = {n: chi[n % 5][0] + kap * chi[n % 5][1] for n in range(1, NN + 1)}
    bb = {1: mp.mpf(0)}
    for n in range(2, NN + 1):
        bb[n] = aa[n] * mp.log(n) - mp.fsum(bb[d] * aa[n // d] for d in range(2, n) if n % d == 0)
    return bb

def build(FORM, L, N, theta=0, deg=7):
    L = mp.mpf(L); theta = mp.mpf(theta); C = legC(N); G = gnodes(deg)
    D = [(-1) ** i for i in range(N)]
    wt = lambda r: mp.cosh(theta * r)
    k1 = lambda z: mp.exp(-z) / mp.sinh(z)
    k2 = lambda z: 2 * mp.exp(-z / 4) * mp.sinh(3 * z / 4) / mp.sinh(z)
    sech = lambda z: 1 / mp.cosh(z / 2)
    dh = FORM == 'DH'
    ker = lambda s: (-(k1(s) + k2(s)) / 2 + (sech(s) / 2 if dh else 0)) * wt(s)
    K = 2 * N
    # Hankel: h^(1)_kl(tau)=k!l!/(k+l+1)! tau^{k+l+1}; piece2 h_ij(sigma)=(-1)^{i+j} h^(1)_ij(2-sigma)
    fac = [mp.mpf(factorial(n)) for n in range(2 * N + 2)]
    H1 = lambda mu: [[fac[k] * fac[l] / fac[k + l + 1] * mu[k + l + 1] for l in range(N)] for k in range(N)]
    mu1 = [mp.fsum(w * ker(L * x) * x ** n for x, w in G) for n in range(K + 1)]
    mu2 = [mp.fsum(w * ker(L * (2 - x)) * x ** n for x, w in G) for n in range(K + 1)]
    B1 = CMCt(C, H1(mu1)); B2 = CMCt(C, H1(mu2))
    Barch = [[L * (B1[i][j] + D[i] * D[j] * B2[i][j]) for j in range(N)] for i in range(N)]
    # Toeplitz: e_kl(rho)=sum_a C(k,a)/(a+l+1) rho^{k-a}(1-rho)^{a+l+1}
    def nu(f, delta0):  # nu[p][q]=int_0^1 f(rho)[delta_{p0}*delta0 - rho^p (1-rho)^q] (delta0 path for k1)
        vals = [(x, w, f(x)) for x, w in G]
        out = [[None] * (K + 2) for _ in range(N)]
        for p in range(N):
            for q in range(K + 2):
                if p == 0 and delta0:
                    out[p][q] = mp.fsum(w * fx * (1 - wtL(x) * (1 - x) ** q) for x, w, fx in vals)
                else:
                    out[p][q] = mp.fsum(w * fx * x ** p * (1 - x) ** q for x, w, fx in vals) * (-1 if delta0 else 1) * (1 if not delta0 else 1)
        return out
    wtL = lambda x: wt(L * x)
    # integrand pieces: I1 = int k1(L r)[e(0) - wt e(r)],  I2 = int (-wt k2 + dh*wt*sech) e(r)
    v1 = [(x, w, k1(L * x)) for x, w in G]
    v2 = [(x, w, wt(L * x) * (-k2(L * x) + (sech(L * x) if dh else 0))) for x, w in G]
    def Tmono(pq):  # contract nu-like table pq[p][q] -> monomial matrix
        return [[mp.fsum(comb(k, a) * pq[k - a][a + l + 1] / (a + l + 1) for a in range(k + 1)) for l in range(N)] for k in range(N)]
    pq1 = [[None] * (K + 2) for _ in range(N)]; pq2 = [[None] * (K + 2) for _ in range(N)]
    for p in range(N):
        for q in range(K + 2):
            if p == 0:
                pq1[p][q] = mp.fsum(w * f * (1 - wtL(x) * (1 - x) ** q) for x, w, f in v1)
            else:
                pq1[p][q] = -mp.fsum(w * f * wtL(x) * x ** p * (1 - x) ** q for x, w, f in v1)
            pq2[p][q] = mp.fsum(w * f * x ** p * (1 - x) ** q for x, w, f in v2)
    Tint = CMCt(C, add(Tmono(pq1), Tmono(pq2)))
    c0 = -mp.euler - mp.log(mp.pi) - mp.log(1 - mp.exp(-2 * L)) + (mp.log(5) if dh else 0)
    Tarch = [[L * Tint[i][j] + (c0 if i == j else 0) for j in range(N)] for i in range(N)]
    Tarch = sym(Tarch); Barch = sym(Barch)
    # primes
    NN = int(mp.floor(mp.exp(2 * L))) + 1
    wts = dh_weights(NN) if dh else {n: mang(n) for n in range(2, NN + 1)}
    Mh = [[mp.mpf(0)] * N for _ in range(N)]; Mh2 = [[mp.mpf(0)] * N for _ in range(N)]; Mt = [[mp.mpf(0)] * N for _ in range(N)]
    for n in range(2, NN + 1):
        if wts.get(n, 0) == 0: continue
        ln = mp.log(n); c = wts[n] / mp.sqrt(n) * wt(ln); sg = ln / L
        if sg < 2:
            if sg <= 1:
                for k in range(N):
                    for l in range(N): Mh[k][l] += c * fac[k] * fac[l] / fac[k + l + 1] * sg ** (k + l + 1)
            else:
                t = 2 - sg
                for k in range(N):
                    for l in range(N): Mh2[k][l] += c * fac[k] * fac[l] / fac[k + l + 1] * t ** (k + l + 1)
        if sg < 1:
            for k in range(N):
                for l in range(N):
                    Mt[k][l] += 2 * c * mp.fsum(comb(k, a) / mp.mpf(a + l + 1) * sg ** (k - a) * (1 - sg) ** (a + l + 1) for a in range(k + 1))
    P1 = CMCt(C, Mh); P2 = CMCt(C, Mh2)
    Bp = sym([[P1[i][j] + D[i] * D[j] * P2[i][j] for j in range(N)] for i in range(N)])
    Tp = sym(CMCt(C, Mt))
    # poles: mom(al) = sqrt(L) int phi(u) e^{+-al L u}
    poles = []
    if FORM == 'ZETA': poles = [(mp.mpf(1) / 2, 2)]
    if FORM == 'FT': poles = [(mp.mpf(1) / 2 + theta, 1), (mp.mpf(1) / 2 - theta, 1)]
    pe, po = [], []
    for al, c in poles:
        ap = matmul(C, [[mp.sqrt(L) * mp.fsum(w * x ** k * mp.exp(al * L * x) for x, w in G)] for k in range(N)])
        am = matmul(C, [[mp.sqrt(L) * mp.fsum(w * x ** k * mp.exp(-al * L * x) for x, w in G)] for k in range(N)])
        s = mp.sqrt(mp.mpf(c) / 2)
        pe.append([s * (ap[i][0] + am[i][0]) for i in range(N)]); po.append([s * (ap[i][0] - am[i][0]) for i in range(N)])
    return dict(Ae=add(Tarch, Barch), Ao=add(Tarch, Barch, 1, -1), Pe=add(Tp, Bp), Po=add(Tp, Bp, 1, -1), pe=pe, po=po)

def eigs(M):
    E = mp.eigsy(mp.matrix(M), eigvals_only=True)
    return sorted(E)
