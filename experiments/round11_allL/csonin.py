"""C-Sonin first test (Turbo 6890778106, README sec 8.2). Float64 numerics, NUMERICAL label.
Space: g = sum_i c_i psi(. - x_i), x_i = uniform grid on [0,L] (n points), c in ker l_c (sum c_i cosh(x_i/2) = 0).
Forms (n x n, Toeplitz in lag): A = arch(x_i-x_j), P = primeS_L(x_i-x_j), Pol = pole(x_i-x_j),
N = Psi(x_i-x_j) = <psi(.-x_i), psi(.-x_j)>_{L2}  (Psi = psi*psi~, so N is the L2 Gram of g).
Reported, on ker l_c:
 (i)  literal: lam_min(A) and ||P||_op w.r.t. l2(c)            [basis dependent]
 (ii) basis-free: lam_min(A;N), ||P||_N, ||P+pole||_N w.r.t. ||g||_L2, generalized eigenproblems
 (iii) ratio rho = max |c.P.c| / c.A.c (needs A>0 on ker), and lam_min of total K=A-P+Pol (K = pole - primeS + arch) (reference)."""
import numpy as np, json, scipy.linalg as sl
from scipy import integrate, special
PsiHat = lambda t: (np.sinc(t/(2*np.pi)))**4
bI2 = 16*np.sinh(0.25)**2; PH_i2 = bI2**2
def triG(y):
    ay = abs(y)
    return 2/3 - y*y + ay**3/2 if ay <= 1 else ((2-ay)**3/6 if ay <= 2 else 0.0)
def pp(m):
    for p in range(2, m+1):
        if m % p == 0:
            k = m
            while k % p == 0: k //= p
            return p if k == 1 else None
def primeS(d, L):
    tot = 0.0
    for m in range(2, int(np.exp(2+L))+2):
        p = pp(m)
        if p: tot += np.log(p)/np.sqrt(m)*(triG(np.log(m)+d) + triG(-np.log(m)+d))
    return tot
w = lambda t: special.digamma(0.25 + 0.5j*t).real - np.log(np.pi)
def arch(d, T=1200.0):
    f = lambda t: PsiHat(t)*w(t)*np.cos(t*d)
    edges = np.arange(0, T+1e-9, 2*np.pi)
    return sum(integrate.quad(f, a, b, limit=200, epsabs=1e-14)[0] for a, b in zip(edges[:-1], edges[1:]))/np.pi
pole = lambda d: 2*PH_i2*np.cosh(d/2)
cache = {}
def ker(f, L, d):
    k = (f.__name__, L if f is primeS else None, round(d, 12))
    if k not in cache: cache[k] = f(d, L) if f is primeS else f(d)
    return cache[k]
def toe(f, L, xs):
    n = len(xs); return np.array([[ker(f, L, abs(xs[i]-xs[j])) for j in range(n)] for i in range(n)])
def gen_eigs(M, B):  # eigenvalues of M relative to B (B>0)
    return sl.eigh(M, B, eigvals_only=True)
rows = []
for L in (0.5, 1.0, 2.0):
    for n in (6, 11, 21, 41):
        xs = np.linspace(0, L, n)
        A = toe(arch, L, xs); P = toe(primeS, L, xs); Pol = toe(pole, L, xs)
        N = np.array([[triG(xs[i]-xs[j]) for j in range(n)] for i in range(n)])
        v = np.cosh(xs/2); Q = sl.null_space(v[None, :])            # orthonormal basis of ker l_c
        r = lambda M: Q.T@M@Q
        Ak, Pk, Polk, Nk = r(A), r(P), r(Pol), r(N)
        eA = np.linalg.eigvalsh(Ak); eP = np.linalg.eigvalsh(Pk)
        gA = gen_eigs(Ak, Nk); gP = gen_eigs(Pk, Nk); gPP = gen_eigs(Pk-Polk, Nk); gK = gen_eigs(Ak-Pk+Polk, Nk)
        row = dict(L=L, n=n, condN=float(np.linalg.cond(Nk)),
                   lit_lam_min_A=float(eA[0]), lit_norm_P=float(np.abs(eP).max()),
                   N_lam_min_A=float(gA[0]), N_lam_max_A=float(gA[-1]), N_norm_P=float(np.abs(gP).max()),
                   N_P_eigs_range=[float(gP[0]), float(gP[-1])], N_norm_P_minus_pole=float(np.abs(gPP).max()),
                   N_PminusPole_range=[float(gPP[0]), float(gPP[-1])], N_lam_min_total=float(gK[0]))
        if gA[0] > 0:
            rr = gen_eigs(Pk, Ak); rr2 = gen_eigs(Pk-Polk, Ak)
            row['rho_P_over_A'] = float(np.abs(rr).max()); row['rho_PminusPole_over_A'] = float(np.abs(rr2).max())
        rows.append(row); print(json.dumps(row), flush=True)
json.dump(dict(rows=rows, label='NUMERICAL (float64, arch quad to T=1200, tail<1.2e-9); no zeros used'),
          open('csonin_results.json', 'w'), indent=1)
