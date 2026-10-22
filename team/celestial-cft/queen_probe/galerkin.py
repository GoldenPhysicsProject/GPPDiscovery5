"""lambda_min of the window Weil form Q_L for Epstein Z_Q (lead review item 2b). NUMERICAL (float64). Half-width L: window [-L,L].
Basis: orthonormal Legendre phi_k(x) = sqrt((2k+1)/2L) P_k(x/L), k < N. Real basis suffices: the kernel is real and symmetric, so
Q(f) = Q(Re f) + Q(Im f) and lambda_min over complex f equals lambda_min over real f.
Matrix (same Levy-Khintchine form as Queen's dh_sinh.py, with polarisation G -> S(u) = sym part of C_ij(u) = int phi_i(x) phi_j(x-u) dx):
  A = Omega0 I + int_0^{2L} (I - S(u)) nu(u) du + I int_{2L}^inf nu - 2 sum_k w_k S(u_k) + (p q^T + q p^T),
  where p_i = <phi_i, e^{-x/2}> and q_i = <phi_i, e^{x/2}>.
C_ij(u) is an exact polynomial integral on [-L+u, L]; it is computed by Gauss-Legendre with N+4 nodes, which is exact.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, math, json, numpy as np
from numpy.polynomial import legendre as Lg
sys.path.insert(0,'.')
from gen_coeffs import atoms

def phis(x, L, N):
    t = x/L; V = Lg.legvander(t, N-1)            # (len x, N)
    return (V*np.sqrt((2*np.arange(N)+1)/(2*L))).T  # (N, len x)

def Cmat(u, L, N, gx, gw):
    a, b = -L+u, L
    x = (b-a)/2*gx + (a+b)/2; w = (b-a)/2*gw
    A = phis(x, L, N); B = phis(x-u, L, N)
    C = (A*w) @ B.T
    return (C + C.T)/2

def matrix(L, N, us, ws, Om0):
    gx, gw = Lg.leggauss(N+4)
    I = np.eye(N)
    # archimedean: nu(u) = 1/sinh(u/2)
    K = 400; ux, uw = Lg.leggauss(K)
    arch = np.zeros((N,N))
    for j in range(8):                     # 8 panels on [0,2L], graded towards 0
        a = 2*L*(j/8)**2; b = 2*L*((j+1)/8)**2
        uu = (b-a)/2*ux+(a+b)/2; ww = (b-a)/2*uw
        for u, w in zip(uu, ww):
            arch += w*(I - Cmat(u, L, N, gx, gw))/math.sinh(u/2)
    tail = 2*math.log(1/math.tanh(2*L/4))   # int_{2L}^inf csch(u/2) du
    A = Om0*I + arch + tail*I
    m = us < 2*L
    for u, w in zip(us[m], ws[m]):
        A -= 2*w*Cmat(u, L, N, gx, gw)
    xx = L*gx; wx = L*gw; F = phis(xx, L, N)
    p = (F*wx) @ np.exp(-xx/2); q = (F*wx) @ np.exp(xx/2)
    A += np.outer(p, q) + np.outer(q, p)
    return A

if __name__ == '__main__':
    Q = tuple(int(v) for v in sys.argv[1].split(',')); Ls = [float(v) for v in sys.argv[2].split(',')]; Ns = [int(v) for v in sys.argv[3].split(',')]
    us, ws, Om0, mmin = atoms(Q, max(Ls))
    out = []
    for L in Ls:
        row = dict(form=Q, L=L, L_convention='half-width [-L,L]', lam_min={})
        for N in Ns:
            lam = np.linalg.eigvalsh(matrix(L, N, us, ws, Om0))[0]
            row['lam_min'][N] = float(lam)
            print(f"{Q} L={L} N={N}: lambda_min = {lam:+.6g}", flush=True)
        out.append(row)
    json.dump(out, open(f"galerkin_{Q[0]}{Q[1]}{Q[2]}.json", 'w'), indent=1)
