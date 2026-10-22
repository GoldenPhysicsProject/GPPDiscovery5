"""Exact time-domain assembly of the zeta Weil window form in a Legendre basis, both parities.
Goddess of Prime Number Theory (GPP), 2026-10-22. RH OPEN. Output is NUMERICAL (Rayleigh-Ritz upper bounds; midpoint eigen-solve).
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

Normalisation = Zhu arXiv:2608.24827 eq. (2)/(17) = team S14:  for real f on [-L,L], g = f * f~ (autocorrelation),
  Q(f) = pole - (gamma + log pi + log(1 - e^{-4L})) g(0) + int_0^{2L} 2 (e^{-2x} g(0) - e^{-x/2} g(x)) / (1 - e^{-2x}) dx
         - sum_{log n < 2L} 2 Lambda(n) n^{-1/2} g(log n),
  pole = +2 (int f cosh(x/2))^2 (even sector), -2 (int f sinh(x/2))^2 (odd sector)  [Zhu Lemma 6.1].
Basis phi_k(u) = sqrt((2k+1)/(2L)) P_k(u/L), k < 2N (N modes per parity). Orthonormal, so lambda = eig of the matrix.
g_jk(x) = int phi_j(u) phi_k(u-x) du is a polynomial in x on [0,2L]; integrated exactly by Gauss (M = 2N nodes).
The x-integral (integrand analytic on [0,2L] after the removable point x=0) uses P panels x nX Gauss nodes.
Usage: python legweil.py L N [prec_bits] [P] [extra]     L may be a fraction like 171/100
"""
import sys, json, time, math
from fractions import Fraction
from flint import arb, arb_mat, ctx
import mpmath as mp

def mang(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else 0
    return 0

def tompf(v):
    m, e = v.mid().man_exp()
    return mp.mpf((int(m), int(e)))

def gauss(n):
    out = [arb.legendre_p_root(n, k, weight=True) for k in range(n)]
    return [t for t, w in out], [w for t, w in out]

def run(Lq, N, prec, P=4, extra=60):
    ctx.prec = prec
    L = arb(Lq.numerator) / Lq.denominator
    D = 2 * N                       # degrees 0..D-1
    norms = [((2 * k + 1) / (2 * L)).sqrt() for k in range(D)]
    def basis(us):                  # rows k, cols points
        rows = [[arb(1)] * len(us)]
        ts = [u / L for u in us]
        rows.append(list(ts))
        for k in range(1, D - 1):
            a, b = rows[k], rows[k - 1]
            rows.append([((2 * k + 1) * t * a[i] - k * b[i]) / (k + 1) for i, t in enumerate(ts)])
        return [[norms[k] * v for v in rows[k]] for k in range(D)]
    tM, wM = gauss(D)               # exact for degree <= 2D-1 >= deg(phi_j phi_k) = 2D-2
    def Gsym(x):
        a = x - L; h = (L - a) / 2; c = (L + a) / 2
        us = [h * t + c for t in tM]; ws = [h * w for w in wM]
        Phi = basis(us); Psi = basis([u - x for u in us])
        A = arb_mat([[Phi[k][i] * ws[i] for i in range(D)] for k in range(D)])
        B = arb_mat(Psi).transpose()
        G = A * B
        return G, G.transpose()
    # x-integral nodes
    nX = 2 * D + extra              # polynomial degree in x <= 2D-1; Gauss nX exact to 2nX-1
    tX, wX = gauss(nX)
    arch = arb_mat(D, D)
    for pnl in range(P):
        lo = 2 * L * pnl / P; hi = 2 * L * (pnl + 1) / P
        h = (hi - lo) / 2; c = (hi + lo) / 2
        for t, w in zip(tX, wX):
            x = h * t + c
            G, Gt = Gsym(x)
            den = 1 - (-2 * x).exp()
            e2 = (-2 * x).exp(); eh = (-x / 2).exp()
            fac = 2 * h * w / den
            M = (G + Gt) * (-(eh * fac) / 2)
            for k in range(D):
                M[k, k] += e2 * fac
            arch += M
    const = -(arb.const_euler() + arb.pi().log() + (1 - (-4 * L).exp()).log())
    # primes
    prim = arb_mat(D, D)
    nmax = int(math.exp(2 * float(Fraction(Lq)))) + 1
    terms = []
    for n in range(2, nmax + 1):
        p = mang(n)
        if p == 0: continue
        ln = arb(n).log()
        if not (ln < 2 * L): continue
        terms.append(n)
        G, Gt = Gsym(ln)
        prim += (G + Gt) * (arb(p).log() / arb(n).sqrt())   # 2*Lambda/sqrt(n) * (G+G^T)/2
    # pole moments (entire integrand), Gauss on [-L, L]
    tP, wP = gauss(D + 40)
    us = [L * t for t in tP]; Phi = basis(us)
    ch = [(u / 2).cosh() * L * w for u, w in zip(us, wP)]
    sh = [(u / 2).sinh() * L * w for u, w in zip(us, wP)]
    cvec = [sum((Phi[k][i] * ch[i] for i in range(len(us))), arb(0)) for k in range(D)]
    svec = [sum((Phi[k][i] * sh[i] for i in range(len(us))), arb(0)) for k in range(D)]
    out = {"L": str(Lq), "N_per_parity": N, "prec_bits": prec, "panels": P, "nX": nX, "prime_powers": terms}
    for par in (0, 1):
        idx = [k for k in range(D) if k % 2 == par]
        mp.mp.prec = prec
        Qm = mp.matrix(N, N)
        for a, j in enumerate(idx):
            for b, k in enumerate(idx):
                v = arch[j, k] - prim[j, k]
                if j == k: v += const
                v += (2 * cvec[j] * cvec[k]) if par == 0 else (-2 * svec[j] * svec[k])
                Qm[a, b] = tompf(v)
        Qm = (Qm + Qm.T) / 2
        ev = mp.eigsy(Qm, eigvals_only=True)
        ev = sorted(ev)
        out["even" if par == 0 else "odd"] = {"lam1": mp.nstr(ev[0], 12), "lam2": mp.nstr(ev[1], 8), "lam3": mp.nstr(ev[2], 6)}
    # cross-parity leakage check (should vanish identically)
    out["max_cross_parity"] = mp.nstr(max(abs(tompf(arch[j, k] - prim[j, k])) for j in range(0, D, 2) for k in range(1, D, 2)), 3)
    out["max_radius"] = mp.nstr(max(tompf(arch[j, k].rad()) for j in range(D) for k in range(D)), 3)
    return out

if __name__ == "__main__":
    Lq = Fraction(sys.argv[1]); N = int(sys.argv[2])
    prec = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    P = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    extra = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    t0 = time.time(); o = run(Lq, N, prec, P, extra); o["sec"] = round(time.time() - t0, 1)
    print(json.dumps(o)); sys.stdout.flush()
