"""Recompute the best verified separation witness (float). NUMERICAL. RH OPEN."""
import math
import numpy as np

def mangoldt_upto(N):
    Lam = np.zeros(N + 1)
    is_p = np.ones(N + 1, dtype=bool); is_p[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if is_p[p]:
            is_p[p*p::p] = False
    for p in range(2, N + 1):
        if is_p[p]:
            q = p; lp = math.log(p)
            while q <= N:
                Lam[q] = lp; q *= p
    return Lam

def yoshida(u, h, k=4):
    out = np.zeros_like(u, dtype=float)
    m = np.abs(u) < h
    out[m] = (1.0 - (u[m] / h) ** 2) ** k
    return out

def Q(g, u, Lam, theta=0.0):
    dens = np.exp(0.5 * u) + np.exp(-0.5 * u)
    if theta:
        dens = dens * (np.exp(theta * u) + np.exp(-theta * u))
    pole = float(np.trapezoid(g * dens, u))
    hmax = float(np.max(np.abs(u)))
    abs_u = np.abs(u); order = np.argsort(abs_u)
    uu, idx = np.unique(abs_u[order], return_index=True)
    gg = g[order][idx]
    prime = 0.0
    Nmax = min(int(math.floor(math.exp(hmax))), len(Lam) - 1)
    for n in range(2, Nmax + 1):
        lam = Lam[n]
        if not lam:
            continue
        ln = math.log(n)
        if ln > hmax:
            break
        gval = float(np.interp(ln, uu, gg, left=0.0, right=0.0))
        w = 1.0 if theta == 0.0 else (n**theta + n**(-theta))
        prime += 2.0 * lam * (n**-0.5) * w * gval
    return pole - prime, pole, prime

def main():
    c, w, k = 9.0, 0.4, 4
    h = c + w + 0.02
    Lam = mangoldt_upto(int(math.exp(h)) + 2)
    u = np.linspace(-h, h, max(8001, int(2000 * h) | 1))
    g = yoshida(u - c, w, k) + yoshida(u + c, w, k)
    qz, pz, prz = Q(g, u, Lam, 0.0)
    qt, pt, prt = Q(g, u, Lam, 0.25)
    print("witness: paired Yoshida c=9 w=0.4 k=4 (NOT positive-type)")
    print(f"Q_zeta  = {qz:+.10e}  (pole={pz:.6e}, prime={prz:.6e})")
    print(f"Q_1/4   = {qt:+.10e}  (pole={pt:.6e}, prime={prt:.6e})")
    print("status: NUMERICAL float; RH OPEN; off Weil class")

if __name__ == "__main__":
    main()
