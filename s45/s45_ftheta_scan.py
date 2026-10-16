"""S45 H45a: F_θ prime-weight deformation of ζ Weil form on S14 exact Arb assembly.
W_θ = Pole_ζ + Arch_ζ − Σ_n 2Λ(n)/√n · (n^θ + n^{-θ}) · F(log n)

GAP: arch+pole stay at ζ. CERTIFIED only for this cosh-weight deformation.
Full F_θ needs shifted gamma — NUMERICAL via weil_ritz / S4.
"""
import sys, json, time
sys.path.insert(0, '.')
from s14_exact import *
import mpmath as mp
from flint import arb, fmpq, ctx, arb_mat, acb

def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return (p if m == 1 else None)
    return None

def build_all(M, h, prec):
    ctx.prec = prec
    ha = arb(h)
    Gs = {}
    for j in range(M):
        for k in range(j, M):
            if (j + k) % 2 == 0:
                g = (G(j, k) + G(k, j)) / 2
                Gs[j, k] = Gs[k, j] = g
    D = max(g.degree() for g in Gs.values()) + 1
    def mom(phi):
        return [acb.integral(lambda z, an, m=m: z ** m * phi(z), 0, 2).real for m in range(D + 1)]
    hc = acb(ha)
    m1 = mom(lambda z: (-hc * z).exp() / (acb(0, 1) * hc * z).sinc())
    m2 = mom(lambda z: (-hc * z / 4).exp() * (acb(0, 1) * hc * z * 3 / 4).sinc() * fmpq(3, 4) / (acb(0, 1) * hc * z).sinc())
    m3 = mom(lambda z: (hc * z / 2).cosh())
    c0 = -arb.const_euler() - arb.pi().log() - (1 - (-4 * ha).exp()).log()
    ns = [n for n in range(2, int(mp.exp(2 * mp.mpf(h))) + 2)
          if vonmangoldt(n) and mp.log(n) < 2 * mp.mpf(h)]
    A = arb_mat(M, M)
    T = {n: arb_mat(M, M) for n in ns}
    for (j, k), g in Gs.items():
        co = [arb(c.p) / c.q for c in g.coeffs()]
        G0 = co[0]
        q = [-c for c in co[1:]]
        pole = 4 * ha * ha * sum(co[m] * m3[m] for m in range(len(co)))
        arch = (ha * c0 * G0 + ha * sum(q[m] * m1[m] for m in range(len(q)))
                - 2 * ha * ha * sum(co[m] * m2[m] for m in range(len(co))))
        A[j, k] = pole + arch
        for n in ns:
            v = arb(n).log() / ha
            T[n][j, k] = ha * sum(co[m] * v ** m for m in range(len(co)))
    return A, T, ns

def tomp(Ma, I):
    return mp.matrix([[mp.mpf(Ma[i, j].mid().str(mp.mp.dps + 20, radius=False)) for j in I] for i in I])

def weight_zeta(n):
    p = vonmangoldt(n)
    return 2 * mp.log(p) / mp.sqrt(n)

def weight_ftheta(n, theta):
    p = vonmangoldt(n)
    return 2 * mp.log(p) / mp.sqrt(n) * (mp.power(n, theta) + mp.power(n, -theta))

def run(M, h, theta, prec):
    t0 = time.time()
    A, T, ns = build_all(M, h, prec)
    maxrad = max(float(A[i, j].rad()) for i in range(M) for j in range(M))
    row = {"M": M, "h": h, "theta": str(theta), "prime_powers": ns, "maxrad_A": maxrad,
           "kind": "ftheta_prime_def", "arch": "zeta",
           "gap_note": "arch+pole stay at zeta; only prime weights (n^θ+n^{-θ})-deformed"}
    for label, wfun in (("zeta", weight_zeta),
                        ("ftheta_def", lambda n: weight_ftheta(n, theta))):
        for nm, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
            a = tomp(A, I)
            Tm = {n: tomp(T[n], I) for n in ns}
            W = a.copy()
            for n in ns:
                W -= wfun(n) * Tm[n]
            ev = mp.eigsy(W, eigvals_only=True)
            lmin = min(ev)
            row[f"{label}_{nm}"] = {"lam_min": mp.nstr(lmin, 12), "negative": bool(lmin < 0)}
    row["sec"] = round(time.time() - t0, 1)
    return row

if __name__ == "__main__":
    # argv: M hs_csv theta [dps] [prec]
    M = int(sys.argv[1])
    hs = [float(x) for x in sys.argv[2].split(',')]
    theta = mp.mpf(sys.argv[3])
    dps = int(sys.argv[4]) if len(sys.argv) > 4 else 110
    prec = int(sys.argv[5]) if len(sys.argv) > 5 else 768
    mp.mp.dps = dps
    for h in hs:
        print(json.dumps(run(M, h, theta, prec)), flush=True)
