"""v2: paired bumps, larger h, and compare renorm vs full Weil (with arch) on same g.
NUMERICAL float. RH OPEN.
"""
from __future__ import annotations
import json, math, os, time
import numpy as np
from scipy.special import digamma

def mangoldt_upto(N: int):
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

def paired_yoshida(u, c, w, k=4):
    """Even: bump at +c and -c, width w (support half-width of each bump)."""
    return yoshida(u - c, w, k) + yoshida(u + c, w, k)

def pole_zeta(u):
    return np.exp(0.5 * u) + np.exp(-0.5 * u)

def pole_theta(u, theta):
    return pole_zeta(u) * (np.exp(theta * u) + np.exp(-theta * u))

def Q_renorm(g, u, Lam, theta=0.0):
    dens = pole_zeta(u) if theta == 0.0 else pole_theta(u, theta)
    pole = float(np.trapezoid(g * dens, u))
    hmax = float(np.max(np.abs(u)))
    Nmax = int(math.floor(math.exp(hmax)))
    # interp g on [0,hmax]
    abs_u = np.abs(u)
    order = np.argsort(abs_u)
    uu, idx = np.unique(abs_u[order], return_index=True)
    gg = g[order][idx]
    prime = 0.0
    for n in range(2, min(Nmax, len(Lam) - 1) + 1):
        lam = Lam[n]
        if lam == 0.0:
            continue
        ln = math.log(n)
        if ln > hmax + 1e-15:
            break
        gval = float(np.interp(ln, uu, gg, left=0.0, right=0.0))
        w = 1.0 if theta == 0.0 else (n**theta + n**(-theta))
        prime += 2.0 * lam * (n**-0.5) * w * gval
    return pole - prime, pole, prime

def arch_term(g, u, theta=0.0):
    """Archimedean contribution matching weil_ritz convention for the LINEAR
    functional on even g (not the quadratic Ritz form).

    For the explicit formula applied to even g itself:
      Arch = (1/2pi) int_{-inf}^{inf} hat(g)(r) Re(psi(1/4+ir/2)-log pi) dr
    with hat(g)(r) = int g(u) e^{-i r u} du  (cosine transform since even).

    For F_theta: sum over shifts? Actually for F_theta the completed factor is
    product of two xi, each with same Gamma; the arch weight for zeros of F_theta
    is more subtle. We use the same digamma weight as zeta for both (degree-2
    for F_theta would double — check).

    Standard: for Xi(s)=xi(s-theta)xi(s+theta), log Xi' has
    arch = digamma from each Gamma factor:
      Gamma((s-theta)/2) and Gamma((s+theta)/2), so
      Re[psi((1/4 - theta/2)+ir/2) + psi((1/4+theta/2)+ir/2)]/2 - log pi
    (roughly). We implement that for theta!=0.
    """
    # cosine transform of g on dense r-grid
    # support of g <= H
    H = float(np.max(np.abs(u)))
    R = 200.0 / max(H, 0.1)
    r = np.linspace(0.0, R, 8001)
    dr = r[1] - r[0]
    # hat_g(r) = 2 int_0^H g(u) cos(r u) du
    # build via trapezoid matrix-free
    u_pos = u[u >= 0]
    g_pos = g[u >= 0]
    # ensure u=0 once
    if u_pos.size == 0:
        return 0.0
    cg = np.zeros_like(r)
    for i, rr in enumerate(r):
        cg[i] = 2.0 * float(np.trapezoid(g_pos * np.cos(rr * u_pos), u_pos))
    if theta == 0.0:
        w = digamma(0.25 + 0.5j * r).real - math.log(math.pi)
    else:
        w = (digamma(0.25 - 0.5 * theta + 0.5j * r).real
             + digamma(0.25 + 0.5 * theta + 0.5j * r).real) / 1.0 - math.log(math.pi)
        # two Gamma factors: sum of digammas, minus log pi once? Completed xi has
        # pi^{-s/2} each, so F_theta has pi^{-(s-theta)/2} pi^{-(s+theta)/2} = pi^{-s},
        # so -log pi  * 2 relative to single? xi has 1/2 s(s-1) pi^{-s/2} Gamma(s/2).
        # For product: pi^{-(s-theta)/2 - (s+theta)/2} = pi^{-s}, vs zeta's pi^{-s/2}.
        # Weight in explicit formula: d/ds log(pi^{-s/2}) = -1/2 log pi, etc.
        # Practical check: at theta=0 should reduce to 2*(digamma(1/4+..)-log pi) if we
        # naively add two identical factors — but F_0 = xi^2, degree 2.
        # For fair F_theta vs zeta comparison of renorm alone we already handle poles/primes.
        # For full Weil of F_theta use sum of two digammas - log(pi) * 2:
        w = (digamma(0.25 - 0.5 * theta + 0.5j * r).real
             + digamma(0.25 + 0.5 * theta + 0.5j * r).real
             - 2.0 * math.log(math.pi))
    # (1/2pi) int_{-inf}^{inf} = (1/pi) int_0^inf since even in r
    wr = np.full_like(r, dr); wr[0] *= 0.5; wr[-1] *= 0.5
    return float(np.sum(cg * w * wr) / math.pi)

def main():
    t0 = time.time()
    HS_pair = [0.4, 0.6, 0.8]
    centers = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    thetas = [0.25, 0.1]
    Nneed = int(math.floor(math.exp(12))) + 2  # c+w <= ~12
    Lam = mangoldt_upto(Nneed)

    seps = []
    zeta_neg = []
    rows = []
    for w in HS_pair:
        for c in centers:
            h = c + w + 0.05
            npts = max(2001, int(300 * h) | 1)
            u = np.linspace(-h, h, npts)
            g = paired_yoshida(u, c, w, 4) if c > 1e-12 else yoshida(u, w, 4)
            if float(np.max(np.abs(g))) < 1e-18:
                continue
            qz, pz, prz = Q_renorm(g, u, Lam, 0.0)
            rec = {"c": c, "w": w, "Q_zeta": qz, "pole_z": pz, "prime_z": prz, "Q_theta": {}}
            if qz < -1e-10:
                zeta_neg.append(rec)
            for th in thetas:
                qt, pt, prt = Q_renorm(g, u, Lam, th)
                rec["Q_theta"][str(th)] = qt
                if qz >= -1e-12 and qt < -1e-10:
                    seps.append({"c": c, "w": w, "theta": th, "Qz": qz, "Qt": qt})
            rows.append(rec)

    # larger centered yoshida up to h=12
    big = []
    for h in [6.0, 8.0, 10.0, 12.0]:
        npts = max(4001, int(200 * h) | 1)
        u = np.linspace(-h, h, npts)
        g = yoshida(u, h, 4)
        qz, _, _ = Q_renorm(g, u, Lam, 0.0)
        qt25, _, _ = Q_renorm(g, u, Lam, 0.25)
        big.append({"h": h, "Qz": qz, "Q25": qt25})
        if qz < -1e-10:
            zeta_neg.append({"c": 0, "w": h, "Q_zeta": qz, "note": "centered_big"})

    # Compare renorm vs full Weil on autocorr bumps at h where F_theta Weil fails (~0.46+)
    # Build F = f*f~ for f=yoshida on [-a,a], support F in [-2a,2a]
    compare = []
    for a in [0.25, 0.35, 0.5, 0.7, 1.0, 1.5]:
        h = 2 * a
        npts = max(2001, int(400 * h) | 1)
        u = np.linspace(-h, h, npts)
        xs = np.linspace(-a, a, 601); dx = xs[1] - xs[0]
        f = yoshida(xs, a, 4)
        F = np.array([np.sum(f * yoshida(xs - uu, a, 4)) * dx for uu in u])
        qz, _, _ = Q_renorm(F, u, Lam, 0.0)
        qt, _, _ = Q_renorm(F, u, Lam, 0.25)
        az = arch_term(F, u, 0.0)
        at = arch_term(F, u, 0.25)
        # full Weil ~ renorm - arch  (sign: explicit formula Sum = Pole - Prime - Arch
        # with Arch = (1/2pi)int hat(g) W(r); if W~ Re psi - log pi, and for g>0 hat(g)(0)>0,
        # arch at low freq is negative since digamma(1/4)-log pi < 0)
        compare.append({
            "a": a, "h": h,
            "Qz_renorm": qz, "Q25_renorm": qt,
            "arch_z": az, "arch_25": at,
            "Weil_z_est": qz - az, "Weil_25_est": qt - at,
        })

    out = {
        "status": "NUMERICAL_float",
        "separations_Qz_ge0_Qt_lt0": seps,
        "zeta_neg_count": len(zeta_neg),
        "zeta_neg_sample": zeta_neg[:20],
        "big_centered": big,
        "renorm_vs_weil_autocorr": compare,
        "n_paired_rows": len(rows),
        "elapsed": time.time() - t0,
        "paired_min_Qz": min(r["Q_zeta"] for r in rows),
        "paired_min_Q25": min(r["Q_theta"]["0.25"] for r in rows),
        "paired_min_Q10": min(r["Q_theta"]["0.1"] for r in rows),
    }
    print("=== v2 paired / large / weil compare ===")
    print(f"separations Qz>=0 Qt<0: {len(seps)}")
    for s in seps[:20]:
        print(" ", s)
    print(f"zeta_neg_count: {len(zeta_neg)}")
    print("big_centered:", big)
    print(f"paired min Qz={out['paired_min_Qz']:+.6e} Q25={out['paired_min_Q25']:+.6e} Q10={out['paired_min_Q10']:+.6e}")
    print("\n--- renorm vs Weil est on autocorr (F=f*f~) ---")
    for c in compare:
        print(f"  a={c['a']:.2f}: renorm Qz={c['Qz_renorm']:+.4e} Q25={c['Q25_renorm']:+.4e} | "
              f"arch z={c['arch_z']:+.4e} 25={c['arch_25']:+.4e} | "
              f"WeilEst z={c['Weil_z_est']:+.4e} 25={c['Weil_25_est']:+.4e}")
    path = os.path.join(os.path.dirname(__file__), "result_v2.json")
    # store compact
    slim_rows = [{"c": r["c"], "w": r["w"], "Qz": r["Q_zeta"],
                  "Q25": r["Q_theta"]["0.25"], "Q10": r["Q_theta"]["0.1"]} for r in rows]
    out["paired_rows"] = slim_rows
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print("wrote", path)

if __name__ == "__main__":
    main()
