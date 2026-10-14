"""Renormalized pole-minus-prime pairing: zeta vs F_theta.
NUMERICAL (float / optional mpmath). Not certified. RH OPEN.

Definitions (even compactly supported test function g):
  Q_zeta(g) = int_R g(u)(e^{u/2}+e^{-u/2}) du
            - 2 sum_{n>=2} Lambda(n) n^{-1/2} g(log n)

  Q_theta(g) = int_R g(u)(e^{u/2}+e^{-u/2})(e^{theta u}+e^{-theta u}) du
             - 2 sum_{n>=2} Lambda(n) n^{-1/2} (n^theta + n^{-theta}) g(log n)

Equivalent (even g, M_prime = 2 sum Lambda(n)n^{-1/2} delta_{log n}):
  Q_zeta = int_0^inf g(u)[2 e^{u/2} du - dM_prime] + 2 int_0^inf g(u) e^{-u/2} du

F_theta analogue replaces dM by dM_theta with a(n)=Lambda(n)(n^theta+n^{-theta}),
and pole density by (e^{u/2}+e^{-u/2})(e^{theta u}+e^{-theta u}).

This is NOT the full Weil form (archimedean / digamma term omitted).
We ask only: can Q_theta go negative on small-support even bumps while Q_zeta stays >=0?

Status labels in output: NUMERICAL (float).
"""
from __future__ import annotations
import json, math, os, sys, time
import numpy as np

THETAS = [0.25, 0.1, 0.05]
# support half-widths h; primes enter when log n <= h, i.e. n <= e^h
HS = [0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]

def mangoldt_upto(N: int):
    Lam = np.zeros(N + 1)
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if is_p[p]:
            is_p[p*p::p] = False
    for p in range(2, N + 1):
        if is_p[p]:
            q = p
            lp = math.log(p)
            while q <= N:
                Lam[q] = lp
                q *= p
    return Lam

def trap_int(fvals, u):
    return float(np.trapezoid(fvals, u))

# ---- even test functions on [-h,h] ----
def yoshida(u, h, k=4):
    """C^k-ish bump (1-(u/h)^2)^k on |u|<h; even."""
    out = np.zeros_like(u, dtype=float)
    m = np.abs(u) < h
    out[m] = (1.0 - (u[m] / h) ** 2) ** k
    return out

def gaussian_trunc(u, h, sigma):
    out = np.exp(-0.5 * (u / sigma) ** 2)
    out[np.abs(u) >= h] = 0.0
    return out

def bump_autocorr(u, h, k=4):
    """Autocorrelation F=f*f~ of Yoshida bump supported on [-h/2,h/2],
    so F supported on [-h,h] and of positive type (Fourier >=0)."""
    # F(u) = int f(x) f(x-u) dx with f = yoshida(., h/2, k)
    # compute by quadrature on a grid
    half = h / 2.0
    xs = np.linspace(-half, half, 801)
    dx = xs[1] - xs[0]
    f = yoshida(xs, half, k)
    Fu = np.zeros_like(u, dtype=float)
    for i, uu in enumerate(u):
        if abs(uu) >= h:
            continue
        Fu[i] = np.sum(f * yoshida(xs - uu, half, k)) * dx
    return Fu

def pole_density_zeta(u):
    return np.exp(0.5 * u) + np.exp(-0.5 * u)

def pole_density_theta(u, theta):
    return pole_density_zeta(u) * (np.exp(theta * u) + np.exp(-theta * u))

def Q_renorm(g_of_u, u_grid, Lam, theta=0.0):
    """g_of_u: array g(u_grid), even; u_grid covers [-h,h] densely."""
    # pole integral over R (grid is symmetric)
    if theta == 0.0:
        dens = pole_density_zeta(u_grid)
    else:
        dens = pole_density_theta(u_grid, theta)
    pole = trap_int(g_of_u * dens, u_grid)
    # prime sum: only n with |log n| <= max|u|
    hmax = float(np.max(np.abs(u_grid)))
    Nmax = int(math.floor(math.exp(hmax)))
    prime = 0.0
    # interpolate g at log n by linear interp on |u| (even => g(log n)=g(|log n|))
    abs_u = np.abs(u_grid)
    # sort unique nonnegative for interp
    order = np.argsort(abs_u)
    u_pos = abs_u[order]
    g_pos = g_of_u[order]
    # dedupe
    uniq_u, idx = np.unique(u_pos, return_index=True)
    uniq_g = g_pos[idx]
    for n in range(2, Nmax + 1):
        lam = Lam[n] if n < len(Lam) else 0.0
        if lam == 0.0:
            continue
        ln = math.log(n)
        if ln > hmax:
            break
        gval = float(np.interp(ln, uniq_u, uniq_g))
        if theta == 0.0:
            w = 1.0
        else:
            w = n**theta + n**(-theta)
        prime += 2.0 * lam * (n**-0.5) * w * gval
    return pole - prime, pole, prime

def scan():
    t0 = time.time()
    Nneed = int(math.floor(math.exp(max(HS)))) + 2
    Lam = mangoldt_upto(Nneed)
    rows = []
    families = []
    # family specs: (name, builder)
    for h in HS:
        # dense grid
        npts = max(2001, int(400 * h) | 1)
        u = np.linspace(-h, h, npts)
        for name, builder in [
            ("yoshida_k4", lambda uu, hh=h: yoshida(uu, hh, 4)),
            ("yoshida_k8", lambda uu, hh=h: yoshida(uu, hh, 8)),
            ("gauss_sig_h3", lambda uu, hh=h: gaussian_trunc(uu, hh, max(hh / 3, 0.05))),
            ("gauss_sig_h2", lambda uu, hh=h: gaussian_trunc(uu, hh, max(hh / 2, 0.05))),
            ("autocorr_yosh4", lambda uu, hh=h: bump_autocorr(uu, hh, 4)),
        ]:
            g = builder(u)
            # skip near-zero functions
            if float(np.max(np.abs(g))) < 1e-18:
                continue
            qz, pz, prz = Q_renorm(g, u, Lam, 0.0)
            entry = {
                "h": h,
                "family": name,
                "Q_zeta": qz,
                "pole_zeta": pz,
                "prime_zeta": prz,
                "Q_theta": {},
            }
            for th in THETAS:
                qt, pt, prt = Q_renorm(g, u, Lam, th)
                entry["Q_theta"][str(th)] = {"Q": qt, "pole": pt, "prime": prt}
            rows.append(entry)
            families.append((name, h, qz, {th: entry["Q_theta"][str(th)]["Q"] for th in THETAS}))

    # summarize separations
    sep = []  # Q_zeta >= 0 and some Q_theta < 0
    both_neg = []
    zeta_neg = []
    for name, h, qz, qts in families:
        if qz < -1e-10:
            zeta_neg.append((name, h, qz, qts))
        for th, qt in qts.items():
            if qz >= -1e-12 and qt < -1e-10:
                sep.append((name, h, th, qz, qt))
            if qz < -1e-10 and qt < -1e-10:
                both_neg.append((name, h, th, qz, qt))

    # also: Ritz-style — for each h, minimize Q over span of yoshida powers by sampling coeffs
    # cheap: evaluate on diagonal family only first; then a small Gram Ritz for yoshida basis
    ritz = []
    for h in HS:
        npts = max(2001, int(400 * h) | 1)
        u = np.linspace(-h, h, npts)
        # basis: even powers (1-(u/h)^2)^k for k=2,4,6,8 and autocorr
        basis = [
            yoshida(u, h, 2),
            yoshida(u, h, 4),
            yoshida(u, h, 6),
            yoshida(u, h, 8),
            gaussian_trunc(u, h, max(h / 3, 0.05)),
        ]
        B = np.stack(basis, axis=0)  # 5 x npts
        # build bilinear forms: Q(sum c_i b_i) = c^T A c, with A_ij = Q(mix of bi,bj)
        # For linear functional L(g)=pole-prime, A_ij = L( (bi+bj)/something )? No — Q is LINEAR in g!
        # Wait: Q_renorm is a LINEAR functional of g, not a quadratic form!
        # So "going negative" just means the functional takes a negative value on some positive bump.
        # Ritz over positive cones: min L(g) for g = sum a_i^2 basis_i or a_i>=0.
        # For linear L, the minimum over a simplex is at a vertex. So scanning the family IS enough
        # for positive combinations. For signed g (even but sign-changing), L can do anything.
        # Weil uses positive-type g; our autocorr family is the right one.
        # Record L on positive spans via vertices only (already in rows).
        pass

    out = {
        "status": "NUMERICAL_float",
        "definition": "Q = pole_int - prime_sum; arch omitted; RH OPEN",
        "thetas": THETAS,
        "hs": HS,
        "n_rows": len(rows),
        "separation_zeta_nonneg_theta_neg": [
            {"family": a, "h": b, "theta": c, "Q_zeta": d, "Q_theta": e}
            for a, b, c, d, e in sep
        ],
        "zeta_negative": [
            {"family": a, "h": b, "Q_zeta": c, "Q_theta": d}
            for a, b, c, d in zeta_neg
        ],
        "both_negative_count": len(both_neg),
        "elapsed_sec": time.time() - t0,
        "rows": rows,
    }
    return out

def main():
    out = scan()
    # print concise report
    print("=== Renormalized pole-minus-prime (arch OMITTED) ===")
    print("status: NUMERICAL (float); RH OPEN")
    print(f"rows={out['n_rows']} elapsed={out['elapsed_sec']:.2f}s")
    print()
    print("--- separations: Q_zeta>=0 and Q_theta<0 ---")
    seps = out["separation_zeta_nonneg_theta_neg"]
    if not seps:
        print("(none)")
    else:
        for s in seps[:40]:
            print(f"  {s['family']:16s} h={s['h']:.2f} theta={s['theta']}: "
                  f"Qe={s['Q_zeta']:+.6e}  Qt={s['Q_theta']:+.6e}")
        if len(seps) > 40:
            print(f"  ... +{len(seps)-40} more")
    print()
    print("--- Q_zeta < 0 (dead-end signal if common) ---")
    zn = out["zeta_negative"]
    if not zn:
        print("(none on scanned family)")
    else:
        for s in zn[:30]:
            print(f"  {s['family']:16s} h={s['h']:.2f}: Qe={s['Q_zeta']:+.6e}  Qt={s['Q_theta']}")
        if len(zn) > 30:
            print(f"  ... +{len(zn)-30} more")
    print()
    # table of Q_zeta and Q_{1/4} for autocorr and yoshida_k4
    print("--- sample table (yoshida_k4, autocorr_yosh4) Q_zeta / Q_0.25 ---")
    for r in out["rows"]:
        if r["family"] in ("yoshida_k4", "autocorr_yosh4") and r["h"] in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0):
            qt = r["Q_theta"]["0.25"]["Q"]
            print(f"  {r['family']:16s} h={r['h']:.1f}: "
                  f"Qz={r['Q_zeta']:+.6e}  Q25={qt:+.6e}  "
                  f"polez={r['pole_zeta']:+.4e}  primz={r['prime_zeta']:+.4e}")
    path = os.path.join(os.path.dirname(__file__), "result.json")
    # shrink rows for json: drop huge detail? keep all, it's small
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {path}")

if __name__ == "__main__":
    main()
