#!/usr/bin/env python3
"""Round 16 PS: Satake/Frobenius pair (alpha, beta) for Ramanujan's Delta, and the
Round 15 Toeplitz Hasse Gram T[j,k] = sum_i omega_i^(k-j) on the normalised pair.

Inputs: integer coefficients tau(n) of Delta = q * prod (1-q^n)^24 (exact integer
arithmetic). No zeta zeros, no L-function zeros are used anywhere.

Checks:
  A. tau(n) sanity: known values, Hecke multiplicativity tau(mn)=tau(m)tau(n) (m,n coprime),
     prime-power recursion tau(p^2) = tau(p)^2 - p^11.
  B. 'trace = sum, norm = product': alpha+beta = tau(p), alpha*beta = p^11.
     Exact integer discriminant tau(p)^2 - 4 p^11 < 0 (Deligne bound) <=> the Round 15
     curve form m^2 - tau(p) m n + p^11 n^2 is positive definite.
  C. Toeplitz Gram G[j,k] = omega^(k-j) + omegabar^(k-j) = 2 cos((k-j) theta_p), sizes 1..8:
     min eigenvalue (expected >= -1e-9, rank <= 2).
  D. Fake: complementary-series-type parameters omega in {p^s, p^-s}, 0 < s < 1/2
     (norm = product kept: omega*omega' = 1; lambda_p = omega + 1/omega real, and below the
     unitary-dual bound p^(1/2) + p^(-1/2)); Gram 2 cosh((k-j) s log p): min eigenvalue.
     Also s = 1/2 (non-generic / Saito-Kurokawa-type Arthur parameter p^{+-1/2}).
  E. Growth form (Round 15 s2.4): |sum_i omega_i^k| <= 2 for all k vs 2 cosh(k s log p).
  G. Stepanov-falsifier sign test (team pairs): Dirichlet coefficients b(n) and log-derivative
     coefficients c(n) (-L'/L = sum c(n) n^-s) for ZETA, EPS (x^2+5y^2, b = r(n)/2) and DH
     (b periodic mod 5: 0,1,kappa,-kappa,-1). 'Count positivity' b(n) >= 0 vs 'von Mangoldt
     positivity' c(n) >= 0.
  F. Sym^n amplification: Sym^n of {p^s,p^-s} has extreme parameter p^{n s}; a unitary
     generic spherical rep of GL_{n+1}(Q_p) needs n s < 1/2, so 'all Sym^n unitary' forces s=0.
"""
import json, math
import numpy as np

NMAX = 200

def tau_coeffs(N):
    # coefficients of prod_{n>=1} (1-q^n)^24 up to q^(N-1), then shift by q
    c = [0] * N
    c[0] = 1
    for n in range(1, N):
        for _ in range(24):
            for m in range(N - 1, n - 1, -1):
                c[m] -= c[m - n]
    tau = [0] * (N + 1)
    for m in range(N):
        tau[m + 1] = c[m]
    return tau  # tau[n] for 1 <= n <= N

def primes(n):
    s = [True] * (n + 1); s[0] = s[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i): s[j] = False
    return [i for i in range(n + 1) if s[i]]

def toeplitz_from_power_sums(ps, K):
    # ps(d) = sum_i omega_i^d, d in Z; Gram G[j,k] = ps(k-j)
    G = np.array([[ps(k - j) for k in range(K)] for j in range(K)], dtype=float)
    return 0.5 * (G + G.T)

def main():
    tau = tau_coeffs(NMAX + 1)
    out = {"inputs": "integer tau(n), Delta = q prod(1-q^n)^24, n <= %d; no zeros used" % NMAX}

    # A. sanity
    known = {1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830, 6: -6048, 7: -16744, 11: 534612, 13: -577738}
    A = {"known_values_ok": all(tau[n] == v for n, v in known.items())}
    mult_fail = 0; mult_tests = 0
    for m in range(2, 15):
        for n in range(2, 15):
            if m * n <= NMAX and math.gcd(m, n) == 1:
                mult_tests += 1
                mult_fail += tau[m * n] != tau[m] * tau[n]
    A["hecke_multiplicative_tests"] = mult_tests
    A["hecke_multiplicative_failures"] = mult_fail
    pp = [p for p in primes(int(NMAX ** 0.5))]
    A["prime_square_recursion_ok"] = all(tau[p * p] == tau[p] ** 2 - p ** 11 for p in pp)
    out["A_sanity"] = A

    P = primes(NMAX)
    Ks = list(range(1, 9))
    rows = []
    worst = {K: 1e9 for K in Ks}
    max_abs_lambda = 0.0
    all_disc_neg = True
    for p in P:
        t = tau[p]
        disc = t * t - 4 * p ** 11            # exact integer
        all_disc_neg &= disc < 0
        lam = t / p ** 5.5                     # = omega + omegabar = 2 cos theta
        max_abs_lambda = max(max_abs_lambda, abs(lam))
        theta = math.acos(lam / 2)
        ps = lambda d, th=theta: 2 * math.cos(d * th)
        mins = {}
        for K in Ks:
            ev = np.linalg.eigvalsh(toeplitz_from_power_sums(ps, K))
            mins[K] = float(ev[0]); worst[K] = min(worst[K], mins[K])
        rows.append({"p": p, "tau_p": t, "disc_tau2_minus_4p11": str(disc),
                     "lambda_p": lam, "theta_p": theta,
                     "min_eig_by_size": {str(K): mins[K] for K in Ks}})
    out["B_trace_sum_norm_product"] = {
        "n_primes": len(P),
        "all_discriminants_negative": bool(all_disc_neg),
        "max_|lambda_p|=|tau(p)|/p^5.5": max_abs_lambda,
        "note": "disc<0 <=> m^2 - tau(p) m n + p^11 n^2 positive definite (Round 15 deg(m-nF) form)"}
    out["C_Delta_toeplitz"] = {"sizes": Ks,
                               "worst_min_eig_over_primes_by_size": {str(K): worst[K] for K in Ks},
                               "PSD_all(>=-1e-9)": bool(min(worst.values()) >= -1e-9),
                               "per_prime": rows}

    # D. complementary-series-type fakes
    fake = []
    sigmas = [0.01, 0.05, 7 / 64, 0.25, 0.49, 0.5]
    for p in [2, 3, 5, 7, 11, 101, 199]:
        for s in sigmas:
            x = s * math.log(p)
            ps = lambda d, x=x: 2 * math.cosh(d * x)
            mins = {str(K): float(np.linalg.eigvalsh(toeplitz_from_power_sums(ps, K))[0]) for K in Ks}
            lam = 2 * math.cosh(x)
            fake.append({"p": p, "sigma": s, "lambda_p": lam,
                         "unitary_dual_bound_p^.5+p^-.5": p ** 0.5 + p ** -0.5,
                         "lambda_below_unitary_bound": bool(lam < p ** 0.5 + p ** -0.5 or s == 0.5),
                         "generic_unitary_complementary_series": bool(s < 0.5),
                         "size2_min_eig_exact_2-2cosh(s log p)": 2 - 2 * math.cosh(x),
                         "min_eig_by_size": mins})
    out["D_fake_complementary"] = fake

    # E. growth form
    E = {"Delta_max_|sum omega^k|_k<=50": max(abs(2 * math.cos(k * r["theta_p"])) for r in rows for k in range(1, 51)),
         "fake_p=2_s=7/64_|sum omega^k|_k=10,50": [2 * math.cosh(k * 7 / 64 * math.log(2)) for k in (10, 50)]}
    out["E_growth"] = E

    # F. Sym^n amplification for the fake
    F = {}
    for s in [0.01, 7 / 64, 0.25, 0.49]:
        nmax = math.ceil(0.5 / s) - 1  # largest n with n*s < 1/2
        F[str(s)] = {"largest_n_with_Sym^n_parameter_in_unitary_generic_range(n*s<1/2)": nmax}
    out["F_symn_amplification"] = F

    # G. sign test
    def logderiv(b, N):
        c = [0.0] * (N + 1)
        for n in range(2, N + 1):
            acc = b[n] * math.log(n)
            for d in range(2, n):
                if n % d == 0:
                    acc -= c[d] * b[n // d]
            c[n] = acc  # b[1] = 1
        return c
    N = NMAX
    bz = [0.0] + [1.0] * N
    r = [0] * (N + 1)
    for x in range(-15, 16):
        for y in range(-7, 8):
            m = x * x + 5 * y * y
            if 0 < m <= N: r[m] += 1
    be = [0.0] + [r[n] / 2 for n in range(1, N + 1)]
    kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
    per = [0.0, 1.0, kap, -kap, -1.0]
    bd = [0.0] + [per[n % 5] for n in range(1, N + 1)]
    G = {"kappa_DH": kap}
    for name, b in (("ZETA", bz), ("EPS", be), ("DH", bd)):
        c = logderiv(b, N)
        negb = [n for n in range(1, N + 1) if b[n] < -1e-12]
        negc = [n for n in range(2, N + 1) if c[n] < -1e-9]
        G[name] = {"b_negative_count_n<=%d" % N: len(negb), "first_b_negative": negb[:5],
                   "c_negative_count_n<=%d" % N: len(negc), "first_c_negative": negc[:8],
                   "first_c_negative_values": [c[n] for n in negc[:4]]}
    out["G_stepanov_sign_falsifier"] = G
    print("G:", {k: (v if k == "kappa_DH" else {kk: vv for kk, vv in v.items() if 'values' not in kk}) for k, v in G.items()})

    with open("r16ps_satake.json", "w") as f:
        json.dump(out, f, indent=1)
    # summary
    print("A:", A)
    print("B: primes", len(P), "all disc<0:", all_disc_neg, "max|lambda_p|:", round(max_abs_lambda, 6))
    print("C: worst min eig by size:", {K: f"{worst[K]:.3e}" for K in Ks})
    for r in fake:
        if r["p"] in (2, 199):
            print("D: p=%d s=%.4f lam=%.4f  min eig K=2: %.4e  K=8: %.4e" % (
                r["p"], r["sigma"], r["lambda_p"], r["min_eig_by_size"]["2"], r["min_eig_by_size"]["8"]))
    print("E:", E)
    print("F:", F)

if __name__ == "__main__":
    main()
