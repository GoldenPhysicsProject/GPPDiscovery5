#!/usr/bin/env python3
"""H41a: Independent CERTIFIED finite-prefix U+ / S / E sieve.

Rigorous outward enclosures via mpmath high-dps + explicit geometric Arch
tails. Independent of Codex Discovery2 diagnostic and of S40 float sieve.
Labels: CERTIFIED only when every event's lower bounds clear zero with
documented error budget; NUMERICAL stretch optional via --numerical-stretch.

Barrier convention (Codex/S40):
  At arrival t=log n with cumw INCLUDING current w:
    S = W - (t*cumw - cumwu)   # continuous; current contrib to A is 0
    slope_after = W' - cumw
    a = W''
    X = sqrt(2 a S)            # requires S>=0, a>0
    Uplus_after = X + slope_after
    margin m_n = Uplus_after   # = Uplus_before - w
  Barrier Uplus_before >= w  <=>  Uplus_after >= 0.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from mpmath import (
    mp, mpf, nstr, fabs, sqrt, pi, log, cosh, sinh, exp, euler,
    atanh, atan, catalan, digamma,
)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def set_dps(dps: int) -> None:
    mp.dps = dps


def omega0() -> mpf:
    """ω₀ = −γ − π/2 − 3 log 2 − log π = ψ(1/4) − log π."""
    o1 = -euler - pi / 2 - 3 * log(2) - log(pi)
    o2 = digamma(mpf("0.25")) - log(pi)
    if fabs(o1 - o2) > mpf("1e-40"):
        raise RuntimeError("omega0 identity failed at working precision")
    return o1


def sum_inv_a2_exact() -> mpf:
    """Σ_j 1/a_j² = π²/4 + 2 Catalan, a_j = 2j + 1/2."""
    return pi**2 / 4 + 2 * catalan


def arch_gamma_S(t: mpf, J: int, r_cap: mpf | None = None):
    """Return (gamma_S_lo, gamma_S_hi, tail) for Σ (1−e^{-a t})/a².

    gamma_S = sum_inv_exact − Σ e^{-a t}/a².
    Truncate Σ e at j<J; geometric tail bound with ratio r=e^{-2t}
    (or r_cap if smaller for safety when e^{-2t}≈1 at tiny t — not used
    for t>=log2).
    """
    sum_inv = sum_inv_a2_exact()
    sum_exp = mpf(0)
    for j in range(J):
        a = 2 * j + mpf("0.5")
        sum_exp += exp(-a * t) / (a * a)
    aJ = 2 * J + mpf("0.5")
    eaJ = exp(-aJ * t)
    r = exp(-2 * t)
    if r_cap is not None and r > r_cap:
        r = r_cap
    # If r >= 1 (t<=0) fail; for t>=log2, r<=1/4.
    if r >= 1:
        raise ValueError("geometric ratio r>=1")
    # Tail of Σ_{j>=J} e^{-a t}/a² ≤ (eaJ/aJ²)/(1-r)
    # (ratio of consecutive terms ≤ e^{-2t} = r, and 1/a² decreases).
    first = eaJ / (aJ * aJ)
    tail = first / (1 - r)
    # sum_exp_true ∈ [sum_exp, sum_exp+tail]
    # gamma_S ∈ [sum_inv - sum_exp - tail, sum_inv - sum_exp]
    return sum_inv - sum_exp - tail, sum_inv - sum_exp, tail


def W_enclosure(t: mpf, omega: mpf, J: int):
    """Enclosure of W(t) = 8(cosh(t/2)-1)+(ω₀/2)t + gamma_S."""
    g_lo, g_hi, tail = arch_gamma_S(t, J)
    base = 8 * (cosh(t / 2) - 1) + (omega / 2) * t
    # Closed base uses mpmath elementary; pad rounding below separately.
    return base + g_lo, base + g_hi, tail


def Wprime_closed(t: mpf, omega: mpf) -> mpf:
    z = exp(-t / 2)
    return 4 * sinh(t / 2) + omega / 2 + atanh(z) + atan(z)


def Wpp_closed(t: mpf) -> mpf:
    return 2 * cosh(t / 2) - exp(-t / 2) / (1 - exp(-2 * t))


def choose_J(t: mpf, target_tail: mpf = mpf("1e-25")) -> int:
    """Pick truncation so geometric Arch tail for S is < target_tail."""
    # Crude: need (e^{-a_J t}/a_J²)/(1-e^{-2t}) < target
    # a_J ≈ 2J, so e^{-2J t} / (1-e^{-2t}) ≲ target * a_J²
    # Use J so that e^{-(2J)t} < target * (1-r) * (something)
    r = exp(-2 * t)
    # Increase J until first/(1-r) < target
    for J in range(20, 5000):
        aJ = 2 * J + mpf("0.5")
        first = exp(-aJ * t) / (aJ * aJ)
        if first / (1 - r) < target_tail:
            return J
    return 5000


def sieve_prime_powers(N: int):
    """Exact integer sieve: list of (n, p, k) with n=p^k <= N, Λ(n)=log p."""
    is_prime = bytearray(b"\x01") * (N + 1)
    is_prime[:2] = b"\x00\x00"
    for p in range(2, int(math.isqrt(N)) + 1):
        if is_prime[p]:
            step = p
            start = p * p
            is_prime[start : N + 1 : step] = b"\x00" * (((N - start) // step) + 1)
    events = []
    for p in range(2, N + 1):
        if not is_prime[p]:
            continue
        n = p
        k = 1
        while n <= N:
            events.append((n, p, k))
            # next power
            if n > N // p:
                break
            n *= p
            k += 1
    events.sort(key=lambda x: x[0])
    return events


def rounding_budget(dps: int, n_ops: int = 50) -> mpf:
    """Conservative absolute rounding pad for elementary closed forms."""
    # ~ n_ops ulps at working precision; use 10^{-(dps-8)} generously
    return mpf(10) ** (-(dps - 8)) * n_ops


def run_certified(T: float, dps: int = 40, target_tail: float = 1e-25):
    set_dps(dps)
    t_wall0 = time.time()
    omega = omega0()
    N = int(math.floor(math.exp(T)))
    events = sieve_prime_powers(N)
    # Restrict to n <= exp(T) exactly
    expT = exp(mpf(T))
    events = [(n, p, k) for (n, p, k) in events if mpf(n) <= expT]

    round_pad = rounding_budget(dps)
    target_tail_m = mpf(str(target_tail))

    # Running sums as mpf (logs computed at working dps)
    cumw = mpf(0)
    cumwu = mpf(0)

    # Track minima of lower bounds
    mins = {
        "S_lo": None,
        "E_lo": None,
        "Uplus_after_lo": None,  # = margin m_n
        "Uplus_before_minus_w_lo": None,  # should equal Uplus_after
    }
    argmins = {}
    first_fail = None
    n_fail = 0
    records_sample = []  # first few + argmins later

    # Precompute J schedule by t-buckets for speed
    # For t>=log2, J decreases as t grows
    J_cache = {}

    def J_for(t):
        # bucket by 0.1
        key = int(float(t) * 10)
        if key not in J_cache:
            J_cache[key] = choose_J(t, target_tail_m)
        return J_cache[key]

    n_events = len(events)
    for idx, (n, p, k) in enumerate(events):
        t = log(mpf(n))
        # w = Λ(n)/√n = log(p) / n^{1/2}
        w = log(mpf(p)) / sqrt(mpf(n))
        # Include current kick in cumsums (A contrib of current is 0 at exact t)
        cumw = cumw + w
        cumwu = cumwu + w * t

        J = J_for(t)
        W_lo, W_hi, arch_tail = W_enclosure(t, omega, J)
        # Pad closed-form base rounding
        W_lo = W_lo - round_pad
        W_hi = W_hi + round_pad

        A = t * cumw - cumwu  # exact in mpf arithmetic of our sums
        # S = W - A; A is sum of mpf logs — pad accumulation
        # Accumulation pad: ~ idx * round_pad for sum of w's (very conservative)
        acc_pad = round_pad * (idx + 1) * mpf("0.01")  # scaled; dps=40 ≫ need
        # Better: relative ~ 2^{-dps*log2(10)} * magnitude * n_ops
        acc_pad = mpf(10) ** (-(dps - 10)) * (idx + 1)

        S_lo = W_lo - A - acc_pad
        S_hi = W_hi - A + acc_pad

        Wp = Wprime_closed(t, omega)
        a = Wpp_closed(t)
        # Pad Wp, a
        Wp_lo, Wp_hi = Wp - round_pad, Wp + round_pad
        a_lo, a_hi = a - round_pad, a + round_pad
        if a_lo <= 0:
            raise RuntimeError(f"a_lo<=0 at n={n}")

        # slope_after = W' - cumw
        slope_lo = Wp_lo - cumw - acc_pad
        slope_hi = Wp_hi - cumw + acc_pad

        # For lower bound of X=sqrt(2aS): need lower of 2aS
        # If S_lo < 0, X_lo must be treated carefully — barrier fails or GAP
        two_a_S_lo = 2 * a_lo * S_lo  # if S_lo>0 and a_lo>0 this is valid lower
        # More careful: 2aS minimized over a∈[a_lo,a_hi], S∈[S_lo,S_hi]
        # If S_lo>=0: min is 2*a_lo*S_lo
        # If S_hi<0: definitely negative S — fail
        # If S straddles 0: S_lo<0 — fail for positivity

        ok = True
        fail_reason = None
        if S_lo < 0:
            ok = False
            fail_reason = "S_lo<0"
            X_lo = mpf(0)  # cannot claim
            E_lo = mpf("-inf")  # sentinel
            Up_after_lo = mpf("-inf")
        else:
            two_a_S_lo = 2 * a_lo * S_lo
            X_lo = sqrt(two_a_S_lo)
            X_hi = sqrt(2 * a_hi * S_hi)
            # E = 2aS - slope^2; lower bound:
            # 2aS_lo - max(slope^2 over slope interval)
            # max slope^2 = max(slope_lo^2, slope_hi^2) if same sign, else could be 0
            s2_candidates = [slope_lo * slope_lo, slope_hi * slope_hi]
            if slope_lo <= 0 <= slope_hi:
                s2_max = max(s2_candidates)  # 0 is min, max is endpoint
            else:
                s2_max = max(s2_candidates)
            E_lo = two_a_S_lo - s2_max
            # Uplus_after = X + slope; lower = X_lo + slope_lo
            Up_after_lo = X_lo + slope_lo
            # Uplus_before - w = (X + slope_after + w) - w = X + slope_after
            # so same as Up_after
            if E_lo < 0:
                ok = False
                fail_reason = "E_lo<0"
            if Up_after_lo < 0:
                ok = False
                fail_reason = (fail_reason or "") + "+Uplus_after_lo<0"

        if not ok:
            n_fail += 1
            if first_fail is None:
                first_fail = {
                    "n": n, "p": p, "k": k, "t": nstr(t, 20),
                    "reason": fail_reason,
                    "S_lo": nstr(S_lo, 20),
                    "E_lo": nstr(E_lo, 20) if E_lo != mpf("-inf") else "-inf",
                    "Uplus_after_lo": nstr(Up_after_lo, 20) if Up_after_lo != mpf("-inf") else "-inf",
                }

        def upd(key, val, n_, t_, extra=None):
            if val == mpf("-inf"):
                return
            cur = mins[key]
            if cur is None or val < cur:
                mins[key] = val
                argmins[key] = {"n": n_, "t": nstr(t_, 20), "value": nstr(val, 30)}
                if extra:
                    argmins[key].update(extra)

        upd("S_lo", S_lo, n, t)
        if E_lo != mpf("-inf"):
            upd("E_lo", E_lo, n, t)
        if Up_after_lo != mpf("-inf"):
            upd("Uplus_after_lo", Up_after_lo, n, t)
            upd("Uplus_before_minus_w_lo", Up_after_lo, n, t)

        if idx < 3 or n in (2, 3, 4, 5, 3109):
            records_sample.append({
                "n": n, "t": nstr(t, 15),
                "S_lo": nstr(S_lo, 20), "S_hi": nstr(S_hi, 20),
                "E_lo": nstr(E_lo, 20) if E_lo != mpf("-inf") else "-inf",
                "Uplus_after_lo": nstr(Up_after_lo, 20) if Up_after_lo != mpf("-inf") else "-inf",
                "w": nstr(w, 15), "J": J, "arch_tail": nstr(arch_tail, 8),
                "a": nstr(a, 15),
            })

        if (idx + 1) % 20000 == 0:
            elapsed = time.time() - t_wall0
            print(f"  ... {idx+1}/{n_events} events, elapsed {elapsed:.1f}s, "
                  f"min U+_lo={nstr(mins['Uplus_after_lo'], 10) if mins['Uplus_after_lo'] is not None else '?'}",
                  flush=True)

    elapsed = time.time() - t_wall0
    all_ok = n_fail == 0 and first_fail is None
    # CERTIFIED iff all lower bounds nonnegative and we used rigorous tails
    label = "CERTIFIED" if all_ok else "GAP"

    def fmt_min(key):
        if mins[key] is None:
            return None
        return {
            "lower_bound": nstr(mins[key], 30),
            "argmin": argmins.get(key),
            "positive": float(mins[key]) > 0,
        }

    out = {
        "label": label,
        "T": T,
        "N_max": N,
        "n_events": n_events,
        "dps": dps,
        "target_arch_tail": target_tail,
        "round_pad": nstr(round_pad, 10),
        "n_fail": n_fail,
        "first_fail": first_fail,
        "min_S_lo": fmt_min("S_lo"),
        "min_E_lo": fmt_min("E_lo"),
        "min_Uplus_after_lo": fmt_min("Uplus_after_lo"),
        "min_margin_m_n": fmt_min("Uplus_after_lo"),  # m_n = U+_after
        "barrier_holds_prefix": all_ok,
        "sample_events": records_sample,
        "Codex_float_mins_T14": {
            "min_S": {"value": 0.0277700093804, "n": 3109},
            "min_E": {"value": 0.0794953498652, "n": 2},
            "min_Uplus": {"value": 0.12188865, "n": 2, "note": "through e^16 float"},
        },
        "elapsed_sec": elapsed,
        "method": (
            "exact integer Λ sieve; mpmath mpf closed W',W''; "
            "Arch Σ via π²/4+2G − truncated exp + geometric tail; "
            "outward pads for rounding/accumulation; independent of Codex/S40 float"
        ),
        "RH": "OPEN",
        "claims": {
            "prefix_Uplus_barrier": label,
            "global_kick_barrier": "OPEN",
        },
    }
    return out


def run_numerical_stretch(T: float):
    """Float64 stretch beyond certified T — NUMERICAL only. Independent impl."""
    import numpy as np
    t0 = time.time()
    set_dps(40)
    omega = float(omega0())
    N = int(math.floor(math.exp(T))) + 1
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, math.isqrt(N) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    lam = np.zeros(N + 1, dtype=np.float64)
    for p in np.flatnonzero(sieve):
        n = int(p)
        lp = math.log(p)
        while n <= N:
            lam[n] = lp
            n *= int(p)
    ns = np.flatnonzero(lam)
    ns = ns[(ns >= 2) & (ns <= math.exp(T))]
    t = np.log(ns.astype(np.float64))
    w = lam[ns] / np.sqrt(ns.astype(np.float64))
    cumw = np.cumsum(w)
    cumwu = np.cumsum(w * t)
    c = float(sum_inv_a2_exact())
    W = 8 * (np.cosh(t / 2) - 1) + (omega / 2) * t + c
    z = np.clip(np.exp(-t / 2), 1e-300, 1 - 1e-15)
    Wprime = 4 * np.sinh(t / 2) + omega / 2 + np.arctanh(z) + np.arctan(z)
    for j in range(120):
        a = 2 * j + 0.5
        W -= np.exp(-a * t) / (a * a)
    curvature = 2 * np.cosh(t / 2) - np.exp(-t / 2) / (1 - np.exp(-2 * t))
    S = W - (t * cumw - cumwu)
    slope = Wprime - cumw
    E = 2 * curvature * S - slope * slope
    X = np.sqrt(np.maximum(2 * curvature * S, 0.0))
    Uplus = X + slope  # after
    iminS = int(np.argmin(S))
    iminE = int(np.argmin(E))
    iminU = int(np.argmin(Uplus))
    return {
        "label": "NUMERICAL",
        "T": T,
        "n_events": int(len(ns)),
        "min_S": {"value": float(S[iminS]), "n": int(ns[iminS]), "t": float(t[iminS])},
        "min_E": {"value": float(E[iminE]), "n": int(ns[iminE]), "t": float(t[iminE])},
        "min_Uplus_after": {"value": float(Uplus[iminU]), "n": int(ns[iminU]), "t": float(t[iminU])},
        "neg_S": int(np.count_nonzero(S < 0)),
        "neg_E": int(np.count_nonzero(E < 0)),
        "neg_Uplus": int(np.count_nonzero(Uplus < 0)),
        "elapsed_sec": time.time() - t0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=14.0, help="certified max t")
    ap.add_argument("--dps", type=int, default=40)
    ap.add_argument("--target-tail", type=float, default=1e-25)
    ap.add_argument("--numerical-stretch", type=float, default=0.0,
                    help="if >T, also run float NUMERICAL up to this T")
    ap.add_argument("--skip-cert", action="store_true")
    args = ap.parse_args()

    results = {"session": "S41", "H41a": {}}

    if not args.skip_cert:
        print(f"H41a CERTIFIED sieve T={args.T}, dps={args.dps}", flush=True)
        cert = run_certified(args.T, dps=args.dps, target_tail=args.target_tail)
        results["H41a"]["certified"] = cert
        path = os.path.join(OUT_DIR, f"h41a_cert_T{int(args.T)}.json")
        with open(path, "w") as f:
            json.dump(cert, f, indent=2)
        print(json.dumps({
            "label": cert["label"],
            "n_events": cert["n_events"],
            "n_fail": cert["n_fail"],
            "min_S_lo": cert["min_S_lo"],
            "min_E_lo": cert["min_E_lo"],
            "min_Uplus_after_lo": cert["min_Uplus_after_lo"],
            "elapsed_sec": cert["elapsed_sec"],
        }, indent=2), flush=True)
        print("Wrote", path, flush=True)

    if args.numerical_stretch and args.numerical_stretch > 0:
        print(f"H41a NUMERICAL stretch T={args.numerical_stretch}", flush=True)
        num = run_numerical_stretch(args.numerical_stretch)
        results["H41a"]["numerical_stretch"] = num
        path = os.path.join(OUT_DIR, f"h41a_numerical_T{int(args.numerical_stretch)}.json")
        with open(path, "w") as f:
            json.dump(num, f, indent=2)
        print(json.dumps(num, indent=2), flush=True)
        print("Wrote", path, flush=True)

    summary_path = os.path.join(OUT_DIR, "h41a_summary.json")
    # compact summary
    summary = {
        "certified_T": args.T if not args.skip_cert else None,
        "certified_label": results.get("H41a", {}).get("certified", {}).get("label"),
        "numerical_T": args.numerical_stretch or None,
        "RH": "OPEN",
    }
    if "certified" in results.get("H41a", {}):
        c = results["H41a"]["certified"]
        summary["min_Uplus_after_lo"] = c.get("min_Uplus_after_lo")
        summary["min_S_lo"] = c.get("min_S_lo")
        summary["min_E_lo"] = c.get("min_E_lo")
        summary["n_events"] = c.get("n_events")
        summary["barrier_holds_prefix"] = c.get("barrier_holds_prefix")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print("Wrote", summary_path, flush=True)


if __name__ == "__main__":
    main()
