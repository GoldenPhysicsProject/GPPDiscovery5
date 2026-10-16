#!/usr/bin/env python3
"""H40a.F: Independent NUMERICAL sieve of prime-power arrivals for min S,E,U+.

Uses W' closed form (atanh+atan) + exact a=W'', float64 for speed with
mpmath omega0. NOT interval-certified. Compare to Codex mins.
"""
from __future__ import annotations
import argparse, json, math, os, time
import numpy as np
from mpmath import mp

OUT_DIR = os.path.dirname(__file__)

def run(T: float):
    mp.dps = 40
    omega0 = float(-mp.euler - mp.pi/2 - 3*mp.log(2) - mp.log(mp.pi))
    N = int(math.floor(math.exp(T))) + 1
    t0 = time.time()
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, math.isqrt(N) + 1):
        if sieve[p]:
            sieve[p*p::p] = False
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

    # W via sum1 closed + truncated exp; for float diagnostic use polygamma constant
    c = float(mp.polygamma(1, mp.mpf("0.25")) / 4)  # sum 1/a_j^2
    W = 8 * (np.cosh(t / 2) - 1) + (omega0 / 2) * t + c
    # W' closed form
    z = np.exp(-t / 2)
    # clip z to (0,1-eps) for atanh
    z = np.clip(z, 1e-300, 1 - 1e-15)
    Wprime = 4 * np.sinh(t / 2) + omega0 / 2 + np.arctanh(z) + np.arctan(z)
    for j in range(80):
        a = 2 * j + 0.5
        W -= np.exp(-a * t) / (a * a)
    # a(t)=W''(t)
    curvature = 2 * np.cosh(t / 2) - np.exp(-t / 2) / (1 - np.exp(-2 * t))
    # At event times: S continuous uses A = t*cumw - cumwu, which INCLUDES
    # the current kick's contribution to B but (t-u)_+=0 for current, so
    # A = sum_{u_n < t} ... + 0 for current. cumw includes current n.
    # Codex: S = W - (t*cumw - cumwu). At exact arrival, (t-u_n)=0 for current,
    # so t*w_n - w_n*t = 0 contribution from current — OK included.
    S = W - (t * cumw - cumwu)
    # slope AFTER kick: W' - cumw (cumw includes current w)
    slope = Wprime - cumw
    E = 2 * curvature * S - slope * slope
    X = np.sqrt(np.maximum(2 * curvature * S, 0.0))
    Uplus = X + slope
    Uminus = X - slope

    iminS = int(np.argmin(S))
    iminE = int(np.argmin(E))
    iminU = int(np.argmin(Uplus))
    iminUm = int(np.argmin(Uminus))
    out = {
        "T": T,
        "n_events": int(len(ns)),
        "last_n": int(ns[-1]),
        "min_S": {"value": float(S[iminS]), "n": int(ns[iminS]), "t": float(t[iminS])},
        "min_E": {"value": float(E[iminE]), "n": int(ns[iminE]), "t": float(t[iminE])},
        "min_Uplus": {"value": float(Uplus[iminU]), "n": int(ns[iminU]), "t": float(t[iminU])},
        "min_Uminus": {"value": float(Uminus[iminUm]), "n": int(ns[iminUm]), "t": float(t[iminUm])},
        "neg_E_count": int(np.count_nonzero(E < 0)),
        "neg_S_count": int(np.count_nonzero(S < 0)),
        "neg_Uplus_count": int(np.count_nonzero(Uplus < 0)),
        "Codex_min_S": {"value": 0.0277700093804, "n": 3109},
        "Codex_min_E": {"value": 0.0794953498652, "n": 2},
        "Codex_min_Uplus_note": "at n=2 ≈0.12188865 (through e^16)",
        "elapsed_sec": time.time() - t0,
        "label": "NUMERICAL",
    }
    # Compare
    out["match_min_S_n"] = int(ns[iminS]) == 3109
    out["match_min_E_n"] = int(ns[iminE]) == 2
    out["S_rel_diff_vs_Codex"] = abs(float(S[iminS]) - 0.0277700093804) / 0.0277700093804
    out["E_rel_diff_vs_Codex"] = abs(float(E[iminE]) - 0.0794953498652) / 0.0794953498652
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=12.0)
    args = ap.parse_args()
    out = run(args.T)
    path = os.path.join(OUT_DIR, f"h40a_sieve_T{int(args.T)}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print("Wrote", path)

if __name__ == "__main__":
    main()
