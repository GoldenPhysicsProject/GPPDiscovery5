#!/usr/bin/env python3
"""H40a.D/E/F-opt: Algebraic identity numerical spot-check + optional Fenchel slack.

D identities are elementary (PROVED in note). Here: float spot-check ΔE and U± jumps.
E Landau logic is audited in the markdown note (no code proof).
"""
from __future__ import annotations
import json, math, os
import numpy as np
from mpmath import mp

OUT = os.path.join(os.path.dirname(__file__), "h40a_algebra_fenchel.json")
mp.dps = 30
omega0 = float(-mp.euler - mp.pi/2 - 3*mp.log(2) - mp.log(mp.pi))

def Wprime(t):
    z = math.exp(-t/2)
    return 4*math.sinh(t/2) + omega0/2 + math.atanh(min(z, 1-1e-15)) + math.atan(z)

def Wsecond(t):
    return 2*math.cosh(t/2) - math.exp(-t/2)/(1 - math.exp(-2*t))

def W_val(t, J=60):
    c = float(mp.polygamma(1, mp.mpf("0.25"))/4)
    W = 8*(math.cosh(t/2)-1) + (omega0/2)*t + c
    for j in range(J):
        a = 2*j + 0.5
        W -= math.exp(-a*t)/(a*a)
    return W

# Spot-check kick identities at a few prime powers using high-res float
T = 8.0
N = int(math.floor(math.exp(T))) + 1
sieve = np.ones(N+1, dtype=bool); sieve[:2]=False
for p in range(2, math.isqrt(N)+1):
    if sieve[p]: sieve[p*p::p]=False
lam = np.zeros(N+1)
for p in np.flatnonzero(sieve):
    n=int(p); lp=math.log(p)
    while n<=N:
        lam[n]=lp; n*=int(p)
ns = [int(n) for n in np.flatnonzero(lam) if 2<=n<=math.exp(T)]

kick_checks = []
# Walk events, verify ΔE = 2w p_minus - w^2 and U+ jump = -w
cumw = 0.0; cumwu = 0.0
prev_p_after = None
for n in ns[:200]:  # first 200 events
    t = math.log(n)
    w = lam[n]/math.sqrt(n)
    # Before kick: B = cumw (not yet including n)
    S = W_val(t) - (t*cumw - cumwu)  # A before adding current (current contrib 0)
    p_minus = Wprime(t) - cumw
    a = Wsecond(t)
    E_minus = 2*a*S - p_minus**2
    X = math.sqrt(max(2*a*S, 0))
    Up_minus = X + p_minus
    Um_minus = X - p_minus
    # After kick
    p_plus = p_minus - w
    E_plus = 2*a*S - p_plus**2
    Up_plus = X + p_plus
    Um_plus = X - p_plus
    dE = E_plus - E_minus
    dE_formula = 2*w*p_minus - w*w
    kick_checks.append({
        "n": n,
        "dE": dE,
        "dE_formula": dE_formula,
        "dE_err": abs(dE - dE_formula),
        "dUp": Up_plus - Up_minus,
        "dUp_plus_w": (Up_plus - Up_minus) + w,
        "dUm_minus_w": (Um_plus - Um_minus) - w,
    })
    cumw += w
    cumwu += w*t

max_dE_err = max(c["dE_err"] for c in kick_checks)
max_dUp_err = max(abs(c["dUp_plus_w"]) for c in kick_checks)
max_dUm_err = max(abs(c["dUm_minus_w"]) for c in kick_checks)

# Optional Fenchel slack at a few t: C_t - W_*(B_t)
# W_*(b) = v*b - W(v) where W'(v)=b, v>log2
def solve_v(b, lo=math.log(2)+1e-12, hi=None):
    if hi is None:
        hi = max(lo+1, 2*math.log(max(b,1)+1)+5)
    for _ in range(80):
        mid = 0.5*(lo+hi)
        if Wprime(mid) < b:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)

fenchel = []
cumw=0.0; cumwu=0.0
sample_ns = [n for n in ns if n in (2,3,5,16,100,1000,3109) or n==ns[min(500,len(ns)-1)]]
# recompute properly at selected events
cumw=0.0; cumwu=0.0
for n in ns:
    t=math.log(n); w=lam[n]/math.sqrt(n)
    cumw += w; cumwu += w*t
    if n not in sample_ns and n not in (2,3,5,16):
        continue
    B=cumw; C=cumwu
    if B <= Wprime(math.log(2))+1e-14:
        continue
    v = solve_v(B)
    Wstar = v*B - W_val(v)
    slack = C - Wstar
    S = W_val(t) - (t*B - C)
    fenchel.append({"n": n, "t": t, "B": B, "C": C, "v": v, "slack": slack, "S": S})

results = {
    "algebra_spotcheck": {
        "n_events_checked": len(kick_checks),
        "max_dE_err": max_dE_err,
        "max_dUp_err": max_dUp_err,
        "max_dUm_err": max_dUm_err,
        "label": "PROVED (elementary); NUMERICAL spot-check OK" if max_dE_err < 1e-10 else "GAP",
    },
    "fenchel_samples": fenchel,
    "fenchel_label": "NUMERICAL",
    "landau_audit": "see note — CONFIRMED skeleton with listed GAPs",
}
with open(OUT,"w") as f:
    json.dump(results, f, indent=2)
print("max_dE_err", max_dE_err, "max_dUp", max_dUp_err)
print("fenchel samples", len(fenchel))
for row in fenchel[:8]:
    print(row["n"], "slack", row["slack"], "S", row["S"])
print("Wrote", OUT)
