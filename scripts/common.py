"""Shared helpers for GPPDiscovery5 (Grok's lab).

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

Everything here is FLOATING-POINT (mpmath) unless a caller states otherwise.
Nothing here is an interval certificate.
"""
from __future__ import annotations
import json, platform, subprocess, time
from mpmath import mp, mpf, mpc, sqrt, pi, gamma, zeta, dirichlet

def xi(s):
    """Riemann xi(s) = 1/2 s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)."""
    s = mpc(s)
    return s * (s - 1) / 2 * pi ** (-s / 2) * gamma(s / 2) * zeta(s)

def F_theta(s, theta):
    """Falsifier F_theta(s) = xi(s-theta) xi(s+theta) (0 < theta < 1/2)."""
    return xi(mpc(s) - theta) * xi(mpc(s) + theta)

# Davenport-Heilbronn (conductor 5, odd): f = ((1-i k)/2) L(s,chi) + ((1+i k)/2) L(s,chibar),
# chi mod 5 with chi(2) = i, k = (sqrt(10-2 sqrt5) - 2)/(sqrt5 - 1).
_CHI = [0, 1, 1j, -1j, -1]
_CHIB = [0, 1, -1j, 1j, -1]

def dh_kappa():
    return (sqrt(10 - 2 * sqrt(5)) - 2) / (sqrt(5) - 1)

def dh(s):
    k = dh_kappa()
    return (1 - 1j * k) / 2 * dirichlet(s, _CHI) + (1 + 1j * k) / 2 * dirichlet(s, _CHIB)

def dh_completed(s):
    """Lambda_DH(s) = (5/pi)^{s/2} Gamma((s+1)/2) f(s); satisfies Lambda(s) = Lambda(1-s)."""
    s = mpc(s)
    return (5 / pi) ** (s / 2) * gamma((s + 1) / 2) * dh(s)

def _git_rev():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return None

def write_result(path, data, *, kind="float", method="", dps=None):
    """Write result.json-style output stamped with precision and float-vs-certified.

    kind: "float" (mpmath/numpy; no rigorous error bound) or "certified" (interval/ball
    arithmetic with stated enclosures). Default float.
    """
    assert kind in ("float", "certified")
    out = {
        "kind": kind,
        "method": method,
        "dps": dps if dps is not None else mp.dps,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git": _git_rev(),
        "python": platform.python_version(),
        "credit": "Daniel Toupin / Golden Physics Project, with AI assistance from Grok.",
    }
    out.update(data)
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    return out
