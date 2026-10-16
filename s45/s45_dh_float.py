"""S45 H45a: NUMERICAL (float) DH Weil Ritz on Legendre·(1-u²)^4 basis — same family as S14/S44.
Reuses experiments/weil_small_support/weil_ritz.py logic (DH arch digamma(3/4)+log(5/π),
c(n) from -f'/f with κ=√(φ+2)-φ from q=5). NOT Arb-CERTIFIED (arch not in S14 exact form).
φ from choice q=5 (INPUT-CHOICE; genuine for conductor-5 DH).
"""
import sys, os, json, numpy as np
from numpy.polynomial import legendre as L
from scipy.special import digamma
from scipy.linalg import eigh

M = int(sys.argv[1]) if len(sys.argv) > 1 else 24
hs = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [0.55, 0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5]

xq0, wq0 = L.leggauss(1600)
# κ = √(φ+2) − φ = (√(10−2√5)−2)/(√5−1); φ from q=5
kap = (np.sqrt(10 - 2 * np.sqrt(5)) - 2) / (np.sqrt(5) - 1)
phi = (1 + np.sqrt(5)) / 2

def vm(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return np.log(p) if m == 1 else 0.0
    return 0.0

def dh_a(n):
    chi = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
    return ((1 - 1j * kap) / 2 * chi[n % 5] + (1 + 1j * kap) / 2 * np.conj(chi[n % 5])).real

Nmax = int(np.exp(2 * max(hs))) + 5
c_dh = {1: 0.0}
for n in range(2, Nmax):
    c_dh[n] = dh_a(n) * np.log(n) - sum(c_dh[d] * dh_a(n // d) for d in range(2, n) if n % d == 0)

def basis(x, h, M):
    u = x / h
    B = np.array([L.legval(u, [0] * j + [1]) for j in range(M)]) * (1 - u ** 2) ** 4
    B[:, np.abs(u) >= 1] = 0
    return B

def form_dh(h, M):
    xq, wq = xq0 * h, wq0 * h
    B = basis(xq, h, M)
    gram = (B * wq) @ B.T
    R = 800.0 / h
    r = np.linspace(0, R, 40001)
    dr = r[1] - r[0]
    wr = np.full_like(r, dr); wr[0] = wr[-1] = dr / 2
    E = lambda z: (B * wq) @ np.exp(np.outer(xq, z))
    Q = np.zeros((M, M))
    w = digamma(0.75 + 0.5j * r).real + np.log(5 / np.pi)
    Er = E(1j * r)
    Q += 2 * np.real((Er * (w * wr)) @ np.conj(Er).T) / (2 * np.pi)
    def Fmat(u):
        Bs = basis(xq - u, h, M)
        A = (B * wq) @ Bs.T
        return (A + A.T) / 2
    for n in range(2, int(np.exp(2 * h)) + 1):
        if c_dh.get(n, 0) != 0:
            Q -= 2 * c_dh[n] / np.sqrt(n) * Fmat(np.log(n))
    return Q, gram

def form_ftheta_full(h, M, theta):
    """Full F_θ float form (shifted E, zeta digamma, poles at ±1/2±θ) — NUMERICAL."""
    xq, wq = xq0 * h, wq0 * h
    B = basis(xq, h, M)
    gram = (B * wq) @ B.T
    R = 800.0 / h
    r = np.linspace(0, R, 40001)
    dr = r[1] - r[0]
    wr = np.full_like(r, dr); wr[0] = wr[-1] = dr / 2
    E = lambda z: (B * wq) @ np.exp(np.outer(xq, z))
    Q = np.zeros((M, M))
    w = digamma(0.25 + 0.5j * r).real - np.log(np.pi)
    def Fmat(u):
        Bs = basis(xq - u, h, M)
        A = (B * wq) @ Bs.T
        return (A + A.T) / 2
    for sg in ([0.0] if theta == 0 else [theta, -theta]):
        Ep, Em = E(1j * r + sg), E(-1j * r - sg)
        A = (Ep * (w * wr)) @ Em.T
        Q += 2 * np.real((A + A.T) / 2) / (2 * np.pi)
        for a in (0.5, -0.5):
            e1 = E(np.array([a + sg]))[:, 0]
            e2 = E(np.array([-a - sg]))[:, 0]
            P = np.outer(e1, e2).real
            Q += (P + P.T) / 2
        for n in range(2, int(np.exp(2 * h)) + 1):
            lv = vm(n)
            if lv:
                Q -= lv / np.sqrt(n) * (n ** sg + n ** (-sg)) * Fmat(np.log(n))
    return Q, gram

def lam_sectors(Q, G, M):
    out = {}
    for nm, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        ev = eigh(Q[np.ix_(I, I)], G[np.ix_(I, I)], eigvals_only=True)
        out[nm] = {"lam_min": float(ev[0]), "negative": bool(ev[0] < 0)}
    return out

rows = []
# DH scan
for h in hs:
    Q, G = form_dh(h, M)
    r = {"kind": "DH", "M": M, "h": h, "label": "NUMERICAL",
         "phi": float(phi), "kappa": float(kap),
         "phi_status": "from choice q=5 (INPUT-CHOICE; genuine for conductor-5 DH)"}
    r.update(lam_sectors(Q, G, M))
    rows.append(r)
    print(json.dumps(r), flush=True)

# Full F_θ float at onset grid (for comparison to prime-def CERTIFIED)
for theta, hgrid in ((0.25, [0.36, 0.40, 0.44, 0.46, 0.48, 0.50, 0.55]),
                     (0.05, [0.40, 0.48, 0.50, 0.52, 0.54, 0.56])):
    for h in hgrid:
        Q, G = form_ftheta_full(h, M, theta)
        r = {"kind": "F_theta_full", "M": M, "h": h, "theta": theta, "label": "NUMERICAL"}
        r.update(lam_sectors(Q, G, M))
        rows.append(r)
        print(json.dumps(r), flush=True)

with open("dh_ftheta_float.json", "w") as f:
    json.dump({"rows": rows, "phi_status": "DH κ uses φ from q=5 (INPUT-CHOICE)"}, f, indent=2)
print("wrote dh_ftheta_float.json", flush=True)
