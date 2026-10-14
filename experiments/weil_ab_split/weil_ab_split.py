"""Weil form A/B split at F_theta failure windows. NUMERICAL (float).

Reuse experiments/weil_small_support/weil_ritz.py normalization exactly.
Q = A + B, where
  A = archimedean/gamma digamma integral
  B = pole terms - prime (von Mangoldt) terms
For F_theta: sum over sig in {+theta,-theta} of the zeta-type formula on F_sig.

At each h: Ritz-minimize total Q on even (and odd) Legendre-bump basis,
then report the signed Rayleigh quotients of A and B on that eigenvector
(and optionally pole-only / primes-only inside B).

Question: at h where F_{1/4} goes negative while zeta stays +, does F_theta
fail because B gets worse, A gets worse, or both?

RH remains OPEN. Ritz negatives are float witnesses; positives are not proofs.
Run: python weil_ab_split.py [M] [hmin] [hmax] [dh]
"""
from __future__ import annotations
import json, os, sys
import numpy as np
from numpy.polynomial import legendre as L
from scipy.special import digamma
from scipy.linalg import eigh

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from common import write_result

M = int(sys.argv[1]) if len(sys.argv) > 1 else 18
hs = (
    np.arange(float(sys.argv[2]), float(sys.argv[3]) + 1e-9, float(sys.argv[4]))
    if len(sys.argv) > 4
    else np.arange(0.40, 0.521, 0.01)
)
xq0, wq0 = L.leggauss(1200)


def vm(n: int) -> float:
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return float(np.log(p)) if m == 1 else 0.0
    return 0.0


def basis(x, h):
    u = x / h
    B = np.array([L.legval(u, [0] * j + [1]) for j in range(M)]) * (1 - u**2) ** 4
    B[:, np.abs(u) >= 1] = 0
    return B


def parts(h: float, theta: float = 0.0):
    """Return (QA, QB_pole, QB_prime, gram). Total Q = QA + QB_pole + QB_prime.
    QB_prime already carries the minus sign (negative definite-looking contribution).
    """
    xq, wq = xq0 * h, wq0 * h
    B = basis(xq, h)
    gram = (B * wq) @ B.T
    R = 600.0 / h
    r = np.linspace(0, R, 30001)
    dr = r[1] - r[0]
    wr = np.full_like(r, dr)
    wr[0] = wr[-1] = dr / 2
    E = lambda z: (B * wq) @ np.exp(np.outer(xq, z))

    def Fmat(u):
        Bs = basis(xq - u, h)
        A = (B * wq) @ Bs.T
        return (A + A.T) / 2

    QA = np.zeros((M, M))
    QBp = np.zeros((M, M))  # pole
    QBpr = np.zeros((M, M))  # primes (with minus already)
    w = digamma(0.25 + 0.5j * r).real - np.log(np.pi)
    sigs = [0.0] if theta == 0.0 else [theta, -theta]
    for sg in sigs:
        Ep, Em = E(1j * r + sg), E(-1j * r - sg)
        Ablk = (Ep * (w * wr)) @ Em.T
        QA += 2 * np.real((Ablk + Ablk.T) / 2) / (2 * np.pi)
        for a in (0.5, -0.5):
            e1 = E(np.array([a + sg]))[:, 0]
            e2 = E(np.array([-a - sg]))[:, 0]
            P = np.outer(e1, e2).real
            QBp += (P + P.T) / 2
        for n in range(2, int(np.exp(2 * h)) + 1):
            lv = vm(n)
            if lv:
                QBpr -= lv / np.sqrt(n) * (n**sg + n ** (-sg)) * Fmat(np.log(n))
    return QA, QBp, QBpr, gram


def rayleigh(Q, v, G, vGv=None):
    if vGv is None:
        vGv = float(v @ G @ v)
    return float(v @ Q @ v) / vGv


def analyze(h: float, theta: float, label: str):
    QA, QBp, QBpr, G = parts(h, theta)
    QB = QBp + QBpr
    QT = QA + QB
    out = {"h": round(float(h), 4), "label": label, "theta": theta}
    for sector, idx in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        Ix = np.ix_(idx, idx)
        # Ritz on total
        evals, evecs = eigh(QT[Ix], G[Ix])
        i_min = int(np.argmin(evals))
        v = evecs[:, i_min]
        # pad back to full M
        vfull = np.zeros(M)
        vfull[idx] = v
        vGv = float(vfull @ G @ vfull)
        a = rayleigh(QA, vfull, G, vGv)
        bp = rayleigh(QBp, vfull, G, vGv)
        bpr = rayleigh(QBpr, vfull, G, vGv)
        b = bp + bpr
        tot = a + b
        # also lambda_min of A and of B alone (diagnostic)
        lamA = float(eigh(QA[Ix], G[Ix], eigvals_only=True)[0])
        lamB = float(eigh(QB[Ix], G[Ix], eigvals_only=True)[0])
        out[sector] = {
            "lam_total": float(evals[i_min]),
            "A_on_min": a,
            "B_on_min": b,
            "B_pole_on_min": bp,
            "B_primes_on_min": bpr,
            "check_A_plus_B": tot,
            "lam_min_A_alone": lamA,
            "lam_min_B_alone": lamB,
        }
    return out


cases = [
    ("zeta", 0.0),
    ("F_1/4", 0.25),
]
# optional F_1/20 for context near the window
if os.environ.get("EXTRA_F120", "0") == "1":
    cases.append(("F_1/20", 0.05))

results = {nm: [] for nm, _ in cases}
for h in hs:
    for nm, th in cases:
        row = analyze(float(h), th, nm)
        results[nm].append(row)
        e = row["even"]
        print(
            f"{nm:7s} h={h:.3f} even: tot={e['lam_total']:+.6e} "
            f"A={e['A_on_min']:+.6e} B={e['B_on_min']:+.6e} "
            f"(pole={e['B_pole_on_min']:+.6e} primes={e['B_primes_on_min']:+.6e})",
            flush=True,
        )

# Cross-eval: at each h, take F_1/4 even minimizer and evaluate zeta's A,B on it, and vice versa.
# Rebuild once per h.
cross = []
for h in hs:
    QA_z, QBp_z, QBpr_z, G = parts(float(h), 0.0)
    QB_z = QBp_z + QBpr_z
    QT_z = QA_z + QB_z
    QA_f, QBp_f, QBpr_f, _ = parts(float(h), 0.25)
    QB_f = QBp_f + QBpr_f
    QT_f = QA_f + QB_f
    ev = list(range(0, M, 2))
    Ix = np.ix_(ev, ev)
    # F_theta minimizer
    evals_f, evecs_f = eigh(QT_f[Ix], G[Ix])
    vf = np.zeros(M)
    vf[ev] = evecs_f[:, int(np.argmin(evals_f))]
    vGv = float(vf @ G @ vf)
    # zeta minimizer
    evals_z, evecs_z = eigh(QT_z[Ix], G[Ix])
    vz = np.zeros(M)
    vz[ev] = evecs_z[:, int(np.argmin(evals_z))]
    vGv_z = float(vz @ G @ vz)
    cross.append(
        {
            "h": round(float(h), 4),
            "on_F14_min": {
                "Q_F14": rayleigh(QT_f, vf, G, vGv),
                "A_F14": rayleigh(QA_f, vf, G, vGv),
                "B_F14": rayleigh(QB_f, vf, G, vGv),
                "Q_zeta": rayleigh(QT_z, vf, G, vGv),
                "A_zeta": rayleigh(QA_z, vf, G, vGv),
                "B_zeta": rayleigh(QB_z, vf, G, vGv),
            },
            "on_zeta_min": {
                "Q_zeta": rayleigh(QT_z, vz, G, vGv_z),
                "A_zeta": rayleigh(QA_z, vz, G, vGv_z),
                "B_zeta": rayleigh(QB_z, vz, G, vGv_z),
                "Q_F14": rayleigh(QT_f, vz, G, vGv_z),
                "A_F14": rayleigh(QA_f, vz, G, vGv_z),
                "B_F14": rayleigh(QB_f, vz, G, vGv_z),
            },
        }
    )
    c = cross[-1]["on_F14_min"]
    print(
        f"CROSS  h={h:.3f} on F14-min: "
        f"A_F={c['A_F14']:+.4e} B_F={c['B_F14']:+.4e} | "
        f"A_z={c['A_zeta']:+.4e} B_z={c['B_zeta']:+.4e} "
        f"dA={c['A_F14']-c['A_zeta']:+.4e} dB={c['B_F14']-c['B_zeta']:+.4e}",
        flush=True,
    )

tag = f"M{M}_{hs[0]:.3f}_{hs[-1]:.3f}"
out_path = os.path.join(os.path.dirname(__file__), f"result_{tag}.json")
write_result(
    out_path,
    {
        "experiment": "weil_ab_split",
        "label": "NUMERICAL (Ritz; A=arch/gamma, B=pole-primes)",
        "M": M,
        "split": {
            "A": "digamma integral (same w as weil_ritz.py; summed over sig for F_theta)",
            "B": "pole outer-products minus von Mangoldt F(log n) terms",
            "normalization": "identical to experiments/weil_small_support/weil_ritz.py",
        },
        "own_minimizer": results,
        "cross_on_shared_testfn": cross,
        "note": (
            "lam_total < 0 is a float witness of failure for that L-function on this window. "
            "Positive lam_total does NOT prove positivity. RH OPEN."
        ),
    },
    kind="float",
    method="Legendre-bump Ritz A/B split; Gauss-Legendre 1200 + trapezoid in r",
)
print(f"wrote {out_path}", flush=True)
