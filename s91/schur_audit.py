"""Grok S91 Task B: audit of Turbo's wave-227 per-parity Schur scalars ('L=3': S_even = 4.22e-10, S_odd = 1.24e-11) in OUR normalisation.
Window x in [e^-L, e^L] (S87 builder, primes n <= e^{2L}; Connes-Consani mu = e^{2L}, L_CC = 2L). R = A - P (pole-free form).
 even: Q_e = R + p p^T,  t_e = p^T R^-1 p,  s_e = 1 + t_e   (Q_e >= 0 <=> s_e <= 0, given inertia(R_e) = (N-1, 1));
 odd : Q_o = R - q q^T,  t_o = q^T R^-1 q,  s_o = 1 - t_o   (Q_o >= 0 <=> s_o >= 0, given R_o > 0).
Relative cancellation rel = |s| / (1 + |t|). Blocks: prime block P (enters as -P), boundary block B = A +/- pole (arch + pole).
mpmath dps 90, N = 24 (and N = 32 at 0.6, 1.0, 1.2, 1.4 for basis dependence). HIGH-PRECISION NUMERICAL, not interval. RH OPEN."""
import json, math, mpmath as mp
import forms as FM
def inertia(M, tol):
    E = mp.eigsy(M, eigvals_only=True); return [sum(1 for x in E if x < -tol), sum(1 for x in E if x > tol)], E
rows = []
L3 = math.log(3) / 2      # mu = lambda^2 = 3 (CC's variable); also L = log 3 (lambda = 3)
cases = [(L3, 24), (L3, 32), (math.log(3), 24), (math.log(3), 32)] + [(L, 24) for L in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)] + [(L, 32) for L in (0.6, 1.0, 1.2, 1.4)]
for L, N in cases:
    d = FM.load("ZETA", 0, L, N); r = {"L": L, "N": N, "mu=e^2L": math.exp(2 * L)}
    for par in "eo":
        A = d["A" + par]; P = d["P" + par]; v = d["p" + par][0]; sg = 1 if par == "e" else -1
        R = A - P; Q = R + sg * v * v.T; B = A + sg * v * v.T
        y = mp.lu_solve(R, v); t = (v.T * y)[0]; s = 1 + t if par == "e" else 1 - t
        EQ = mp.eigsy(Q, eigvals_only=True); qmin, qmax = min(EQ), max(EQ)
        tol = mp.mpf(10) ** -60
        iP, _ = inertia(P, tol); iB, _ = inertia(B, tol); iR, _ = inertia(R, tol); iA, _ = inertia(A, tol)
        r[par] = {"s": mp.nstr(s, 8), "t=pole'R^-1pole": mp.nstr(t, 12), "rel=|s|/(1+|t|)": mp.nstr(abs(s) / (1 + abs(t)), 6),
                  "qmin": mp.nstr(qmin, 6), "qmin/||Q||": mp.nstr(qmin / qmax, 6), "|y|^2": mp.nstr((y.T * y)[0], 6),
                  "first_order -s/|y|^2 (even) or s/|y|^2 (odd)": mp.nstr(-sg * s / (y.T * y)[0], 6),
                  "inertia(-,+) P (prime block)": iP, "inertia B=A+/-pole (boundary)": iB, "inertia R=A-P": iR, "inertia A": iA,
                  "P_indefinite": iP[0] > 0 and iP[1] > 0, "B_indefinite": iB[0] > 0 and iB[1] > 0}
    rows.append(r)
    print(L, N, "| E s=%s rel=%s qmin=%s q/|Q|=%s P%s B%s | O s=%s rel=%s qmin=%s q/|Q|=%s P%s B%s" % tuple(
        [r["e"][k] for k in ("s", "rel=|s|/(1+|t|)", "qmin", "qmin/||Q||", "inertia(-,+) P (prime block)", "inertia B=A+/-pole (boundary)")] +
        [r["o"][k] for k in ("s", "rel=|s|/(1+|t|)", "qmin", "qmin/||Q||", "inertia(-,+) P (prime block)", "inertia B=A+/-pole (boundary)")]), flush=True)
# scaling fits on N=24 rows, L in [0.6, 1.4]; and log-linear interpolation of the L where rel equals Turbo's numbers
import numpy as np
fit = {}
for par in "eo":
    R24 = sorted([r for r in rows if r["N"] == 24 and 0.5 <= r["L"] <= 1.5], key=lambda r: r["L"])
    Lv = np.array([r["L"] for r in R24]); muv = np.exp(2 * Lv)
    lg = np.array([math.log10(float(mp.mpf(r[par]["rel=|s|/(1+|t|)"]))) for r in R24])
    m = (Lv >= 0.6) & (Lv <= 1.4)
    bL = np.polyfit(Lv[m], lg[m], 1); bmu = np.polyfit(muv[m], lg[m], 1)
    tgt = math.log10(4.22e-10 if par == "e" else 1.24e-11)
    k = next((i for i in range(len(lg) - 1) if lg[i] >= tgt >= lg[i + 1]), None)
    Lx = None if k is None else float(Lv[k] + (tgt - lg[k]) * (Lv[k + 1] - Lv[k]) / (lg[k + 1] - lg[k]))
    fit[par] = {"log10 rel ~ a + b L, N=24, L 0.6-1.4 (b,a)": bL.tolist(), "log10 rel ~ a + b mu (b,a)": bmu.tolist(),
                "decades per 0.1 in L (N=24)": float(-bL[0] / 10), "L_ours where N=24 rel equals Turbo value (interp)": Lx,
                "-> mu = e^{2L}": None if Lx is None else math.exp(2 * Lx)}
print(json.dumps(fit, indent=1))
json.dump({"rows": rows, "fits": fit, "turbo_claim": {"L": 3, "S_even": 4.22e-10, "S_odd": 1.24e-11}}, open("schur_audit.json", "w"), indent=1)
