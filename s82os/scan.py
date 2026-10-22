"""S82-OS scan (NUMERICAL, float64): OS/Hankel form B vs Weil parity sectors T+-B, for ZETA and DH; archimedean sign;
UV limit of B_min (Carleman); Perron vector of the prime operator vs the Weil near-null vector. RH OPEN."""
import json, math, sys
import numpy as np
from os_hankel import build
from numpy.polynomial import legendre as leg
out = []
# identity check: k1+k2 == e^{s/2}/sinh s == 2 sum_k e^{-(2k+1/2)s}
ss = np.linspace(0.05, 3, 40)
k12 = np.exp(-ss) / np.sinh(ss) + 2 * np.exp(-ss / 4) * np.sinh(3 * ss / 4) / np.sinh(ss)
cm = 2 * sum(np.exp(-(2 * k + 0.5) * ss) for k in range(400))
out.append({"check": "k1+k2 = e^{s/2}/sinh s = 2 sum_k e^{-(2k+1/2)s}", "max_rel_err_1": float(np.max(abs(k12 / (np.exp(ss / 2) / np.sinh(ss)) - 1))),
            "max_rel_err_2": float(np.max(abs(k12 / cm - 1)))})
def basisfun(L, N, x):
    t = 2 * x / L - 1
    return np.array([math.sqrt((2 * k + 1) / L) * leg.legval(t, [0] * k + [1]) for k in range(N)])
for FORM in ("ZETA", "DH"):
    for L in (0.2, 0.3, 0.5, 0.7, 0.85, 1.0):
        for N in ((24, 40) if L in (0.5, 1.0) else (24,)):
            T, B = build(FORM, L, N); c = build.comp
            eB = np.linalg.eigvalsh(B); eBa = np.linalg.eigvalsh(c['B_arch']); eBp = np.linalg.eigvalsh(c['B_prime']) if np.any(c['B_prime']) else np.array([0.0])
            eP, VP = np.linalg.eigh(T + B); eM = np.linalg.eigvalsh(T - B)
            r = {"kind": "float", "form": FORM, "L": L, "N": N, "B_min": eB[0], "B_max": eB[-1],
                 "B_neg": int((eB < -1e-8).sum()), "B_pos": int((eB > 1e-8).sum()),
                 "B_arch_max(<=0 expected)": eBa[-1], "B_arch_min": eBa[0],
                 "B_prime_min": eBp[0], "B_prime_max": eBp[-1],
                 "even_min": eP[0], "odd_min": eM[0]}
            eT = np.linalg.eigvalsh(T)
            if eT[0] > 0:
                Li = np.linalg.inv(np.linalg.cholesky(T)); eC = np.linalg.eigvalsh(Li @ B @ Li.T)
                r.update({"C_min": eC[0], "C_max": eC[-1]})
            # Perron analysis (even sector prime operator P_even = T_prime + B_prime, positivity preserving for ZETA)
            Pe = c['T_prime'] + c['B_prime']
            if np.any(Pe):
                lp, V = np.linalg.eigh(0.5 * (Pe + Pe.T)); vP = V[:, -1]
                xs = np.linspace(0, L, 401); fP = vP @ basisfun(L, N, xs); fP *= np.sign(fP.sum())
                v0 = VP[:, 0]; ap = c['ap'] + c['am']   # even pole vector  ~ cosh(u/2)
                r.update({"P_even_top": lp[-1], "P_even_bottom": lp[0], "Perron_fn_min/max": float(fP.min() / fP.max()),
                          "overlap(nullvec,Perron)": float(abs(v0 @ vP)), "overlap(nullvec,cosh/|cosh|)": float(abs(v0 @ ap) / np.linalg.norm(ap)),
                          "overlap(Perron,cosh/|cosh|)": float(abs(vP @ ap) / np.linalg.norm(ap))})
            r = {k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in r.items()}
            out.append(r); print(json.dumps(r), flush=True)
json.dump(out, open("scan_results.json", "w"), indent=1)
