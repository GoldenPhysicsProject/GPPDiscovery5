"""Look-elsewhere search: how many expressions of comparable complexity hit 1/alpha?
Daniel Toupin / Golden Physics Project. Seat: Heart of the Cosmos / Grok (GPPDiscovery5).

For each grammar (family) we enumerate every expression, count exact hits within
tolerance eps of T = 137.035999177 (CODATA 2022), and compare with the number expected
by chance: E = N_local * (2 eps / W), where N_local is the number of family values in the
window [T - W/2, T + W/2] (W = 0.2). A hit is only evidence if observed >> expected,
i.e. P(>=1 chance hit) = 1 - exp(-E) is small. Float64 is ample (eps >= 1e-9 absolute);
all reported hits are re-verified at 50 digits with mpmath.
"""
import numpy as np, itertools, json, math
from fractions import Fraction
from mpmath import mp, mpf
mp.dps = 50
T = 137.035999177
TOLS = {"1sigma_CODATA (2.1e-8)": 2.1e-8, "Cs-Rb spread (1.6e-7)": 1.6e-7, "1 ppm (1.37e-4)": 1.37e-4}
W = 0.2
phi = (1 + 5 ** 0.5) / 2

def stats(vals, labels_fn, name, size_note):
    vals = np.asarray(vals)
    N = vals.size
    local = int(np.count_nonzero(np.abs(vals - T) < W / 2))
    res = {"family": name, "N": int(N), "note": size_note, "N_in_window_0.2": local, "tols": {}}
    order = np.argsort(np.abs(vals - T))
    for tn, eps in TOLS.items():
        hits = int(np.count_nonzero(np.abs(vals - T) < eps))
        E = local * 2 * eps / W
        res["tols"][tn] = {"observed": hits, "expected_by_chance": E, "P_at_least_one": 1 - math.exp(-E)}
    res["best"] = [{"expr": labels_fn(int(i)), "value": float(vals[i]), "diff": float(vals[i] - T)} for i in order[:8]]
    return res

out = []

# ---- Family A: golden-angle 3-term grammar (Heyrovska/Sherbon/Pellis shape) ----
# 360 phi^-2 + s1 * a phi^-n + s2 * t3, t3 in {b phi^-m} U {(b phi)^-m}, a,b in 1..9, n,m in 1..9
t2 = [(s * a * phi ** -n, f"{'+' if s>0 else '-'}{a}phi^-{n}") for s in (1, -1) for a in range(1, 10) for n in range(1, 10)]
t3 = [(s * b * phi ** -m, f"{'+' if s>0 else '-'}{b}phi^-{m}") for s in (1, -1) for b in range(1, 10) for m in range(1, 10)]
t3 += [(s * (b * phi) ** -m, f"{'+' if s>0 else '-'}({b}phi)^-{m}") for s in (1, -1) for b in range(2, 10) for m in range(1, 10)]
A2 = np.array([x for x, _ in t2]); A3 = np.array([x for x, _ in t3])
valsA = (360 * phi ** -2 + A2[:, None] + A3[None, :]).ravel()
labA = lambda i: f"360phi^-2 {t2[i // len(t3)][1]} {t3[i % len(t3)][1]}"
out.append(stats(valsA, labA, "A golden-angle + 2 golden terms", f"{len(t2)} x {len(t3)}"))

# ---- Family B: GGC grammar 360/phi^2 (1 - pi/n1 phi^-12 - pi/n2 phi^-16 - pi/n3 phi^-20), n_i in 1..1000 ----
# counted analytically per (n1,n2): values are monotone in n3, so count n3 with |v-T|<eps directly.
n = np.arange(1, 1001, dtype=float)
base = 360 * phi ** -2
c1 = np.pi / n * phi ** -12; c2 = np.pi / n * phi ** -16; k3 = np.pi * phi ** -20
NB = 1000 ** 3
resB = {"family": "B GGC 3-correction grammar (pi/n_i, n_i<=1000)", "N": NB, "note": "1000^3", "tols": {}}
S12 = (1 - c1[:, None] - c2[None, :])  # shape 1000x1000
# value(n3) = base*(S12 - k3/n3); need n3 in [k3/(S12 - (T-eps)/base) , k3/(S12-(T+eps)/base)]
def count_B(eps):
    lo_den = S12 - (T - eps) / base; hi_den = S12 - (T + eps) / base
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(lo_den > 0, k3 / lo_den, np.inf); b = np.where(hi_den > 0, k3 / hi_den, np.inf)
    lo = np.minimum(a, b); hi = np.maximum(a, b)
    lo = np.clip(np.ceil(lo), 1, 1001); hi = np.clip(np.floor(hi), 0, 1000)
    return int(np.sum(np.maximum(0, hi - lo + 1)))
localB = count_B(W / 2)
for tn, eps in TOLS.items():
    E = localB * 2 * eps / W
    resB["tols"][tn] = {"observed": count_B(eps), "expected_by_chance": E, "P_at_least_one": 1 - math.exp(-E)}
resB["N_in_window_0.2"] = localB
# list the 1-sigma hits explicitly
hitsB = []
for i1 in range(1000):
    need = k3 / (S12[i1] - T / base)
    for i2 in np.nonzero((need > 0.5) & (need < 1000.5))[0]:
        n3 = int(round(need[i2])); v = base * (S12[i1, i2] - k3 / n3)
        if abs(v - T) < 2.1e-8: hitsB.append({"n1": i1 + 1, "n2": int(i2) + 1, "n3": n3, "diff": v - T})
resB["best"] = sorted(hitsB, key=lambda h: abs(h["diff"]))[:12]
resB["n_1sigma_hits_listed"] = len(hitsB)
out.append(resB)

# ---- Family C: Daniel-framework monomials r * x^a * y^b ----
from mpmath import zeta, catalan, euler, log, pi as mpi, sqrt as msqrt, nsum, inf
L2m3 = float(nsum(lambda k: 1 / (3 * k + 1) ** 2 - 1 / (3 * k + 2) ** 2, [0, inf]))  # L(2, chi_-3)
K = {"pi": math.pi, "phi": phi, "e": math.e, "sqrt2": 2 ** .5, "sqrt3": 3 ** .5, "sqrt5": 5 ** .5, "sqrt7": 7 ** .5,
     "zeta3": float(zeta(3)), "zeta5": float(zeta(5)), "G": float(catalan), "ln2": math.log(2), "gammaE": float(euler),
     "L2chi-3": L2m3}
names = list(K); vals = np.array([K[k] for k in names])
exps = [-4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4]
rats = sorted({Fraction(p, q) for p in range(1, 401) for q in range(1, 13)})
R = np.array([float(r) for r in rats])
mon = []; monlab = []
for i, j in itertools.combinations_with_replacement(range(len(names)), 2):
    for a in exps:
        for b in exps:
            if i == j and b != exps[0]: continue  # single-constant case once
            v = vals[i] ** a * (vals[j] ** b if i != j else 1.0)
            mon.append(v); monlab.append(f"{names[i]}^{a}" + (f"*{names[j]}^{b}" if i != j else ""))
mon = np.array(mon)
valsC = (R[:, None] * mon[None, :]).ravel()
labC = lambda k: f"({rats[k // len(mon)]}) * {monlab[k % len(mon)]}"
out.append(stats(valsC, labC, "C r * x^a * y^b, 13 framework constants, r=p/q p<=400 q<=12",
                 f"{len(rats)} rationals x {len(mon)} monomials"))

# ---- Family D: 137 + r * x^a * y^b (the 'anomaly' 0.035999177) ----
rats2 = sorted({Fraction(p, q) for p in range(1, 13) for q in range(1, 13)})
R2 = np.array([float(r) for r in rats2])
valsD = (137 + np.concatenate([R2, -R2])[:, None] * mon[None, :]).ravel()
signs = [("+", r) for r in rats2] + [("-", r) for r in rats2]
labD = lambda k: f"137 {signs[k // len(mon)][0]} ({signs[k // len(mon)][1]}) * {monlab[k % len(mon)]}"
out.append(stats(valsD, labD, "D 137 +/- r * x^a * y^b, r=p/q p,q<=12", f"{2*len(rats2)} x {len(mon)}"))

# verify best hits at 50 digits for family A
from mpmath import sqrt as S
PHI = (1 + S(5)) / 2
def evalA(lbl):
    return eval(lbl.replace("phi^", "PHI**").replace(")^", ")**").replace("phi", "*PHI").replace("(", "(").replace("360*PHI", "360*PHI"),
                {"PHI": PHI})
json.dump(out, open("search_results.json", "w"), indent=1, default=float)
for r in out:
    print(f"\n== {r['family']}  N={r['N']:.3g} ({r['note']}), values within +-0.1: {r['N_in_window_0.2']}")
    for tn, t in r["tols"].items():
        print(f"   {tn:26s} observed={t['observed']:8d}  expected={t['expected_by_chance']:.3g}  P(>=1 by chance)={t['P_at_least_one']:.3g}")
    for b in r["best"][:8]:
        print("     ", b)
