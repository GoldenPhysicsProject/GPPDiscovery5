"""Empirical look-elsewhere control with decoy targets.
For each grammar, build the set of DISTINCT values (equal rewrites, e.g. golden identities
phi^-1 = phi^-2 + phi^-3, count once), then count how many land within eps of
  (i) the real target T = 137.035999177 and
  (ii) 2000 decoy targets drawn uniformly from [T-0.01, T+0.01] (excluding |d|<1e-5).
The decoy mean is the honest chance rate: it captures the non-uniform, hierarchical clustering
that tiny correction terms create (a uniform-density estimate understates it).
Daniel Toupin / Golden Physics Project. Seat: Heart of the Cosmos / Grok (GPPDiscovery5)."""
import numpy as np, itertools, math, json
from fractions import Fraction
from mpmath import zeta, catalan, euler, nsum, inf
T = 137.035999177; phi = (1 + 5 ** .5) / 2
rng = np.random.default_rng(20261121)
dec = T + rng.uniform(-0.01, 0.01, 4000); dec = dec[np.abs(dec - T) > 1e-5][:2000]
EPS = {"1sigma 2.1e-8": 2.1e-8, "Cs-Rb 1.6e-7": 1.6e-7}

def distinct(v):
    v = v[np.isfinite(v)]
    return np.unique(np.round(v, 11))

def hits(sorted_v, t, eps):
    return np.searchsorted(sorted_v, t + eps) - np.searchsorted(sorted_v, t - eps)

fam = {}
# A: 360phi^-2 + a phi^-n + third term (as in search.py)
t2 = np.array([s * a * phi ** -n for s in (1, -1) for a in range(1, 10) for n in range(1, 10)])
t3 = np.array([s * b * phi ** -m for s in (1, -1) for b in range(1, 10) for m in range(1, 10)]
              + [s * (b * phi) ** -m for s in (1, -1) for b in range(2, 10) for m in range(1, 10)])
fam["A golden angle + 2 golden terms"] = distinct((360 * phi ** -2 + t2[:, None] + t3[None, :]).ravel())
# A4: golden angle + 3 distinct terms from t3
idx = np.array(list(itertools.combinations(range(len(t3)), 3)))
fam["A4 golden angle + 3 golden terms"] = distinct(360 * phi ** -2 + t3[idx].sum(axis=1))
del idx
# B: GGC grammar restricted to n_i <= 400 (6.4e7 values) for memory
n = np.arange(1, 401, dtype=float); base = 360 * phi ** -2
S12 = (1 - (np.pi / n * phi ** -12)[:, None] - (np.pi / n * phi ** -16)[None, :]).ravel()
vB = (base * (S12[:, None] - (np.pi * phi ** -20 / n)[None, :])).ravel()
vB = vB[np.abs(vB - T) < 0.02]
fam["B GGC grammar (pi/n_i, n_i<=400)"] = distinct(vB)
# C, D: framework monomials
L2m3 = float(nsum(lambda k: 1 / (3 * k + 1) ** 2 - 1 / (3 * k + 2) ** 2, [0, inf]))
K = [math.pi, phi, math.e, 2 ** .5, 3 ** .5, 5 ** .5, 7 ** .5, float(zeta(3)), float(zeta(5)), float(catalan),
     math.log(2), float(euler), L2m3]
exps = [-4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4]
mon = []
for i, j in itertools.combinations_with_replacement(range(len(K)), 2):
    for a in exps:
        for b in exps:
            if i == j and b != exps[0]: continue
            mon.append(K[i] ** a * (K[j] ** b if i != j else 1.0))
mon = np.array(mon)
R = np.array(sorted({float(Fraction(p, q)) for p in range(1, 401) for q in range(1, 13)}))
vC = (R[:, None] * mon[None, :]).ravel(); vC = vC[np.abs(vC - T) < 0.02]
fam["C r*x^a*y^b (13 constants)"] = distinct(vC)
R2 = np.array(sorted({float(Fraction(p, q)) for p in range(1, 13) for q in range(1, 13)}))
vD = (137 + np.concatenate([R2, -R2])[:, None] * mon[None, :]).ravel(); vD = vD[np.abs(vD - T) < 0.02]
fam["D 137 +/- r*x^a*y^b"] = distinct(vD)

out = []
for name, v in fam.items():
    row = {"family": name, "distinct_values_near_T(|d|<0.02)": int(np.count_nonzero(np.abs(v - T) < 0.02))}
    for en, eps in EPS.items():
        real = int(hits(v, T, eps)); d = np.array([hits(v, t, eps) for t in dec])
        lam = d.mean(); frac = np.mean(d >= max(real, 1))
        row[en] = {"real_target_hits": real, "decoy_mean_hits": float(lam),
                   "decoy_frac_with_>=real_hits": float(frac),
                   "poisson_P(>=real)": float(1 - sum(math.exp(-lam) * lam ** k / math.factorial(k) for k in range(real))) if real else 1.0}
    out.append(row)
    print(json.dumps(row))
json.dump(out, open("control_targets.json", "w"), indent=1)
