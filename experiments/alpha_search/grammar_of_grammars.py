"""Global look-elsewhere: many sibling grammars of the same shape as Sherbon/Pellis.
Lead L in {360/phi^2, 137, 2pi*... no}: L + s1*a*c2^-n + s2*t3, where
  c2, c3 in {phi, pi, e, sqrt2, sqrt3, sqrt5}, t3 in {b*c3^-m (b=1..9)} or {(b*c3)^-m (b=2..9)}, a,n,m in 1..9.
Each (lead, c2, c3, shape) is one grammar of ~5e4 expressions. We count, per grammar, distinct
values within eps of T, and the decoy chance rate per grammar.
Daniel Toupin / Golden Physics Project. Seat: Heart of the Cosmos / Grok (GPPDiscovery5)."""
import numpy as np, math, json
T = 137.035999177; phi = (1 + 5 ** .5) / 2
C = {"phi": phi, "pi": math.pi, "e": math.e, "sqrt2": 2 ** .5, "sqrt3": 3 ** .5, "sqrt5": 5 ** .5}
LEADS = {"360/phi^2": 360 * phi ** -2, "137": 137.0, "4pi^3+pi^2+pi": 4 * math.pi ** 3 + math.pi ** 2 + math.pi}
rng = np.random.default_rng(11)
dec = T + rng.uniform(-0.01, 0.01, 20000); dec = dec[np.abs(dec - T) > 1e-5]
eps = 2.1e-8
rows = []; total_hits = 0; total_rate = 0.0
for ln, L in LEADS.items():
    for c2n, c2 in C.items():
        t2 = np.array([s * a * c2 ** -n for s in (1, -1) for a in range(1, 10) for n in range(1, 10)])
        lab2 = [f"{'+' if s>0 else '-'}{a}*{c2n}^-{n}" for s in (1, -1) for a in range(1, 10) for n in range(1, 10)]
        for c3n, c3 in C.items():
            for shape in ("b*c^-m", "(b*c)^-m"):
                if shape == "b*c^-m":
                    t3 = [(s * b * c3 ** -m, f"{'+' if s>0 else '-'}{b}*{c3n}^-{m}") for s in (1, -1) for b in range(1, 10) for m in range(1, 10)]
                else:
                    t3 = [(s * (b * c3) ** -m, f"{'+' if s>0 else '-'}({b}*{c3n})^-{m}") for s in (1, -1) for b in range(2, 10) for m in range(1, 10)]
                a3 = np.array([x for x, _ in t3])
                V = L + t2[:, None] + a3[None, :]
                flat = V.ravel(); sel = np.nonzero(np.abs(flat - T) < eps)[0]
                hitvals = {}
                for k in sel:
                    hitvals.setdefault(round(float(flat[k]), 11), f"{ln} {lab2[k // len(t3)]} {t3[k % len(t3)][1]}")
                u = np.unique(np.round(flat, 11))
                h = np.searchsorted(u, dec + eps) - np.searchsorted(u, dec - eps)
                rate = float((h > 0).mean())
                total_hits += len(hitvals); total_rate += rate
                rows.append({"lead": ln, "c2": c2n, "c3": c3n, "shape": shape, "N": int(flat.size),
                             "hits_1sigma": list(hitvals.values()), "chance_rate": rate})
print("grammars:", len(rows), " total distinct 1-sigma hits:", total_hits, " expected by chance (sum of rates):", round(total_rate, 3))
for r in rows:
    if r["hits_1sigma"]: print(r["lead"], r["c2"], r["c3"], r["shape"], r["hits_1sigma"], "rate", r["chance_rate"])
json.dump({"grammars": len(rows), "total_hits": total_hits, "expected": total_rate, "rows": rows}, open("grammar_of_grammars.json", "w"), indent=1)
