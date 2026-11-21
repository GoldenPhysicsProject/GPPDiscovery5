"""Family A4: golden angle + THREE golden correction terms (one more term than Sherbon/Pellis).
Shows how fast 1-sigma hits appear once the grammar grows by one term.
Daniel Toupin / Golden Physics Project. Seat: Heart of the Cosmos / Grok (GPPDiscovery5)."""
import numpy as np, math, json
T = 137.035999177; phi = (1 + 5 ** .5) / 2
terms = [(s * a * phi ** -n, f"{'+' if s>0 else '-'}{a}phi^-{n}") for s in (1, -1) for a in range(1, 10) for n in range(1, 10)]
terms += [(s * (b * phi) ** -m, f"{'+' if s>0 else '-'}({b}phi)^-{m}") for s in (1, -1) for b in range(2, 10) for m in range(1, 10)]
v = np.array([x for x, _ in terms]); L = [l for _, l in terms]; n = len(v)
base = 360 * phi ** -2
# unordered triples i<j<k: for each (i,j) find k>j with |base+vi+vj+vk-T|<eps using a sorted array
order = np.argsort(v); vs = v[order]
res = {}
for name, eps in {"1sigma (2.1e-8)": 2.1e-8, "Cs-Rb spread (1.6e-7)": 1.6e-7, "window 0.1": 0.1}.items():
    cnt = 0; hits = []
    for i in range(n):
        for j in range(i + 1, n):
            need = T - base - v[i] - v[j]
            lo = np.searchsorted(vs, need - eps); hi = np.searchsorted(vs, need + eps)
            for kk in range(lo, hi):
                k = order[kk]
                if k > j:
                    cnt += 1
                    if eps < 1e-6: hits.append((L[i], L[j], L[k], base + v[i] + v[j] + v[k] - T))
    res[name] = (cnt, hits)
N = n * (n - 1) * (n - 2) // 6
local = res["window 0.1"][0]
out = {"family": "A4 golden angle + 3 golden terms", "N": N, "N_in_window_0.2": local, "tols": {}}
for name in ("1sigma (2.1e-8)", "Cs-Rb spread (1.6e-7)"):
    eps = 2.1e-8 if name.startswith("1sigma") else 1.6e-7
    E = local * 2 * eps / 0.2
    distinct = sorted({round(h[3], 11) for h in res[name][1]})
    out["tols"][name] = {"observed": res[name][0], "distinct_values": distinct, "expected_by_chance": E, "hits": [list(map(str, h)) for h in res[name][1][:40]]}
    print(name, "observed", res[name][0], "distinct values", distinct, "expected", round(E, 2))
    for h in res[name][1][:15]: print("   360phi^-2", *h)
print("N", N, "local", local)
json.dump(out, open("search_golden4.json", "w"), indent=1)
