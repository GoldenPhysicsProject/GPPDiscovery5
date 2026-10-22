"""S82-OS: node count and pole-channel overlap of existing Arb-certified Ritz witnesses (s75). Shape analysis is NUMERICAL (float);
the Rayleigh quotients themselves are the CERTIFIED s75 values. f(u) = sum_m v_m cos((2m+1) pi u / 2L) on [-L, L] (even)."""
import json, numpy as np
from fractions import Fraction
res = []
for fn in ["cosvec_DH_even_L3421_2000_K300.json", "cosvec_DH_even_L137_80_K129.json", "cosvec_ZETA_even_L137_80_K129.json"]:
    d = json.load(open("../s75/" + fn)); L = float(Fraction(d["L"].replace('_', '/')))
    v = np.array([float(Fraction(q)) for q in d["vec"]]); K = len(v)
    a = (2 * np.arange(K) + 1) * np.pi / (2 * L); xs = np.linspace(0, L, 20001)
    w = np.full(xs.size, xs[1] - xs[0]); w[0] = w[-1] = w[0] / 2
    f = np.cos(np.outer(xs, a)) @ v; f /= f[np.argmax(abs(f))]
    ov = lambda g: float(abs((w * f * g).sum()) / np.sqrt((w * f * f).sum() * (w * g * g).sum()))
    r = {"kind": "float shape of certified witness", "file": fn, "L": L, "K": K, "RQ_certified": d["RQ"][:45],
         "sign_changes_on_[0,L]": int(((f[:-1] * f[1:]) < 0).sum()), "min/max": float(f.min()),
         "overlap_cosh(u/2)": ov(np.cosh(xs / 2)), "overlap_first_cosine": ov(np.cos(np.pi * xs / (2 * L)))}
    res.append(r); print(json.dumps(r))
json.dump(res, open("witness_shape_results.json", "w"), indent=1)
