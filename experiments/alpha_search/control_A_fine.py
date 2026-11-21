"""Pin down the chance rate for grammar A (Heyrovska/Sherbon/Pellis shape) with 200000 decoys."""
import numpy as np, json
T = 137.035999177; phi = (1 + 5 ** .5) / 2
t2 = np.array([s * a * phi ** -n for s in (1, -1) for a in range(1, 10) for n in range(1, 10)])
t3 = np.array([s * b * phi ** -m for s in (1, -1) for b in range(1, 10) for m in range(1, 10)]
              + [s * (b * phi) ** -m for s in (1, -1) for b in range(2, 10) for m in range(1, 10)])
v = np.unique(np.round((360 * phi ** -2 + t2[:, None] + t3[None, :]).ravel(), 11))
rng = np.random.default_rng(7)
out = {}
for W in (0.01, 0.5):
    d = T + rng.uniform(-W, W, 200000); d = d[np.abs(d - T) > 1e-5]
    for eps in (2.1e-8, 1.6e-7):
        h = np.searchsorted(v, d + eps) - np.searchsorted(v, d - eps)
        out[f"W={W} eps={eps}"] = {"decoys": int(d.size), "decoys_with_hit": int((h > 0).sum()), "rate": float((h > 0).mean())}
print(json.dumps(out, indent=1)); json.dump(out, open("control_A_fine.json", "w"), indent=1)
