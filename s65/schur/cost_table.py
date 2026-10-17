"""Grok S65 H65a: cost model, our one-stroke lift pipeline vs Mori/Liu bounded comparison with a given prime bound s.
Ours: T = 2 pi e^{A+beta}, K/parity ~ 0.815 L T (S58/S61/S63 runs).
Mori route: gamma = s + m, Omega solves log(Omega/2pi) - 1/Omega - gamma = C with C = 1.02*u0*E0(c), E0=(c+1)/(pi c^2), u0 = 428 (Mori VIII
values at log 15; assumed L-insensitive), c = Omega L; N_total ~ 1.67 c (Mori's empirical pass ratio), modes/parity = N/2. ESTIMATE only."""
import math, json
def omega(gam, L):
    O = 100.0
    for _ in range(200):
        c = O*L; C = 1.02*428*(c+1)/(math.pi*c*c); O = 2*math.pi*math.exp(gam + C + 1/O)
    return O, C
rows = []
for L, lab, A_ours, beta, s_quad, s_per, lam in [(1.4, '7/5', 5.0536, 0.26, 4.3813, 3.71150, 3.70611), (1.42, '71/50', 6.05, 0.20, None, 3.81418, None),
                                                  (math.log(15)/2, 'log15/2', 4.5928, 0.2, 4.0854, 3.49385, 3.48085)]:
    r = {"L": lab, "ours_T": 2*math.pi*math.exp(A_ours+beta), "ours_K": 0.815*L*2*math.pi*math.exp(A_ours+beta)}
    for nm, s in [("mori_w", 4.0915 if lab=='log15/2' else (4.6041 if lab=='7/5' else 4.5934)), ("quad", s_quad), ("perron", s_per)]:
        if s is None: continue
        for m in (0.5, 0.25):
            O, C = omega(s+m, L); r[f"{nm}_m{m}"] = {"s": s, "Omega": round(O), "C": round(C, 4), "modes_per_parity": round(1.67*O*L/2)}
    rows.append(r); print(json.dumps(r))
json.dump(rows, open('cost_table.json', 'w'), indent=1)
