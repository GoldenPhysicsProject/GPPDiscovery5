"""S82-OS H82-OS-b probe (NUMERICAL, float64). Pole-free Weil form Q0 = Q - pole on the even and odd sectors of [-L,L].
Questions: (i) how many negative eigenvalues has Q0 (even/odd)? (ii) is the ground state of Q0_even a positive function
(Krein-Rutman / Perron, as Lambda>=0 + Markov archimedean resolvent would suggest)? (iii) overlap of that ground state with
the pole vector cosh(u/2). Also DH (no pole; signed b(n)). RH OPEN."""
import json, math
import numpy as np
from numpy.polynomial import legendre as leg
from os_hankel import build
def fvals(L, N, v, xs):
    t = 2 * xs / L - 1
    return sum(v[k] * math.sqrt((2 * k + 1) / L) * leg.legval(t, [0] * k + [1]) for k in range(N))
res = []
for FORM in ("ZETA", "DH"):
    for L in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
        N = 28
        T, B = build(FORM, L, N); c = build.comp
        ape = c['ap'] + c['am']; apo = c['ap'] - c['am']
        Qe = T + B; Qo = T - B
        Q0e = Qe - np.outer(ape, ape) if FORM == 'ZETA' else Qe
        Q0o = Qo + np.outer(apo, apo) if FORM == 'ZETA' else Qo
        ee, Ve = np.linalg.eigh(Q0e); eo = np.linalg.eigvalsh(Q0o)
        xs = np.linspace(0, L, 801); g = fvals(L, N, Ve[:, 0], xs); g *= np.sign(g[np.argmax(abs(g))])
        r = {"kind": "float", "form": FORM, "L": L, "N": N,
             "Q0_even_eigs_low3": [float(x) for x in ee[:3]], "Q0_even_neg": int((ee < -1e-12).sum()),
             "Q0_odd_eigs_low2": [float(x) for x in eo[:2]], "Q0_odd_neg": int((eo < -1e-12).sum()),
             "ground_state_min/max": float(g.min() / g.max()),
             "overlap(ground,cosh)": float(abs(Ve[:, 0] @ ape) / np.linalg.norm(ape)),
             "full_even_min": float(np.linalg.eigvalsh(Qe)[0]), "full_odd_min": float(np.linalg.eigvalsh(Qo)[0])}
        res.append(r); print(json.dumps(r), flush=True)
json.dump(res, open("perron_pole_results.json", "w"), indent=1)
