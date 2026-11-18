"""Quadrature convergence of the one-shift block at zeta's tightest shift (s = 5.78) and FTAU's (7.34): npan 24 vs 64."""
import numpy as np, json, hasse_core as hc
out = []
for name, s in [('ZETA', 5.78), ('FTAU', 7.34), ('ZETA', 4.0)]:
    mo = hc.MODELS[name](); l = 0.6; M = hc.Mof(l, 45.0)
    W = hc.Bs(mo, l, M, 0.0)['B']; W = 0.5*(W+W.T)
    r = {}
    for npan in (24, 64):
        B = hc.Bs(mo, l, M, s, npan=npan)['B']; B = 0.5*(B+B.T)
        r[npan] = (float(np.linalg.eigvalsh(W+B)[0]), float(np.linalg.eigvalsh(W-B)[0]))
    out.append(dict(model=name, s=s, npan24=r[24], npan64=r[64], diff=float(max(abs(a-b) for a, b in zip(r[24], r[64])))))
    print(out[-1])
json.dump(out, open('check_quad.json', 'w'), indent=1)
