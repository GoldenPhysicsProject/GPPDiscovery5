"""Round 14E driver. PART=pair: Sylvester / inertia test on the fakes (quotient by OWN window-truncated range E_F).
PART=zeta: zeta lambda_min on S_eps^perp vs L, eps (concentration cut) and fmax (basis size)."""
import os, sys, json, math, time, numpy as np
from gap import *
PART = sys.argv[1]; out = []
def M_of(F, L): return int(math.ceil(F*L/math.pi)) + 8
if PART == 'pair':
    jobs = [('EPS', md.eps(), 'EPS', [4.0, 6.0], 60, {}), ('FTH', md.fth(0.2), 'FTH', [2.0, 4.0], 45, {}),
            ('DH', md.dh(), 'DH', [4.0, 6.0], 125, dict(widths=(1.5, 1.0), kf=1.0)),
            ('BEUR2', md.beur2(), 'ZETA', [4.5, 6.0], 60, {}),     # E_B not admissible (see md); quotient by zeta's range
            ('ZETA', md.zeta(), 'ZETA', [4.0, 6.0], 60, {})]
    for lab, mo, ename, Ls, F, kw in jobs:
        for L in Ls:
            M = M_of(F, L); W = wmatrix(mo, L, M)
            for eps in [1e-10, 1e-8, 1e-6]:
                t0 = time.time(); r = analyse(ename, mo, L, M, F, eps, W=W, **kw); r.update(model=lab, E_from=ename, sec=round(time.time()-t0, 1))
                out.append(r); print(json.dumps({k: r[k] for k in ['model','E_from','L','M','eps','dimS','n_neg_W','n_neg_perp','lam_min_W','lam_min_perp','bilinear_resid_rel','radical_resid_rel','sec'] if k in r}), flush=True)
            json.dump(out, open('pair.json', 'w'), indent=1)
if PART == 'zeta':
    mo = md.zeta()
    for F in [45, 90]:
        for L in [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
            M = M_of(F, L); W = wmatrix(mo, L, M)
            for eps in [1e-10, 1e-8, 1e-6]:
                t0 = time.time(); r = analyse('ZETA', mo, L, M, F, eps, W=W, widths=(1.5, 1.0), kf=1.1); r.update(sec=round(time.time()-t0, 1))
                out.append(r); print(json.dumps({k: r[k] for k in ['F','L','M','eps','dimS','n_tiny_W','dim_perp','lam_min_W','lam_min_perp','lam2_perp','bilinear_resid_rel','sec'] if k in r}), flush=True)
            json.dump(out, open('zscan.json', 'w'), indent=1)
