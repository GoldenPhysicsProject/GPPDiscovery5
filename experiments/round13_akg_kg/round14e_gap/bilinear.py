"""14D redo (Golden Eye 3.2.1-2): bilinear radical test. For S = window-truncated range of E_F (own) and of a same-band
non-radical control E_F' with F' = F - 0.1*2^{-s} (i.e. a(2) -> a(2) - 0.1), report quadratic residual max|eig(S^T W S)|/|W|
and bilinear residual ||W S||_2/|W| = sup_{v in S, phi in V_M} |W(v,phi)|/(|v||phi||W|).  Separation = control/own."""
import json, math, numpy as np, gap
from gap import *
jobs = [('ZETA', md.zeta(), 4.0, 60, {}), ('ZETA', md.zeta(), 6.0, 60, {}), ('EPS', md.eps(), 6.0, 60, {}),
        ('FTH', md.fth(0.2), 4.0, 45, {}), ('DH', md.dh(), 4.0, 125, dict(widths=(1.5, 1.0), kf=1.0))]
out = []
for name, mo, L, F, kw in jobs:
    t, Fl = line(name); gap._line[name+'_P'] = np.vstack([t, Fl - 0.1*2.0**(-0.5-1j*t.real)])
    M = int(math.ceil(F*L/math.pi))+8; W = wmatrix(mo, L, M); nW = np.abs(np.linalg.eigvalsh(W)).max()
    for eps in [1e-8, 1e-6]:
        row = dict(model=name, L=L, M=M, eps=eps)
        for lab, en in [('own', name), ('control_a2', name+'_P')]:
            S, info = radical_space(en, mo, L, M, F, eps, **kw)
            row[lab+'_dimS'] = int(S.shape[1])
            if S.shape[1]:
                row[lab+'_quad'] = float(np.abs(np.linalg.eigvalsh(S.T @ W @ S)).max()/nW)
                row[lab+'_bil'] = float(np.linalg.norm(W @ S, 2)/nW)
        if row.get('own_bil') and row.get('control_a2_bil'): row['separation_bil'] = row['control_a2_bil']/row['own_bil']; row['separation_quad'] = row['control_a2_quad']/row['own_quad']
        out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open('bilinear.json', 'w'), indent=1)
