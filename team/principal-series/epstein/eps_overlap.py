"""Nodes of the lowest even Galerkin vector of the Epstein (x^2+5y^2) window form and its overlaps with the EVEN parts
Im f_Q = -sin(g0 x) sinh(b0 x) of Queen-style vectors at Celestial's 13 CERTIFIED off-line zeros (Re f_Q is odd: overlap 0 exactly).
L = half-width. NUMERICAL."""
import json, glob, sys, numpy as np
Z = [((d['box'][0]+d['box'][1])/2 - 0.5, (d['box'][2]+d['box'][3])/2) for d in json.load(open('certified_zero_boxes_celestial_copy.json')) if d['Q'] == [1,0,5] and d['count'] == 1]
for fn in sorted(glob.glob('lkvec_EPS_*.json'), key=lambda f: (float(f.split('_L')[1].split('_K')[0]), int(f.split('_K')[1][:-5]))):
    d = json.load(open(fn)); L = float(d['L']); v = np.array([float(t) for t in d['vec']]); K = len(v)
    a = (2*np.arange(K)+1)*np.pi/(2*L); x = np.linspace(-L, L, 200001); f = np.cos(np.outer(x, a)) @ v
    core = np.abs(x) < L - 1e-4; s = np.sign(f[core]); nodes = int(np.sum(s[1:] != s[:-1]))
    ip = lambda p, q: np.trapezoid(p*q, x); nf = np.sqrt(ip(f, f))
    ov = []
    for b0, g0 in Z:
        gE = -np.sin(g0*x)*np.sinh(b0*x); ov.append((round(g0, 2), float('%.3g' % (ip(f, gE)/nf/np.sqrt(ip(gE, gE))))))
    best = max(ov, key=lambda t: abs(t[1]))
    pk = int(np.argmax(np.abs(v))); mass = v**2/np.sum(v**2)
    print(json.dumps(dict(L=L, K=K, ev=d['ev'][:3], nodes=nodes, peak_freq=round(float(a[pk]), 2), mass_top_modes=[round(float(t), 3) for t in mass[:4]],
          mass_last10=float('%.2g' % mass[-10:].sum()), best_overlap=best, overlaps=ov)))
