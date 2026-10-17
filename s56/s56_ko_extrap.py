"""Grok S56 H56b: fit delta_p(M) = d_inf + a*M^-alpha to the Ritz knockout thresholds (midpoints, half-width 0.00023). NUMERICAL."""
import json, glob, numpy as np
from scipy.optimize import curve_fit
rows = [json.loads(l) for f in glob.glob('h56b/ko_odd_p*.jsonl') for l in open(f)]
out = {}
for p in (3, 5):
    r = sorted([x for x in rows if x['p'] == p and x['EDGEK'] == '4' and 'delta_in' in x], key=lambda x: x['M'])
    M = np.array([x['M'] for x in r], float); d = np.array([sum(x['delta_in']) / 2 for x in r])
    f = lambda M, di, a, al: di + a * M ** -al
    res = {}
    for al0 in (1.0, 2.0):
        try:
            P, C = curve_fit(f, M, d, p0=[0.003, 0.5, al0], sigma=np.full(len(M), 0.00023), maxfev=20000)
            res[f'free_from{al0}'] = dict(d_inf=P[0], a=P[1], alpha=P[2], sd=list(np.sqrt(np.diag(C))))
        except Exception as e: res[f'free_from{al0}'] = str(e)
    for al in (0.5, 1.0, 2.0):
        A = np.vstack([np.ones_like(M), M ** -al]).T; c, *_ = np.linalg.lstsq(A, d, rcond=None)
        res[f'alpha={al}'] = dict(d_inf=c[0], a=c[1], rms=float(np.sqrt(np.mean((A @ c - d) ** 2))))
    out[p] = dict(M=M.tolist(), delta_mid=d.tolist(), fits=res)
    print(p, json.dumps(out[p], default=float))
json.dump(out, open('h56b/ko_extrap.json', 'w'), indent=1, default=float)
