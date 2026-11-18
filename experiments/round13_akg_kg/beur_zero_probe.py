"""Check the Re s > 0 continuation of log zeta_B'' used in debranges.py BEUR: |zeta_B| must vanish at the quartet."""
import numpy as np, math, mpmath as mp, json
import harness as h_
Zm = h_.zeta_model(); B = h_.beurling_model(x0=55.0); m = Zm['n'] < 55
n = np.concatenate([Zm['n'][m], B['n']]); c = np.concatenate([Zm['c'][m], B['c']]); wts = c/np.log(n); logn = np.log(n); ell = math.log(n.max())
def logzB(s): return np.sum(wts*np.exp(-s*logn)) + complex(mp.e1((s-1)*ell)) - sum(complex(mp.e1((s-r)*ell)) for r in B['quartet'])
out = {}
for sig in (0.3, 0.4, 0.5, 0.6, 0.7):
    ts = np.arange(26.0, 30.0, 0.02); v = [logzB(complex(sig, t)).real for t in ts]; i = int(np.argmin(v))
    out[str(sig)] = dict(min_log_abs_zeta=float(v[i]), t=float(ts[i]), median=float(np.median(v)))
print(json.dumps(out)); json.dump(out, open('beur_zero_probe.json', 'w'), indent=1)
