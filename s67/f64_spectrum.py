"""Grok S67: float64 spectrum diagnostic of S66 dumps (NUMERICAL).  lam_near0 in S66 returns the eigenvalue
NEAREST ZERO (inverse power iteration), not lambda_min.  Here: eigvalsh in float64 of A+J and T; any eigenvalue
<< -1e-12*||M|| is a genuine negative eigenvalue (float64 resolves it); values near 0 are below float resolution."""
import sys, pickle, json, numpy as np
for fn in sys.argv[1:]:
    D = pickle.load(open(fn, 'rb')); out = {'dump': fn.split('/')[-1]}
    for k in ('T', 'A', 'J'):
        D[k] = np.array([[float(x) for x in r] for r in D[k]])
    for name, M in (('T', D['T']), ('A+J', D['A'] + D['J'])):
        M = (M + M.T) / 2; ev = np.linalg.eigvalsh(M); nrm = abs(ev).max()
        neg = ev[ev < -1e-12 * nrm]
        out[name] = dict(norm=float(nrm), lowest5=[float(x) for x in ev[:5]], n_clearly_negative=int(len(neg)))
    print(json.dumps(out), flush=True); del D
