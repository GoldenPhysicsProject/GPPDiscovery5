"""Grok S33: floor lambda_max(P_L) (s28_pmax Galerkin, NUMERICAL) at L=1.28, 1.3. RH OPEN."""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's28'))
from s28_pmax import Pmat, mangoldt_terms
res = {}
for L in [1.25, 1.28, 1.3]:
    vals = [float(np.linalg.eigvalsh(Pmat(L, N))[-1]) for N in (40, 80, 120)]
    res[str(L)] = dict(lam_max=vals[-1], conv=vals, A_L_raw=sum(2*lp/math.sqrt(n) for n, lp in mangoldt_terms(L)), status='NUMERICAL')
    print(L, vals, flush=True)
json.dump(res, open('h33_pmax.json', 'w'), indent=1)
