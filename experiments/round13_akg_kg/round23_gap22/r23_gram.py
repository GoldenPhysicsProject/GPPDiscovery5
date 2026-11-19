"""Round 23, route (a): Gram/Toeplitz test of the gamma = 0 slice.
G(a,b) = D(a) + D(b) - D(|a-b|), D(0) = 0.  G PSD on every finite set  <=>  D is of negative-definite type
<=>  (Schoenberg/Levy-Khintchine, with D(0)=0, D even, the Laplace identity)  D(L) = sum_rho (1-cos gamma L)/gamma^2 with all gamma real
<=> RH.  (PROVED equivalence modulo the standard Schoenberg theorem; the finite-grid test itself is NUMERICAL.)
Off-line zeros contribute terms ~ cosh((beta-1/2)L)-type growth which are not of negative-definite type; the test must fail for fakes.
D is built only from the integers (closed form / explicit formula); no zeros are inputs."""
import json, math, numpy as np
from r22_fejer import rho
from r22_models import lam_table
from r23_slice import D_zeta
out = {}
Lmax = math.log(1e5); N = 10**5 + 1
for name in ('ZETA', 'FTH', 'EPS', 'DH'):
    lam = lam_table(name, N)
    out[name] = {}
    for h in (0.25, 0.1, 0.05):
        grid = np.arange(1, int(Lmax/h) + 1)*h
        allL = np.arange(0, int(Lmax/h) + 1)*h
        Dv = np.zeros(len(allL))
        for i, L in enumerate(allL[1:], 1):
            Dv[i] = (np.pi*L/2*rho(name, math.exp(L), np.array([0.0]), lam)[0])
        n = len(grid); I = np.arange(1, n + 1)
        G = Dv[I][:, None] + Dv[I][None, :] - Dv[np.abs(I[:, None] - I[None, :])]
        w = np.linalg.eigvalsh(G)
        # smallest leading principal block that is already indefinite
        first = None
        for m in range(2, n + 1):
            if np.linalg.eigvalsh(G[:m, :m])[0] < -1e-10*max(1, abs(G[:m, :m]).max()):
                first = m; break
        out[name][str(h)] = dict(n=n, min_eig=float(w[0]), max_eig=float(w[-1]),
                                 first_indefinite_block=first, L_at_first=float(grid[first - 1]) if first else None)
        print(name, h, out[name][str(h)], flush=True)
json.dump(out, open('r23_gram.json', 'w'), indent=1)
