"""Round 6: Schur complement of the E_12 block against modes 13..12+K (float64, same spectral quadrature as p2_boundary_float.py).
NUMERICAL test of the leading-order prediction  M11 - M12 M22^{-1} M21 -> M11 - u u^T / g  (g = G(l,l))."""
import numpy as np, json
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
N = 12
res = []
for Kh in [20, 50, 100, 200, 400]:
    J = list(range(1, N+Kh+1))
    F = np.array([S(j) for j in J]); s = np.array([sig(j) for j in J])
    M = (F*meas) @ F.T - 2*bI2**2*np.outer(s, s)
    M11, M12, M22 = M[:N,:N], M[:N,N:], M[N:,N:]
    ev22 = np.linalg.eigvalsh(M22)
    Sc = M11 - M12 @ np.linalg.solve(M22, M12.T)
    wv = np.array([w_(j) for j in J[N:]])
    D = np.diag(M22) - wv**2*g
    res.append(dict(K=Kh, lam_min_M22=ev22[0], lam_min_full=np.linalg.eigvalsh(M)[0], lam_min_schur=np.linalg.eigvalsh(Sc)[0],
                    diag_excess_first_last=[D[0], D[-1]], diag_excess_times_om4_over_logom=[D[0]*om(N+1)**4/np.log(om(N+1)), D[-1]*om(N+Kh)**4/np.log(om(N+Kh))]))
    print(res[-1], flush=True)
uu = np.linalg.eigvalsh(M11 - np.outer(u,u)/g)[0]
print("leading-order limit lam_min(M11 - uu^T/g) =", uu)
json.dump(dict(sweep=res, leading_limit=uu), open('p2_schur_sweep_results.json','w'), indent=1, default=float)
