"""Round 6: probe of the named remainder quantities (float64; NUMERICAL, not certified).
H = M22 - g w w^T (high block minus boundary rank-one), E = M12 - u w^T (cross block minus boundary term), D_j = 3 w(om_j)/om_j^4
(window-symbol diagonal: (2/l)(1/pi) * avg Psi^2 = 6/t^4 * w(t) * ∫S_j^2 = pi l/4)."""
import numpy as np, json
from scipy import special
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
N = 12; res = []
for Kh in [50, 100, 200, 400]:
    J = list(range(1, N+Kh+1))
    F = np.array([S(j) for j in J]); s = np.array([sig(j) for j in J])
    M = (F*meas) @ F.T - 2*bI2**2*np.outer(s, s)
    M11, M12, M22 = M[:N,:N], M[:N,N:], M[N:,N:]
    oms = np.array([om(j) for j in J[N:]]); wv = np.array([w_(j) for j in J[N:]])
    D = 3*(special.digamma(0.25+0.5j*oms).real - np.log(np.pi))/oms**4
    H = M22 - g*np.outer(wv, wv); E = M12 - np.outer(u, wv)
    Dm = 1/np.sqrt(D)
    kap = np.linalg.eigvalsh(Dm[:,None]*H*Dm[None,:])
    EH = E @ np.linalg.solve(H, E.T)
    res.append(dict(K=Kh, kappa_min=kap[0], kappa_max=kap[-1], normsq_E_Dinvhalf=np.linalg.norm(E*Dm[None,:], 2)**2,
                    lam_max_E_Hinv_Et=np.linalg.eigvalsh(EH)[-1],
                    lam_min_schur=np.linalg.eigvalsh(M11 - M12 @ np.linalg.solve(M22, M12.T))[0]))
    print(res[-1], flush=True)
json.dump(res, open('p2_remainder_probe_results.json','w'), indent=1, default=float)
