"""Float NUMERICAL: three-block assembly X' = E_12 + delta + Y(psi_13..72) exact, Z = H_73 (psi_j, j >= 73) replaced by its uniform
lower bound lam_Z * ||g||^2 (g-units: psi_j <-> phi_j / omega_j^2). Compare Schur onto E_12 with exact."""
import numpy as np, json
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
N = 12; J = list(range(13, 413))
Fx = [S(i) for i in range(1, N+1)] + [Fd]; sx = [sig(i) for i in range(1, N+1)] + [sd]
Fp = [S(j) - w_(j)*Fd for j in J]; sp = [sig(j) - w_(j)*sd for j in J]
F = np.array(Fx + Fp); s = np.array(sx + sp)
M = (F*meas) @ F.T - 2*bI2**2*np.outer(s, s)
n0 = N + 1; ny = 60            # X' = first 13 + 60
oms = np.array([om(j) for j in J])
def schurE(A):
    k = A.shape[0]; R = A[N:, N:]; Cc = A[:N, N:]
    return np.linalg.eigvalsh(A[:N,:N] - Cc @ np.linalg.solve(R, Cc.T))[0]
exact = schurE(M)
iX = n0 + ny
MZ = M[iX:, iX:]; scl = oms[ny:]**2      # g-normalised Z basis: phi_j = omega_j^2 psi_j
MZg = MZ*np.outer(scl, scl)
lamZ_float = np.linalg.eigvalsh(MZg)[0]
res = dict(schur_exact_K400=exact, lamZ_float_g_units_73_412=lamZ_float, diagZ_g=[float(MZg[0,0]), float(MZg[-1,-1])])
C = M[:iX, iX:]*scl[None, :]
for lamZ in [lamZ_float, 5.0, 3.0, 1.0]:
    A = M[:iX, :iX] - C @ C.T/lamZ
    res[f'schur_with_lamZ={lamZ:.3f}'] = schurE(A)
    res[f'YdeltaPD_with_lamZ={lamZ:.3f}'] = float(np.linalg.eigvalsh(A[N:,N:]*np.outer(1/np.sqrt(np.diag(A[N:,N:])),1/np.sqrt(np.diag(A[N:,N:]))))[0])
# size of Y-Z coupling in g-units
Cg = C[n0:, :]*np.array([om(j)**2 for j in range(13,73)])[:,None]
res['max_row_norm_YZ_g'] = float(np.sqrt((Cg**2).sum(1)).max()); res['Y_diag_g_min'] = float((np.diag(M[n0:iX,n0:iX])*np.array([om(j)**4 for j in range(13,73)])).min())
print(json.dumps(res, indent=1, default=float)); json.dump(res, open('assembly_float_results.json','w'), indent=1, default=float)
