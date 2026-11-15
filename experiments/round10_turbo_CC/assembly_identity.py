"""Float NUMERICAL: Schur onto E_12 of A_X' - CC^T/lam_Z with CC^T from the exact identity Gram(k) - HH^T (full j>=73 tail)."""
import numpy as np, json
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
N = 12; JY = range(13, 73)
w_ = lambda j: np.sqrt(2/l)*(-1)**(j+1)/om(j)
Fx = [S(i) for i in range(1, N+1)] + [Fd] + [S(k) - w_(k)*Fd for k in JY]
sx = [sig(i) for i in range(1, N+1)] + [sd] + [sig(k) - w_(k)*sd for k in JY]
F = np.array(Fx); s = np.array(sx)
A = (F*meas) @ F.T - 2*bI2**2*np.outer(s, s)
CC = np.load('CCt_identity.npy'); Csp = np.load('Csp.npy')
def schurE(M):
    R = M[N:, N:]; Cc = M[:N, N:]
    return float(np.linalg.eigvalsh(M[:N, :N] - Cc @ np.linalg.solve(R, Cc.T))[0])
def ydPD(M):
    R = M[N:, N:]; d = 1/np.sqrt(np.diag(R)); return float(np.linalg.eigvalsh(R*np.outer(d, d))[0])
res = {}
for lam in [11.694, 5.24628, 1.0]:
    res[f'lam={lam}'] = dict(schur_identity=schurE(A - CC/lam), YdPD_identity=ydPD(A - CC/lam),
                             schur_spectral412=schurE(A - Csp @ Csp.T/lam))
res['schur_A_alone'] = schurE(A)
print(json.dumps(res, indent=1)); json.dump(res, open('assembly_identity_results.json', 'w'), indent=1)
