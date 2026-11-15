"""Round 7 (b) float study (NUMERICAL): three-block route X=E_12, delta, Psi = span{psi_j = phi_j - w_j delta, j >= 13}.
(T1) coercivity c of B on Psi in D-units, D_j = 3 w(omega_j)/omega_j^4;  (T2) couplings v_j = (B(phi_i,psi_j))_i, B(delta,psi_j)) and omega^3 |v_j|;
(T3) Schur onto X of [[M11,u],[u^T,g]] - sum_j v_j v_j^T/(c D_j)."""
import numpy as np, json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
wsym = lambda x: special.digamma(0.25+0.5j*x).real - np.log(np.pi)
N = 12; Kh = 400; J = list(range(N+1, N+Kh+1))
Fx = np.array([S(i) for i in range(1, N+1)] + [Fd]); sx = np.array([sig(i) for i in range(1, N+1)] + [sd])
Gx = (Fx*meas) @ Fx.T - 2*bI2**2*np.outer(sx, sx)
Fp = np.array([S(j) - w_(j)*Fd for j in J]); sp = np.array([sig(j) - w_(j)*sd for j in J])
Gp = (Fp*meas) @ Fp.T - 2*bI2**2*np.outer(sp, sp)
V = (Fx*meas) @ Fp.T - 2*bI2**2*np.outer(sx, sp)          # 13 x K couplings
Dd = np.array([3*wsym(om(j))/om(j)**4 for j in J]); Ds = 1/np.sqrt(Dd)
c = np.linalg.eigvalsh(Gp*np.outer(Ds, Ds))[0]
oms = np.array([om(j) for j in J])
env = np.abs(V)*oms**3
def schurX(A): return np.linalg.eigvalsh(A[:N,:N] - np.outer(A[:N,N], A[N,:N])/A[N,N])[0]
exact = np.linalg.eigvalsh((Gx - V @ np.linalg.solve(Gp, V.T))[:N,:N] - 0)[0]  # not Schur onto X yet
Mfull = np.block([[Gx, V], [V.T, Gp]])
lam_exact = np.linalg.eigvalsh(Mfull[:N,:N] - Mfull[:N,N:] @ np.linalg.solve(Mfull[N:,N:], Mfull[N:,:N]))[0]
Lc = (V/(c*Dd)) @ V.T
lam_cD = schurX(Gx - Lc)
lam_lead = schurX(Gx)
out = dict(K=Kh, c_all_j_ge_13=float(c), lam_exact=float(lam_exact), lam_leading=float(lam_lead), lam_with_cD_bound=float(lam_cD),
           env_omega3_absv_max_per_row=[float(x) for x in env.max(axis=1)], env_at_j=[[int(J[k]), float(env[:,k].max())] for k in [0, 10, 59, 60, 100, 200, 399]])
# tail-only loss beyond 72 with cD bound
k72 = 60; L72 = (V[:, k72:]/(c*Dd[k72:])) @ V[:, k72:].T
out['tail_beyond72_loss_norm_cD'] = float(np.linalg.norm(L72, 2))
print(json.dumps(out, indent=1)); json.dump(out, open('r7_tailbound_float_results.json','w'), indent=1)
