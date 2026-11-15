"""Round 7 (b): is the tail coupling v_j rank-one to leading order, v_j ~ (w_j/omega_j^2) q, q = second-boundary vector (NUMERICAL)?
Residual envelope omega^4 |v_j - (w_j/omega_j^2) q| tests the O(omega^-4) remainder (interior term after Riemann-Lebesgue)."""
import numpy as np, json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
wsym = lambda x: special.digamma(0.25+0.5j*x).real - np.log(np.pi)
N = 12; J = list(range(73, 473))
Fx = np.array([S(i) for i in range(1, N+1)] + [Fd]); sx = np.array([sig(i) for i in range(1, N+1)] + [sd])
Fp = np.array([S(j) - w_(j)*Fd for j in J]); sp = np.array([sig(j) - w_(j)*sd for j in J])
V = (Fx*meas) @ Fp.T - 2*bI2**2*np.outer(sx, sp)
oms = np.array([om(j) for j in J]); wv = np.array([w_(j) for j in J])
scal = wv/oms**2
q = (V @ scal)/(scal @ scal)            # least-squares rank-one direction
R = V - np.outer(q, scal)
Dd = 3*wsym(oms)/oms**4
out = dict(q=q.tolist(), rel_resid=float(np.linalg.norm(R)/np.linalg.norm(V)),
           resid_env_omega4=[[int(J[k]), float(np.abs(R[:,k]).max()*oms[k]**4)] for k in [0, 50, 100, 200, 399]],
           kappa_rank1_beyond72_c1=float(np.sum(scal**2/Dd)), resid_loss_beyond72_c1=float(np.linalg.norm((R/Dd) @ R.T, 2)))
print(json.dumps(out, indent=1)); json.dump(out, open('r7_rank1tail_float_results.json','w'), indent=1)
