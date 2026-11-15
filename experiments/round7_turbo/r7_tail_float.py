"""Round 7 (b) structural float test (NUMERICAL): is the boundary-corrected tail psi_j = phi_j - w_j delta_l coercive at the window rate?
psi_j transform: c2 S_j(t) - w_j sin(t l);  pole coefficient sig_j - w_j sinh(l/2).  D_j = 3 w(omega_j)/omega_j^4 (window-symbol scale)."""
import numpy as np, json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
wsym = lambda x: special.digamma(0.25+0.5j*x).real - np.log(np.pi)
res = {}
for (j0, Kh) in [(13, 60), (73, 100), (73, 200), (73, 400)]:
    J = list(range(j0, j0+Kh))
    Fp = np.array([S(j) - w_(j)*Fd for j in J]); sp = np.array([sig(j) - w_(j)*sd for j in J])
    Gp = (Fp*meas) @ Fp.T - 2*bI2**2*np.outer(sp, sp)
    Fh = np.array([S(j) for j in J]); sh = np.array([sig(j) for j in J])
    Gh = (Fh*meas) @ Fh.T - 2*bI2**2*np.outer(sh, sh)
    wv = np.array([w_(j) for j in J]); H = Gh - g*np.outer(wv, wv)
    Dd = np.array([3*wsym(om(j))/om(j)**4 for j in J]); Ds = 1/np.sqrt(Dd)
    r = dict(lam_psi_Dnorm=float(np.linalg.eigvalsh(Gp*np.outer(Ds,Ds))[0]), lam_H_Dnorm=float(np.linalg.eigvalsh(H*np.outer(Ds,Ds))[0]),
             diag_psi_over_D=[float(Gp[0,0]/Dd[0]), float(Gp[-1,-1]/Dd[-1])], lam_psi=float(np.linalg.eigvalsh(Gp)[0]))
    res[f"{j0}+{Kh}"] = r; print(j0, Kh, r, flush=True)
json.dump(res, open('r7_tail_float_results.json','w'), indent=1)
