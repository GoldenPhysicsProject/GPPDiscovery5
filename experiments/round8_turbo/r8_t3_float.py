"""Round 8: (T3) float assembly with tail in closed-form shape (NUMERICAL).
Psi = span{psi_j, j>=13} as ONE block with B|_Psi >= c D (c = 0.6892, float, (T1) GAP).  Exact v_j for 13..72 (accurate run, T=2e5);
tail j>=73: v_j = a_j q + rho_j, a_j = w_j/omega_j^2, |rho_j| <= R(omega)/omega^4 with R = 2.1 log(omega) (envelope from sweep).
vv^T <= (1+eta) a^2 qq^T + (1+1/eta) |rho|^2 I."""
import numpy as np, json
exec(open('r8_t2.py').read().split("out = dict")[0].replace("JJ = [73, 123, 173, 273, 472]", "JJ = list(range(13, 73))"))
q, r4, v = run(2e5, chunk=20000)
src = open('../round6_turbo/p2_boundary_float.py').read().split("out = dict")[0]
ns = {}; exec(src, ns)
Mx = np.array([[ns['B'](ns['S'](i),ns['sig'](i),ns['S'](k),ns['sig'](k)) for k in range(1,13)] for i in range(1,13)])
u = np.array([ns['B'](ns['S'](i),ns['sig'](i),ns['Fd'],ns['sd']) for i in range(1,13)]); g = ns['g']
Gx = np.block([[Mx, u[:,None]],[u[None,:], np.array([[g]])]])
wsym = lambda x: special.digamma(0.25+0.5j*x).real - np.log(np.pi)
c = 0.6892
D = lambda j: 3*wsym(om(j))/om(j)**4
def schurX(A): return np.linalg.eigvalsh(A[:N,:N] - np.outer(A[:N,N], A[N,:N])/A[N,N])[0]
Lfin = sum(np.outer(v[:,k], v[:,k])/(c*D(j)) for k,j in enumerate(JJ))
jt = np.arange(73, 2_000_001, dtype=float); omt = jt*np.pi/l; wt_ = special.digamma(0.25+0.5j*omt).real - np.log(np.pi)
Dt = 3*wt_/omt**4; a2 = (2/l)/omt**6
kappa = float(np.sum(a2/Dt)); kappa += (2/l)/(3*wt_[-1])*(l/np.pi)**2/jt[-1]      # integral tail beyond 2e6 (w slowly varying)
R = 2.1*np.log(omt); epsr = float(np.sum((R/omt**4)**2/Dt))
res = {}
for eta in [0.05, 0.1, 0.2, 0.5]:
    A = Gx - Lfin - (1+eta)*kappa*np.outer(q,q)/c - (1+1/eta)*epsr/c*np.eye(N+1)
    res[str(eta)] = float(schurX(A))
lam_lead = float(schurX(Gx)); lam_fin = float(schurX(Gx - Lfin))
out = dict(c=c, kappa_73_inf=kappa, eps_rho_73_inf=epsr, lam_leading=lam_lead, lam_after_13_72=lam_fin, lam_T3_by_eta=res,
           q=q.tolist(), g=float(g))
print(json.dumps(out, indent=1)); json.dump(out, open('r8_t3_float_results.json','w'), indent=1)
# Row-weighted remainder bound: |rho_jx| <= R_x(omega)/omega^4, rho rho^T <= (sum_x R_x) diag(R_x) / omega^8  (Cauchy-Schwarz).
# R_delta = 2.1 log(omega); R_phi = 0.045 * omega_73/omega (envelopes from r8_sweep, NUMERICAL).
Rd = 2.1*np.log(omt); Rp = 0.045*om(73)/omt
Wsum = Rd + N*Rp
dg_phi = float(np.sum(Wsum*Rp/omt**8/Dt)); dg_del = float(np.sum(Wsum*Rd/omt**8/Dt))
res2 = {}
for eta in [0.02, 0.05, 0.1, 0.2, 0.5]:
    A = Gx - Lfin - (1+eta)*kappa*np.outer(q,q)/c - (1+1/eta)/c*np.diag([dg_phi]*N + [dg_del])
    res2[str(eta)] = float(schurX(A))
out.update(rowweighted_diag_phi=dg_phi, rowweighted_diag_delta=dg_del, lam_T3_rowweighted_by_eta=res2)
print(json.dumps({k: out[k] for k in ['rowweighted_diag_phi','rowweighted_diag_delta','lam_T3_rowweighted_by_eta']}, indent=1))
json.dump(out, open('r8_t3_float_results.json','w'), indent=1)
