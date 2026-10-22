"""Tail diagnostics for the zeta Weil-form certificate at a window [-L,L] (L = HALF-width), even sector.
Input: a cosritz_dump W (unnormalised cos basis, Gram = L*I) of size Kbig, and the normalised soft vector v of the
leading K block (vecs.py). Outputs (all NUMERICAL, float64 except where stated):
  mu(K0)   = lambda_min(W[K0:,K0:])/L        (Rayleigh floor on the high block, team units Q/||f||^2)
  eta(K0)  = ||W[K0:,:K0]||_2 / L             (coupling norm)
  eta^2/mu vs lambda_K0                        -> the Schur lift lemma (PROVED) is usable iff lambda_K0 >= eta^2/mu
  d2       = r^T W_HH^{-1} r / L, r = W_HL v  (mpmath, dps from dump; second-order shift of the soft eigenvalue
              when the modes K..Kbig-1 are added; NOT a bound)
usage: schur_tail.py Wdump.json vec.json "K0,K0,.." out.json
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, numpy as np, mpmath as mp
D=json.load(open(sys.argv[1])); V=json.load(open(sys.argv[2])); K0s=[int(x) for x in sys.argv[3].split(',')]
W=D["W"]; N=len(W); Ls=D["L"]; p,q=(Ls.split('/')+['1'])[:2]; Lf=float(p)/float(q)
Wf=np.array([[float(x[0]) for x in r] for r in W])
res={"dump":sys.argv[1],"L":Ls,"Kbig":N,"rows":[]}
for K0 in K0s:
    H=Wf[K0:,K0:]; C=Wf[K0:,:K0]
    mu=np.linalg.eigvalsh(H)[0]/Lf; eta=np.linalg.norm(C,2)/Lf
    res["rows"].append({"K0":K0,"mu":mu,"eta":eta,"eta2_over_mu":eta**2/mu if mu>0 else None,
        "lam_min_LL_float":np.linalg.eigvalsh(Wf[:K0,:K0])[0]/Lf})
    print(res["rows"][-1],flush=True)
# second-order shift with the converged soft vector of the leading K block
mp.mp.dps=200; K=len(V["coef"]); v=[mp.mpf(c) for c in V["coef"]]
WHL=[[mp.mpf(W[i][j][0]) for j in range(K)] for i in range(K,N)]
r=[mp.fsum(WHL[i][j]*v[j] for j in range(K)) for i in range(N-K)]
WHH=mp.matrix([[mp.mpf(W[i][j][0]) for j in range(K,N)] for i in range(K,N)])
x=mp.lu_solve(WHH,mp.matrix(r))
d2=mp.fsum(r[i]*x[i] for i in range(N-K))/Lf
res.update({"K_soft":K,"Q_soft_K":V["Q_over_norm2"],"norm_r_over_L":mp.nstr(mp.sqrt(mp.fsum(t*t for t in r))/Lf,12),
            "second_order_shift":mp.nstr(d2,12),"predicted_Q_at_Kbig":mp.nstr(mp.mpf(V["Q_over_norm2"])-d2,12)})
print({k:res[k] for k in ("K_soft","Q_soft_K","norm_r_over_L","second_order_shift","predicted_Q_at_Kbig")})
json.dump(res,open(sys.argv[4],"w"),indent=1)
