"""S12: does local (per-prime) control + pole/arch coercivity prove Weil positivity past the prime-2 threshold? NUMERICAL (float, Ritz).
Compares lam0(h)=lam_min(pole+arch) on g in [-h,h] against the exact local bound sum_n 2 Lambda(n) n^-1/2 rho, rho<=1/2 (a>=h)."""
import sys, os, numpy as np
from scipy.linalg import eigh
sys.argv=[sys.argv[0],"14","0.3","0.3","1"]
sys.path.insert(0,"../experiments/weil_small_support")
import importlib.util
spec=importlib.util.spec_from_file_location("wr","../experiments/weil_small_support/weil_ritz.py")
src=open(spec.origin).read().split("cases = [")[0]
ns={"__file__":spec.origin}; exec(src,ns)
M=ns["M"]
for h in [0.35,0.40,0.45,0.50,0.54]:
    Q,G=ns["form"](h,"Z"); 
    ns_vm=ns["vm"]
    # rebuild prime part directly
    xq,wq=ns["xq0"]*h,ns["wq0"]*h; B=ns["basis"](xq,h)
    def Fmat(u):
        A=(B*wq)@ns["basis"](xq-u,h).T; return (A+A.T)/2
    Pr=sum(2*ns_vm(n)/np.sqrt(n)*Fmat(np.log(n)) for n in range(2,int(np.exp(2*h))+1) if ns_vm(n))
    Q0=Q+Pr
    for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
        w,v=eigh(Q0[np.ix_(I,I)],G[np.ix_(I,I)]); l0=w[0]; x=v[:,0]
        wf,_=eigh(Q[np.ix_(I,I)],G[np.ix_(I,I)])
        bound=sum(2*ns_vm(n)/np.sqrt(n)*0.5 for n in range(2,int(np.exp(2*h))+1) if ns_vm(n))
        pr_on=x@Pr[np.ix_(I,I)]@x/(x@G[np.ix_(I,I)]@x)
        print(f"h={h:.2f} {nm}: lam0(pole+arch)={l0:.3e}  localbound={bound:.3f}  prime-term on bottom vec={pr_on:+.3e}  lam_full={wf[0]:+.3e}")
