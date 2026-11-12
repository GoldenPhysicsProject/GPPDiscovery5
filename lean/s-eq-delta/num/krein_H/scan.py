# Cheapest falsifier for H (triangle psi): min of Q on ker l_c ∩ span{tau_{x_i} psi}, Q from the ZERO-FREE side.
import numpy as np, json
import importlib.util, sys
sys.argv=["x","1"]
spec=importlib.util.spec_from_file_location("h","hef_tri.py")
src=open("hef_tri.py").read().split("mp.mp.dps=20")[0]
ns={}; exec(src,ns)
h=0.25; N=24
Kp=np.array([ns["K_prime"](k*h) for k in range(N+1)])
np.save("Kp.npy",Kp)
res=[]
for n in range(2,N+2):
    x=np.arange(n)*h
    M=np.array([[Kp[abs(i-j)] for j in range(n)] for i in range(n)])
    a=np.cosh(x/2); a/=np.linalg.norm(a)
    P=np.eye(n)-np.outer(a,a)
    # basis of a-perp
    U,_,_=np.linalg.svd(P); B=U[:,:n-1]
    ev=np.linalg.eigvalsh(B.T@M@B)
    evfull=np.linalg.eigvalsh(M)
    res.append({"n":n,"L_window":1+ (n-1)*h,"min_eig_kerlc":float(ev[0]),"min_eig_full":float(evfull[0]),"max_eig":float(evfull[-1])})
    print(res[-1],flush=True)
json.dump({"h":h,"Kp":Kp.tolist(),"rows":res},open("scan.json","w"),indent=1)
