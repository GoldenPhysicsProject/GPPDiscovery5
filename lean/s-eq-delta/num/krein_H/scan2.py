import numpy as np, json
Kp=np.load("Kp.npy"); h=0.25; N=24
def run(n,shift):
    x=(np.arange(n)-(n-1)/2)*h+shift
    M=np.array([[Kp[abs(i-j)] for j in range(n)] for i in range(n)])
    a=np.cosh(x/2); a/=np.linalg.norm(a)
    U,_,_=np.linalg.svd(np.eye(n)-np.outer(a,a)); B=U[:,:n-1]
    w,V=np.linalg.eigh(B.T@M@B)
    return float(w[0]), (B@V[:,0]).tolist(), x.tolist()
rows=[]
for n in range(2,N+2):
    best=None
    for s in np.arange(0,3.01,0.25):
        m,v,x=run(n,s)
        L=1+max(abs(x[0]),abs(x[-1]))   # window [-L,L] containing all supports
        r={"n":n,"shift":float(s),"L":float(L),"min_Q_on_ker_lc":m}
        if best is None or m<best["min_Q_on_ker_lc"]: best=dict(r,v=v,x=x)
        rows.append(r)
    print(n, best["shift"], best["L"], best["min_Q_on_ker_lc"])
json.dump({"h":h,"rows":rows},open("scan2.json","w"),indent=1)
