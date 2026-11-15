import numpy as np, json
exec(open('t2rig.py').read().split("out=dict")[0])
Q=lambda w: R(w)+Pi(w)
ref73=[-1.76852820e-02,5.37143649e-03,-3.62446009e-02,2.92743995e-02,-1.68687287e-02,1.35280159e-02,5.17525298e-03,2.86638282e-02,-3.55069657e-02,1.65808977e-02,-1.54761376e-02,-8.05459830e-04]
j=73; wj=om(j); Qj=Q(wj); res=[]
for i in range(1,13):
    wi=om(i); v=(16/np.pi)*(2/l)*(-1)**(i+j)*wi*wj/(wj**2-wi**2)*(Qj/wj-Q(wi)/wi)-(bI2**2/8)*sig(i)*sig(j)
    res.append([i,float(v),ref73[i-1]]); print(i,v,ref73[i-1])
json.dump(res,open('phicheck_results.json','w'),indent=1)
