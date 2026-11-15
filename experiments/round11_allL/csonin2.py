import numpy as np, json, scipy.linalg as sl
exec(open('csonin.py').read().split('rows = []')[0])
out=[]
for L in (0.5,1.0,2.0):
    n=41; xs=np.linspace(0,L,n)
    A=toe(arch,L,xs); P=toe(primeS,L,xs); Pol=toe(pole,L,xs)
    N=np.array([[triG(xs[i]-xs[j]) for j in range(n)] for i in range(n)])
    Q=sl.null_space(np.cosh(xs/2)[None,:]); r=lambda M:Q.T@M@Q
    ev,V=sl.eigh(r(A),r(N)); c=Q@V[:,0]
    ts=np.linspace(0,40,4001); F=np.abs((c[None,:]*np.exp(1j*ts[:,None]*xs[None,:])).sum(1))**2*(np.sinc(ts/(2*np.pi)))**4
    tpk=float(ts[np.argmax(F)]); frac_low=float(F[ts<6.29].sum()/F.sum())
    # cost of the arch-negative direction in the other terms (Rayleigh quotients wrt N)
    nN=c@N@c
    out.append(dict(L=L,n=n,neg_inertia_A=int((ev<0).sum()),lam_A=[float(e) for e in ev[:4]],
        bottom_vec_peak_t=tpk,bottom_vec_weight_below_6p29=frac_low,
        on_bottom_vec=dict(A=float(c@A@c/nN),minus_P=float(-c@P@c/nN),pole=float(c@Pol@c/nN),total=float(c@(A-P+Pol)@c/nN))))
    print(json.dumps(out[-1]),flush=True)
json.dump(out,open('csonin2_results.json','w'),indent=1)
