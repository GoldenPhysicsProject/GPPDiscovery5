# R2.1: J2 closed form (statement sheet §2) at eps=0.03, no zero data (zeta continuation only).
import mpmath as mp, json
mp.mp.dps=50
def Lxi(w): return mp.zeta(w,derivative=1)/mp.zeta(w)+1/w+1/(w-1)-mp.log(mp.pi)/2+mp.digamma(w/2)/2
def dLxi(w):
    z=mp.zeta(w); zp=mp.zeta(w,derivative=1); zpp=mp.zeta(w,derivative=2)
    return zpp/z-(zp/z)**2-1/w**2-1/(w-1)**2+mp.polygamma(1,w/2)/4
th=mp.mpf(1)/4
Ah=lambda q,w,e: e**2/(e**2-(w-q)**2)
g1=mp.im(mp.zetazero(1))   # used only as a LABEL
models={'rho1':(Lxi,dLxi,mp.mpc(0.5,g1)),
        's_theta':(lambda w:Lxi(w-th)+Lxi(w+th),lambda w:dLxi(w-th)+dLxi(w+th),mp.mpc(0.75,g1))}
def J2(L,dL,q,e):
    tq=1-mp.conj(q)
    B=-((e**2/4)*(dL(q+e)+dL(q-e))+(e/4)*(L(q-e)-L(q+e)))
    if abs(tq-q)<mp.mpf(10)**-40: A=B
    else:
        A=-(-(e/2)*Ah(tq,q+e,e)*L(q+e)+(e/2)*Ah(tq,q-e,e)*L(q-e)-(e/2)*Ah(q,tq+e,e)*L(tq+e)+(e/2)*Ah(q,tq-e,e)*L(tq-e))
    D=2*mp.re(A)-2*mp.re(B); det=A*A-abs(B)**2
    return A,B,D,det
flint={'rho1':{'A':'1.0000000004475105965','D':'0'},
       's_theta':{'A':'-0.0072260127662419','ReB':'1.0000130546933103','D':'-2.0144781349191','det':'-0.9999738942965'}}
flint_ratio={'0.12':'-2.2435717508','0.06':'-2.0584290569','0.03':'-2.0144518369'}
out={'label':'NUMERICAL','method':'J2 four-point formula, zeta continuation, dps 50, no zero data','rows':[]}
for nm,(L,dL,q) in models.items():
    for es in ['0.12','0.06','0.03']:
        e=mp.mpf(es); A,B,D,det=J2(L,dL,q,e)
        r={'case':nm,'eps':es,'A':mp.nstr(A,20),'B':mp.nstr(B,20),'D':mp.nstr(D,20),'det':mp.nstr(det,20)}
        if nm=='s_theta':
            r['D_over_ReB']=mp.nstr(D/mp.re(B),14); r['flint_ratio']=flint_ratio[es]
            r['ratio_diff']=mp.nstr(abs(D/mp.re(B)-mp.mpf(flint_ratio[es])),3)
        if es=='0.03':
            f=flint[nm]
            r['diff_A']=mp.nstr(abs(A-mp.mpf(f['A'])),3)
            if 'ReB' in f:
                r['diff_ReB']=mp.nstr(abs(mp.re(B)-mp.mpf(f['ReB'])),3)
                r['diff_D']=mp.nstr(abs(D-mp.mpf(f['D'])),3); r['diff_det']=mp.nstr(abs(mp.re(det)-mp.mpf(f['det'])),3)
            else: r['diff_B']=mp.nstr(abs(B-mp.mpf(f['A'])),3)
        out['rows'].append(r); print(r)
json.dump(out,open('j2_eps003_results.json','w'),indent=1)
