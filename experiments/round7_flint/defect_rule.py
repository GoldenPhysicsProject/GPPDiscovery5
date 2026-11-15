import mpmath as mp, json
mp.mp.dps=30
g1=mp.im(mp.zetazero(1))
LG=lambda w: -mp.log(mp.pi)/2+mp.digamma(w/2)/2
Phi=lambda w: LG(w)+LG(1-w)
Lxi=lambda w: mp.zeta(w,derivative=1)/mp.zeta(w)+1/w+1/(w-1)+LG(w)
Ah=lambda q,w,e: e**2/(e**2-(w-q)**2)
def W(w,e,uk,uj): return 2*Ah(uk,w,e)*Ah(uj,w,e)       # = (1/2) A_k A_j (peak-2 corpus weight)
def res(f,p,r=mp.mpf('1e-3')):
    return mp.quad(lambda th: f(p+r*mp.expj(th))*r*mp.expj(th),[0,2*mp.pi])/(2*mp.pi)
def truesum(e,uk,uj):   # J1: sum_w m_w W(w) = -sum_poles Res(W*Lxi)
    f=lambda w: W(w,e,uk,uj)*Lxi(w)
    ps={uk+e,uk-e,uj+e,uj-e}
    return -sum(res(f,p) for p in ps)
def defect_pred(e,uk,uj,s):
    f=lambda w: W(w,e,uk,uj)*Phi(w)
    ps=[]; tot=0
    for p in {uk+e,uj+e}:
        if mp.re(p)<0.5: tot+=res(f,p)
    for p in {uk-e,uj-e}:
        if mp.re(p)>0.5: tot+=s*res(f,p)
    return tot
e=mp.mpf('0.12')
cases={'iG1_diag':(mp.mpc(0,g1),mp.mpc(0,g1),None),
 'T1_half+20i':(mp.mpc(0.5,20),mp.mpc(0.5,20),mp.mpc('0.00018544496360264223',0)),
 'T2_tau_pair':(mp.mpc(0,g1),mp.mpc(1,g1),mp.mpc('0.00557958855898004609',0)),
 'T3_quarter':(mp.mpc(0.25,g1),mp.mpc(0.25,g1),mp.mpc('0.0653062039433421984','-0.000530785479988416697')),
 'rho1':(mp.mpc(0.5,g1),mp.mpc(0.5,g1),mp.mpc('1.0000001145048385914',0))}
out={}
for k,(uk,uj,asm) in cases.items():
    T=truesum(e,uk,uj)/2  # G units = C/2
    row={'true_G':mp.nstr(T,14)}
    for s in (1,-1):
        P=T+defect_pred(e,uk,uj,s)/2
        row['pred_assembly_G_s%+d'%s]=mp.nstr(P,14)
        if asm is not None: row['abs(pred-assembly)_s%+d'%s]=mp.nstr(abs(P-asm),3)
    if asm is not None: row['assembly_G']=mp.nstr(asm,14); row['abs(true-assembly)']=mp.nstr(abs(T-asm),3)
    out[k]=row; print(k,row)
json.dump(out,open('defect_rule_results.json','w'),indent=1)
