import mpmath as mp
mp.mp.dps=30
g1=mp.im(mp.zetazero(1))
LG=lambda w: -mp.log(mp.pi)/2+mp.digamma(w/2)/2
Phis={'LG(w)':LG,'LG(1-w)':lambda w:LG(1-w),'LG(w)+LG(1-w)':lambda w:LG(w)+LG(1-w),'LG(w)-LG(1-w)':lambda w:LG(w)-LG(1-w)}
def res(u,e,p,Phi):  # Res_p of 2*Ahat_u^2*Phi ; principal part of Ahat^2 at u+-e: (e^2/4)x^-2 -+ (e/4) x^-1
    sgn = -1 if p=='+' else 1
    z = u+e if p=='+' else u-e
    return 2*((e**2/4)*mp.diff(Phi,z) + sgn*(e/4)*Phi(z))
data=[ (mp.mpc(0,g1),mp.mpf('0.12'), mp.mpc(-0.0411966679324811299/2-0.0037358301982215840, -0.00212250234432643148/2-2.9372743746e-8)),
       (mp.mpc(0,g1),mp.mpf('0.06'), mp.mpc(-0.0239080333601461565/2-0.00021347032370099371, -0.0010612301742963499/2-1.8369330028e-9)),
       (mp.mpc(0.05,g1),mp.mpf('0.12'), mp.mpc(-0.024330597,-0.00095522707)),
       (mp.mpc(-0.1,g1),mp.mpf('0.12'), mp.mpc(-0.024342422,-0.0012733052))]
best=[]
for nm,Phi in Phis.items():
  for S in (['+'],['-'],['+','-']):
    for k in (1,-1,mp.mpf(1)/2,-mp.mpf(1)/2,2,-2):
      err=max(abs(k*sum(res(u,e,p,Phi) for p in S)-gap) for u,e,gap in data)
      best.append((float(err),nm,S,float(k)))
best.sort(); [print(b) for b in best[:6]]
for u,e,gap in data: print(mp.nstr(gap/e,8))
