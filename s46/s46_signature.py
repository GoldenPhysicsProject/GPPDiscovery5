"""Signature check (NUMERICAL, mp eig): planted quartet matrix P = int F_jk 4cos(gamma x)cosh(delta x) dx alone,
per parity sector. Identity (PROVED): 4cos(g(t-s))cosh(d(t-s)) = 4[ C(t)C(s) + S'(t)S'(s) ] - 4[ Sh(t)Sh(s) + Ch(t)Ch(s) ]
with even sector: +cos(gt)cosh(dt), -sin(gt)sinh(dt);  odd sector: +sin(gt)cosh(dt), -cos(gt)sinh(dt).
So each sector has signature exactly (1,1) for delta>0 (rank 2), total (2,2) [cf. GPPVerify NoGo/QuartetSignature.lean]; delta=0 -> rank 1 PSD per sector."""
import sys, json; sys.path.insert(0,'.')
from s46_planted import *
ctx.prec=1024; mp.mp.dps=60; M=24; Gs=Gpolys(M); out=[]
for h in ['0.5','1.3']:
  for g,d in [('17.578','0.308'),('85.699','0.308'),('85.699','0.0')]:
    P=delta_mat(M,h,Gs,[(arb(2),acb(0,1)*arb(g)+arb(d)),(arb(2),acb(0,1)*arb(g)-arb(d))],zerr=None)
    row={"h":h,"gamma":g,"delta":d}
    for nm,I in (("even",list(range(0,M,2))),("odd",list(range(1,M,2)))):
        ev=mp.eigsy(tomp(P,I),eigvals_only=True); mx=max(abs(e) for e in ev)
        row[nm]={"n_pos":sum(1 for e in ev if e>mx*1e-25),"n_neg":sum(1 for e in ev if e<-mx*1e-25),"neg_eig":mp.nstr(min(ev),5),"pos_eig":mp.nstr(max(ev),5)}
    print(json.dumps(row)); out.append(row)
json.dump(out,open('signature.json','w'),indent=1)
