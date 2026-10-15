"""S26 float fine-check of an LP multiplier m on a fine grid (UNTRUSTED)."""
import json,sys,numpy as np
sys.path.insert(0,'.'); from s26_lp_odd import sigA,wDC,wPL_columns
d=json.load(open(sys.argv[1])); h=d['h']; step=float(sys.argv[2]) if len(sys.argv)>2 else 0.0005
u=2*h+d['D']*np.arange(d['K']); a=np.array(d['a']); v=np.array(d['v_pl']); knots=np.array(d['pl_knots'])
def m(x):
    A0,Ap=wPL_columns(x,h,knots); w=sigA(x)+v[0]*A0+sum(vp*A for vp,A in zip(v[1:],Ap))
    for C,Dd in zip(d['C'],d['Dd_list'] if d.get('useDC') else []): w+=C*wDC(x,h,Dd)
    return w-np.cos(np.outer(x,u))@a
mins=[]
for lo in np.arange(0,d['Xi']+20,20.0):
    x=np.arange(lo,lo+20,step); mv=m(x); i=mv.argmin(); mins.append((float(mv[i]),float(x[i])))
mn=min(mins); print(json.dumps(dict(file=sys.argv[1],c_LP=d['c_LP'],fine_min=mn[0],at=mn[1],step=step,C0=float(v[0]),v=v.tolist(),knots=knots.tolist())))
