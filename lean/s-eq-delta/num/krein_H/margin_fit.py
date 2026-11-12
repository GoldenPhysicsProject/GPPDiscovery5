# Turbo 6821237527 item 2: fit decay of the ker-l_c margin (margin.json best_per_n).
import numpy as np, json
d=json.load(open("margin.json")); b=d["best_per_n"]; h=d["h"]
n=np.array([r["n"] for r in b]); m=np.array([r["min_Q_on_ker_lc"] for r in b]); L=np.array([r["L"] for r in b])
S=(n-1)*h  # span of the translate grid; supports fill a window of length S+2
out={"fits":{}}
def fit(name,X,mask):
    A=np.vstack([np.ones(mask.sum()),X[mask]]).T; y=np.log(m[mask])
    c,res,_,_=np.linalg.lstsq(A,y,rcond=None); r=y-A@c
    out["fits"][name]={"intercept":c[0],"slope":c[1],"rms_resid_log":float(np.sqrt(np.mean(r**2))),"max_resid_log":float(np.max(np.abs(r))),"n_range":[int(n[mask][0]),int(n[mask][-1])]}
    print(name,out["fits"][name])
for lo,hi in [(13,57),(13,25),(25,41),(41,57),(25,57)]:
    k=(n>=lo)&(n<=hi)
    fit(f"exp_vs_span_n{lo}-{hi}",S,k); fit(f"power_vs_span_n{lo}-{hi}",np.log(S),k)
k=(n>=13); fit("exp_vs_L_bestshift_n13-57",L,k)
# shift-0 (symmetric grid, L = 1+S/2) series, so L and S are tied
rows=[r for r in d["rows"] if r["shift"]==0.0]
m0=np.array([r["min_Q_on_ker_lc"] for r in rows]); L0=np.array([r["L"] for r in rows]); n0=np.array([r["n"] for r in rows])
k=n0>=13; A=np.vstack([np.ones(k.sum()),L0[k]]).T; c=np.linalg.lstsq(A,np.log(m0[k]),rcond=None)[0]; r=np.log(m0[k])-A@c
out["fits"]["exp_vs_L_shift0_n13-57"]={"slope":c[1],"rms_resid_log":float(np.sqrt(np.mean(r**2)))}; print("exp_vs_L_shift0",c[1],np.sqrt(np.mean(r**2)))
A=np.vstack([np.ones(k.sum()),np.log(L0[k])]).T; c=np.linalg.lstsq(A,np.log(m0[k]),rcond=None)[0]; r=np.log(m0[k])-A@c
out["fits"]["power_vs_L_shift0_n13-57"]={"slope":c[1],"rms_resid_log":float(np.sqrt(np.mean(r**2)))}; print("power_vs_L_shift0",c[1],np.sqrt(np.mean(r**2)))
# local rates over 4-point blocks (span +1)
loc=[]
for i in range(len(n)-4):
    if n[i]>=13: loc.append({"n":int(n[i]),"n+4":int(n[i+4]),"exp_rate_per_unit_span":float(-np.log(m[i+4]/m[i])/(S[i+4]-S[i])),"power_exponent":float(-np.log(m[i+4]/m[i])/np.log(S[i+4]/S[i]))})
out["local"]=loc
for l in loc[::4]: print(l)
json.dump(out,open("margin_fit.json","w"),indent=1)
