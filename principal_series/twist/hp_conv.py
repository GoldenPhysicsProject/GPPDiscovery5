"""basis-convergence check of zeta lam_even, kappa at N=40, dps 80 (HIGH-PRECISION NUMERICAL)."""
import sys; sys.argv=['x','40']; sys.path.insert(0,'/workspace/newbot/twist')
import hp_onset_lib as H, mpmath as mp, json
res=[]
for L in ['0.6','0.8','1.0','1.2']:
    le,ka,lo=H.zeta_data(L); r={"L":L,"N":40,"lam_even":mp.nstr(le,10),"kappa":mp.nstr(ka,10),"lam_odd":mp.nstr(lo,10)}; print(r,flush=True); res.append(r)
json.dump(res,open('hp_conv_N40.json','w'),indent=1)
