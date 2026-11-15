# R2.2: split cross(rho1,rho2,eps) into the two centre terms (T3-type, O(eps^2)) and the rest (T2', O(eps^4)).
import mpmath as mp, json
mp.mp.dps=40
g1=mp.im(mp.zetazero(1)); g2=mp.im(mp.zetazero(2))   # labels only
flint={'0.12':mp.mpf('0.0006070883996287'),'0.06':mp.mpf('0.0001517827770008')}
LG=lambda w: -mp.log(mp.pi)/2+mp.digamma(w/2)/2
out={'rows':[]}
for es,c in flint.items():
    e=mp.mpf(es); d=g2-g1
    near=2*e**2/(e**2+d**2)
    rest=c-near
    out['rows'].append({'eps':es,'cross_flint':mp.nstr(c,16),'centre_terms_2e2/(e2+d2)':mp.nstr(near,16),'rest':mp.nstr(rest,6),'rest/eps^4':mp.nstr(rest/e**4,6)})
Phi=lambda w: LG(w)+LG(1-w)
out['minus_half_Phi_at_i_gamma1']=mp.nstr(-Phi(mp.mpc(0,g1))/2,10)
print(json.dumps(out,indent=1)); json.dump(out,open('cross_split_results.json','w'),indent=1)
