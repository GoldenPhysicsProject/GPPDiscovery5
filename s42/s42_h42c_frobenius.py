"""Grok S42 H42c: re-express S36/S38 PD residual bounds in FROBENIUS form for GPPVerify #371
(GppFrobeniusResidual.abs_quadratic_le_of_sum_sq_le: sum E_ij^2 <= eps^2 => |x^T E x| <= eps |x|^2, no symmetry needed).
From stored certificate numbers only (no recomputation):
  cholres residual: normE = max_i sum_j |E_ij| (Arb upper)  =>  sum_ij E_ij^2 <= sum_i (sum_j|E_ij|)^2 <= K normE^2  => eps_res = sqrt(K) normE
  string error   : delta_str = K * max|entry err|          =>  ||.||_F <= K max  (already Frobenius-valid)
  quadrature     : normE_bound = K * quad_err_entry          =>  Frobenius-valid
  eps_B          : analytic Cauchy-Schwarz (HS-type) cross-block bound, unchanged (outside #371 scope).
Upward safety factor (1+1e-9) on every float read back from JSON. Directed-rounding interval arithmetic via mpmath iv.
Status: CONFIRMED (arithmetics on certificate numbers); the numbers themselves remain Arb/python-flint, NOT Lean. RH OPEN."""
import json, mpmath as mp
iv=mp.iv; iv.dps=40; up=lambda x: iv.mpf(x)*iv.mpf('1.000000001')
runs={'S36':'../s36/cert_L13_T700_K740_FINAL.json','S38':'../s38/run/cert_L13_T700_K740_D8x2_FINAL.json'}
out={}
for tag,fn in runs.items():
    d=json.load(open(fn)); K=d['K']
    for par in ['even','odd']:
        e=d[par]; c=e['cholres']; l0=iv.mpf(str(c['l0']))
        eps_res=iv.sqrt(iv.mpf(K))*up(c['normE']); eps_F=eps_res+up(c['delta_str'])+up(e['normE_bound'])
        lamA=l0-eps_F; lamD=iv.mpf(str(e['lamD_lower']))*iv.mpf('0.999999999')
        cert=mp.mpf(min(lamA.a,lamD.a))-mp.mpf(up(e['eps_B']).b)
        out[f'{tag}_{par}']={'K':K,'l0':str(c['l0']),'normE_rowsum':c['normE'],'eps_res_Frobenius':mp.nstr(eps_res.b,6),
           'eps_F_total(res+str+quad)':mp.nstr(eps_F.b,6),'eps_B':e['eps_B'],'old_certified_lower':e['certified_lower'],
           'Frobenius_certified_lower':mp.nstr(cert,16),'positive':bool(cert>0),
           'reported_floor_still_holds':bool(cert>=mp.mpf(str(c['l0']))*mp.mpf('0.99999'))}
json.dump(out,open('h42c_frobenius.json','w'),indent=1); print(json.dumps(out,indent=1))
