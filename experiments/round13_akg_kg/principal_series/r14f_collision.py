"""Round 14F: on-line zeros of Lambda_lam = zeta_K + lam*L_psi near t~14-16 for small lam (locate the collision)."""
import json, mpmath as mp
mp.mp.dps = 20
exec(open("r14f_quotient.py").read().split("# --- 1.")[0])
Lam_l = lambda lam, t: mp.re((mp.sqrt(20)/(2*mp.pi))**(0.5+1j*t)*mp.gamma(0.5+1j*t)*(zK(0.5+1j*t)+lam*Lpsi(0.5+1j*t)))
def roots(lam, a=13.8, b=16.0, n=440):
    ts = [a+(b-a)*k/n for k in range(n+1)]; v = [Lam_l(lam, t) for t in ts]
    return [round((ts[k]+ts[k+1])/2, 3) for k in range(n) if v[k]*v[k+1] < 0]
grid = {str(l): roots(mp.mpf(l)) for l in [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06]}
res = json.load(open("r14f_quotient.json")); res.pop("collision", None)
res["online_zeros_small_lambda_t_in_[13.8,16]"] = grid
json.dump(res, open("r14f_quotient.json", "w"), indent=1, default=str)
print(json.dumps(grid, indent=1))
