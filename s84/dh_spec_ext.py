"""S84: DH zero-side decomposition of the CERTIFIED DH negative witnesses. Q_DH(f) = sum over DH zeros rho=1/2+iz of fhat(z)^2
(real even f): on-line gamma>0 -> 2 fhat(gamma)^2 >= 0; off-line quadruple 1/2+-d+-i g -> 4 Re fhat(g - i d)^2 (can be < 0).
DH zeros from dh_zeros.py (T=250, mpmath). Tail beyond T from the 1/t^2 asymptotics. NUMERICAL."""
import json, mpmath as mp
from spec_decomp import fhat_factory
mp.mp.dps = 60
Z = json.load(open("dh_zeros_T250.json")); on = [mp.mpf(x) for x in Z["online"]]
import glob
for fz in sorted(glob.glob("dh_zeros_[0-9]*_*.json")): on += [mp.mpf(x) for x in json.load(open(fz))]
Z["T"] = max(int(fz.split("_")[-1].split(".")[0]) for fz in glob.glob("dh_zeros_[0-9]*_*.json")) if glob.glob("dh_zeros_[0-9]*_*.json") else Z["T"]
print("online zeros", len(on), "T", Z["T"]); off = [(mp.mpf(a), mp.mpf(b)) for a, b in Z["offline"]]
res = []
for fn in ["../s75/cosvec_DH_even_L137_80_K129.json", "../s75/cosvec_DH_even_L3421_2000_K300.json", "../s75/cosvec_DH_even_L1711_1000_K240.json"]:
    d = json.load(open(fn)); L = d['L'].replace('_','/')
    f, n2, Lm, v, fr, moms = fhat_factory(d['vec'], L, 'even')
    ton = [2*f(g)**2/n2 for g in on]
    toff = [4*mp.re(f(mp.mpc(g, -(s-mp.mpf(1)/2)))**2)/n2 for s, g in off]
    M1 = moms[0]; T = mp.mpf(Z["T"])
    tail = 2*mp.quad(lambda t: 2*M1**2/t**4*mp.log(5*t/(2*mp.pi))/(2*mp.pi), [T, mp.inf])/n2
    rq = mp.mpf(d['RQ'].strip('[').split('+/-')[0])
    S_on = mp.fsum(ton); S_off = mp.fsum(toff)
    order = sorted(range(len(ton)), key=lambda i: -ton[i])
    r = {"file": fn.split('/')[-1], "L": L, "RQ_DH_certified": mp.nstr(rq, 10),
         "online_sum(all zeros<T)": mp.nstr(S_on, 6), "offline_terms(85.70,114.16,166.48,176.70)": [mp.nstr(x, 6) for x in toff],
         "tail_est": mp.nstr(tail, 3), "spec_total": mp.nstr(S_on+S_off+tail, 8), "spec/geom": mp.nstr((S_on+S_off+tail)/rq, 6),
         "top_online(index,gamma,term)": [(i+1, mp.nstr(on[i], 8), mp.nstr(ton[i], 4)) for i in order[:5]],
         "online_first5_terms": [mp.nstr(x, 4) for x in ton[:5]],
         "|fhat|/|f| at offline zeros": [mp.nstr(abs(f(mp.mpc(g, -(s-mp.mpf(1)/2))))/mp.sqrt(n2), 4) for s, g in off]}
    res.append(r); print(json.dumps(r, indent=1))
json.dump(res, open("dh_spec_results_ext.json", "w"), indent=1)
