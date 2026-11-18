"""Round 14F follow-up: coefficient identity (normalised), Hecke test on a=(r_Q1)/2, collision lambda*."""
import json, mpmath as mp
mp.mp.dps = 25
exec(open("r14f_quotient.py").read().split("# --- 1.")[0])
Q1 = lambda x, y: x*x+5*y*y; Q2 = lambda x, y: 2*x*x+2*x*y+3*y*y
def rQ(Q, n, R=30):
    return sum(1 for x in range(-R, R+1) for y in range(-R, R+1) if Q(x, y) == n)
b1 = lambda n: sum(chi_m20(d) for d in range(1, n+1) if n % d == 0)
b2 = lambda n: sum(chi_m4(d)*chi_5(n//d) for d in range(1, n+1) if n % d == 0)
res = json.load(open("r14f_quotient.json"))
res["epstein_coeff_identity_n<=200_failures"] = [n for n in range(1, 201) if rQ(Q1, n) != b1(n)+b2(n) or rQ(Q2, n) != b1(n)-b2(n)]
res["epstein_coeff_identity_note"] = "r_Q1(n)=b1(n)+b2(n), r_Q2(n)=b1(n)-b2(n); b1=sum_{d|n}chi_-20(d), b2=sum_{d|n}chi_-4(d)chi_5(n/d)"
a = lambda n: mp.mpf(rQ(Q1, n))/2  # normalised a(1)=1
def hd(b, p, n): return b(p)*b(n) - b(p*n) - (chi_m20(p)*b(n//p) if n % p == 0 else 0)
res["hecke_defect"] = {f"p={p},n={n}": {"b1 (zeta_K)": hd(b1, p, n), "b2 (L_psi)": hd(b2, p, n), "r_Q1/2": float(hd(a, p, n))}
                       for p, n in [(3, 3), (3, 7), (7, 3), (2, 3), (7, 7), (23, 3), (3, 23)]}
# collision lambda*: count on-line zeros of Lam_lam in t in [14.6,16.2]
Lam_l = lambda lam, t: mp.re((mp.sqrt(20)/(2*mp.pi))**(0.5+1j*t)*mp.gamma(0.5+1j*t)*(zK(0.5+1j*t)+lam*Lpsi(0.5+1j*t)))
def count(lam, a=14.6, b=16.2, n=320):
    v = [Lam_l(lam, a+(b-a)*k/n) for k in range(n+1)]
    return sum(1 for k in range(n) if v[k]*v[k+1] < 0)
lo, hi = mp.mpf(0.0), mp.mpf(0.2)
cl, ch = count(lo), count(hi)
for _ in range(14):
    mid = (lo+hi)/2
    if count(mid) == cl: lo = mid
    else: hi = mid
res["collision"] = {"count_at_0": cl, "count_at_0.2": ch, "lambda*_bracket": [float(lo), float(hi)]}
json.dump(res, open("r14f_quotient.json", "w"), indent=1, default=str)
print(json.dumps({k: res[k] for k in ["epstein_coeff_identity_n<=200_failures", "hecke_defect", "collision"]}, indent=1, default=str))
