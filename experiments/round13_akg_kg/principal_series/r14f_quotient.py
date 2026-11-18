"""Round 14F checks (single self-contained script; output does not depend on run order):
Epstein x^2+5y^2 genus decomposition, Euler-product / Hecke tests, DH root number and Z5^x non-equivariance,
interference at off-line zeros, homotopy Z_lam = zeta_K + lam L_psi and the collision bracket.
Zeros are test data only. Run: python3 r14f_quotient.py -> r14f_quotient.json"""
import json, mpmath as mp
mp.mp.dps = 20
out = {}
def L(s, q, ch):
    return q**(-s)*mp.fsum(ch(a)*mp.zeta(s, mp.mpf(a)/q) for a in range(1, q+1) if ch(a) != 0)
chi_m4 = lambda a: [0, 1, 0, -1][a % 4]
chi_5 = lambda a: [0, 1, -1, -1, 1][a % 5]
chi_m20 = lambda a: chi_m4(a)*chi_5(a)
zK = lambda s: mp.zeta(s)*L(s, 20, chi_m20)            # zeta_K, K=Q(sqrt-5): 1 [+] chi_-20
Lpsi = lambda s: L(s, 4, chi_m4)*L(s, 5, chi_5)        # genus character psi: chi_-4 [+] chi_5
Q1 = lambda x, y: x*x+5*y*y; Q2 = lambda x, y: 2*x*x+2*x*y+3*y*y
# --- 1. constants: Z_Q1 = zK + Lpsi, Z_Q2 = zK - Lpsi
def lattice(Q, s, R=400):
    return sum(Q(x, y)**(-s) for x in range(-R, R+1) for y in range(-R, R+1) if (x, y) != (0, 0))
out["epstein_constants_s=3"] = {"Z_Q1 lattice(R=400)": lattice(Q1, 3), "zK+Lpsi": float(zK(3)+Lpsi(3)),
                                "Z_Q2 lattice(R=400)": lattice(Q2, 3), "zK-Lpsi": float(zK(3)-Lpsi(3))}
def rQ(Q, n, R=30):
    return sum(1 for x in range(-R, R+1) for y in range(-R, R+1) if Q(x, y) == n)
b1 = lambda n: sum(chi_m20(d) for d in range(1, n+1) if n % d == 0)
b2 = lambda n: sum(chi_m4(d)*chi_5(n//d) for d in range(1, n+1) if n % d == 0)
out["epstein_coeff_identity"] = {"statement": "r_Q1(n)=b1(n)+b2(n), r_Q2(n)=b1(n)-b2(n); b1=1*chi_-20, b2=chi_-4*chi_5",
    "failures_n<=200": [n for n in range(1, 201) if rQ(Q1, n) != b1(n)+b2(n) or rQ(Q2, n) != b1(n)-b2(n)]}
a = lambda n: rQ(Q1, n)/2                                # normalised a(1)=1
out["epstein_non_euler"] = {"a(2)": a(2), "a(3)": a(3), "a(6)": a(6), "a(2)a(3)": a(2)*a(3)}
def hd(b, p, n): return b(p)*b(n) - b(p*n) - (chi_m20(p)*b(n//p) if n % p == 0 else 0)
out["hecke_defect b(p)b(n)-b(pn)-chi_-20(p)b(n/p)"] = {f"p={p},n={n}": {"b1": hd(b1, p, n), "b2": hd(b2, p, n), "r_Q1/2": hd(a, p, n)}
    for p, n in [(3, 3), (3, 7), (7, 3), (2, 3), (7, 7), (23, 3), (3, 23)]}
# --- 2. DH
chi = {0: 0, 1: 1, 2: 1j, 4: -1, 3: -1j}
tau = mp.fsum(chi[k]*mp.exp(2j*mp.pi*k/5) for k in range(1, 5))
eps = tau/(1j*mp.sqrt(5)); epsb = -mp.conj(tau)/(1j*mp.sqrt(5))   # tau(chibar) = chi(-1) conj(tau(chi))
kap = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
al, be = (1-1j*kap)/2, (1+1j*kap)/2
out["DH_rootnumbers"] = {"eps(chi)": str(mp.nstr(eps, 12)), "beta/alpha": str(mp.nstr(be/al, 12)),
                         "1/eps(chibar)": str(mp.nstr(1/epsb, 12)), "|beta/alpha|": float(abs(be/al))}
Lc = lambda s: L(s, 5, lambda k: chi[k % 5]); Lcb = lambda s: L(s, 5, lambda k: mp.conj(chi[k % 5]))
Lam = lambda f, s: (mp.pi/5)**(-(s+1)/2)*mp.gamma((s+1)/2)*f(s)
s0 = mp.mpc(0.3, 4.0)
out["FE_chi_check_rel"] = float(abs(Lam(Lc, s0) - eps*Lam(Lcb, 1-s0))/abs(Lam(Lc, s0)))
c = {k: (al*chi[k]+be*mp.conj(chi[k])).real for k in range(5)}
out["DH_coeffs"] = {k: float(c[k]) for k in range(5)}
out["DH_translation_ratio_c(2a)/c(a)"] = {k: float(c[(2*k) % 5]/c[k]) for k in range(1, 5)}
DH = lambda s: al*Lc(s)+be*Lcb(s)
# --- 3. interference at off-line zeros (equal moduli FORCED by |beta/alpha|=1; the sizes are the content)
r = mp.findroot(DH, mp.mpc(0.8085, 85.6993))
out["DH_offline_zero"] = {"rho": str(mp.nstr(r, 12)), "|L(chi,rho)|": float(abs(Lc(r))), "|L(chibar,rho)|": float(abs(Lcb(r))),
                          "ratio": str(mp.nstr(Lc(r)/Lcb(r), 12)), "-beta/alpha": str(mp.nstr(-be/al, 12))}
re_ = mp.findroot(lambda s: zK(s)+Lpsi(s), mp.mpc(0.93297, 15.66825))
out["EPS_offline_zero"] = {"rho": str(mp.nstr(re_, 12)), "|zeta_K(rho)|": float(abs(zK(re_))), "|L_psi(rho)|": float(abs(Lpsi(re_))),
                           "ratio": str(mp.nstr(zK(re_)/Lpsi(re_), 12))}
z1 = mp.zetazero(1)
out["common_zero_test |L_psi| at first zeta zero"] = float(abs(Lpsi(z1)))
# --- 4. homotopy and collision
track, rho, lam = [], re_, mp.mpf(1)
while lam > mp.mpf(0.1449):
    rho = mp.findroot(lambda s: zK(s)+lam*Lpsi(s), rho)
    track.append([float(lam), float(rho.real), float(rho.imag)])
    lam -= mp.mpf(1)/20 if lam > 0.25 else mp.mpf(1)/200
out["homotopy_offline_zero [lam, Re, Im]"] = track
LamR = lambda lam, t: mp.re((mp.sqrt(20)/(2*mp.pi))**(0.5+1j*t)*mp.gamma(0.5+1j*t)*(zK(0.5+1j*t)+lam*Lpsi(0.5+1j*t)))
def online(lam, a=14.5, b=16.2, h=0.005):
    n = int(round((b-a)/h)); ts = [a+h*k for k in range(n+1)]; v = [LamR(lam, t) for t in ts]
    return [round(ts[k]+h/2, 4) for k in range(n) if v[k]*v[k+1] < 0]
out["online_zeros_t_in_[14.5,16.2]_step0.005"] = {str(l): online(mp.mpf(l)) for l in
    ["0", "0.06", "0.1", "0.13", "0.135", "0.14", "0.1425", "0.145", "0.15", "0.16"]}
json.dump(out, open("r14f_quotient.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
