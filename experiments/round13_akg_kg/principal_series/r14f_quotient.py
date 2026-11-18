"""Round 14F checks: isotypic decomposition of the Epstein (x^2+5y^2) and Davenport-Heilbronn quotients.
Zeros are test data only. Run: python3 r14f_quotient.py -> r14f_quotient.json"""
import json, mpmath as mp
mp.mp.dps = 30
out = {}
def L(s, q, ch):
    return q**(-s)*mp.fsum(ch(a)*mp.zeta(s, mp.mpf(a)/q) for a in range(1, q+1) if ch(a) != 0)
chi_m4 = lambda a: [0, 1, 0, -1][a % 4]
chi_5 = lambda a: [0, 1, -1, -1, 1][a % 5]
chi_m20 = lambda a: chi_m4(a)*chi_5(a)
zK = lambda s: mp.zeta(s)*L(s, 20, chi_m20)            # zeta_K, K=Q(sqrt-5): isobaric 1 [+] chi_-20
Lpsi = lambda s: L(s, 4, chi_m4)*L(s, 5, chi_5)        # genus character psi: chi_-4 [+] chi_5
# --- 1. exact constants: Z_Q1 = zK + Lpsi, Z_Q2 = zK - Lpsi (w=2, h=2: Z_Q = w*zeta_K(s,C) = zK +/- Lpsi)
def lattice(Q, s, R=400):
    return mp.fsum(Q(x, y)**(-s) for x in range(-R, R+1) for y in range(-R, R+1) if (x, y) != (0, 0))
Q1 = lambda x, y: x*x+5*y*y; Q2 = lambda x, y: 2*x*x+2*x*y+3*y*y
s = 3
out["epstein_constants_s=3"] = {"Z_Q1 lattice(R=400)": float(lattice(Q1, s)), "zK+Lpsi": float(zK(s)+Lpsi(s)),
                                "Z_Q2 lattice(R=400)": float(lattice(Q2, s)), "zK-Lpsi": float(zK(s)-Lpsi(s))}
# coefficient identity r_Q1(n)/2 = b1(n)+b2(n), r_Q2(n)/2 = b1-b2, n<=200
def rQ(Q, n, R=30):
    return sum(1 for x in range(-R, R+1) for y in range(-R, R+1) if Q(x, y) == n)
b1 = lambda n: sum(chi_m20(d) for d in range(1, n+1) if n % d == 0)
b2 = lambda n: sum(chi_m4(d)*chi_5(n//d) for d in range(1, n+1) if n % d == 0)
bad = [n for n in range(1, 201) if rQ(Q1, n) != 2*(b1(n)+b2(n)) or rQ(Q2, n) != 2*(b1(n)-b2(n))]
out["epstein_coeff_identity_n<=200_failures"] = bad
# --- 2. Hecke-eigen test (weight-1 level-20 nebentypus chi_-20): b(p)b(n) = b(pn) + chi_-20(p) b(n/p)
def hecke_defect(b, p, n):
    return b(p)*b(n) - b(p*n) - (chi_m20(p)*b(n//p) if n % p == 0 else 0)
a1 = lambda n: b1(n)+b2(n)
rows = {}
for p, n in [(3, 3), (3, 7), (7, 3), (2, 3), (3, 9), (7, 7), (29, 3)]:
    rows[f"p={p},n={n}"] = {"zeta_K piece": hecke_defect(b1, p, n), "L_psi piece": hecke_defect(b2, p, n), "Z_Q1/2": hecke_defect(a1, p, n)}
out["hecke_defect"] = rows
# --- 3. DH: chi mod 5, chi(2)=i. root number, symmetric combination
chi = {0: 0, 1: 1, 2: 1j, 4: -1, 3: -1j}
cx = lambda a: chi[a % 5]; cxb = lambda a: mp.conj(chi[a % 5])
tau = mp.fsum(chi[a]*mp.exp(2j*mp.pi*a/5) for a in range(1, 5))
eps = tau/(1j*mp.sqrt(5)); epsb = mp.conj(tau)/(1j*mp.sqrt(5)) * (-1)  # tau(chibar)=chi(-1)conj(tau(chi)) = -conj(tau)
kap = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
al, be = (1-1j*kap)/2, (1+1j*kap)/2
out["DH_rootnumbers"] = {"eps(chi)": str(mp.nstr(eps, 12)), "eps(chibar)": str(mp.nstr(epsb, 12)),
                         "beta/alpha": str(mp.nstr(be/al, 12)), "1/eps(chibar)": str(mp.nstr(1/epsb, 12)), "|beta/alpha|": float(abs(be/al))}
Lc = lambda s: L(s, 5, cx); Lcb = lambda s: L(s, 5, cxb)
DH = lambda s: al*Lc(s)+be*Lcb(s)
# FE check of completed pieces: Lam(chi,s) = eps Lam(chibar,1-s)
Lam = lambda f, s: (mp.pi/5)**(-(s+1)/2)*mp.gamma((s+1)/2)*f(s)
s0 = mp.mpc(0.3, 4.0)
out["FE_chi_check"] = float(abs(Lam(Lc, s0) - eps*Lam(Lcb, 1-s0))/abs(Lam(Lc, s0)))
# non-eigen under (Z/5)^x translation: c(2a)/c(a)
c = {a: (al*chi[a]+be*mp.conj(chi[a])).real for a in range(5)}
out["DH_translation_ratio_c(2a)/c(a)"] = {a: float(c[(2*a) % 5]/c[a]) for a in range(1, 5)}
out["chi_translation_ratio"] = {a: str(chi[(2*a) % 5]/chi[a]) for a in range(1, 5)}
# --- 4. interference at off-line zeros
r = mp.findroot(DH, mp.mpc(0.8085, 85.6993))
out["DH_offline_zero"] = {"rho": str(mp.nstr(r, 12)), "|L(chi,rho)|": float(abs(Lc(r))), "|L(chibar,rho)|": float(abs(Lcb(r))),
                          "L(chi)/L(chibar) at rho": str(mp.nstr(Lc(r)/Lcb(r), 12)), "-beta/alpha": str(mp.nstr(-be/al, 12))}
re = mp.findroot(lambda s: zK(s)+Lpsi(s), mp.mpc(0.93297, 15.66825))
out["EPS_offline_zero"] = {"rho": str(mp.nstr(re, 12)), "|zeta_K(rho)|": float(abs(zK(re))), "|L_psi(rho)|": float(abs(Lpsi(re))),
                           "zeta_K/L_psi at rho": str(mp.nstr(zK(re)/Lpsi(re), 12))}
# pieces at their own zeros: is Z_Q zero there? (common-zero test)
z1 = mp.zetazero(1)
out["common_zero_test"] = {"|L_psi| at first zeta zero": float(abs(Lpsi(z1))),
                           "|Z_Q1| at first zeta zero": float(abs(zK(z1)+Lpsi(z1)))}
# --- 5. homotopy Z_lam = zeta_K + lam*L_psi (FE, gamma, conductor, real coeffs preserved for real lam)
Lam_l = lambda lam, s: (mp.sqrt(20)/(2*mp.pi))**s*mp.gamma(s)*(zK(s)+lam*Lpsi(s))
track = []; rho = re
lam = mp.mpf(1)
while lam > 0.0:
    try:
        rho = mp.findroot(lambda s: zK(s)+lam*Lpsi(s), rho)
    except Exception:
        break
    track.append([float(lam), float(rho.real), float(rho.imag)])
    if abs(rho.real-0.5) < 2e-3: break
    lam -= mp.mpf(1)/50 if rho.real-0.5 > 0.08 else mp.mpf(1)/400
out["homotopy_offline_zero"] = track
# on-line zeros of Lam_lam in t in [12,19] for several lam (real on the line)
def online(lam, a=12, b=19, n=1400):
    ts = [a+(b-a)*k/n for k in range(n+1)]
    v = [mp.re(Lam_l(lam, mp.mpf(0.5)+1j*t)) for t in ts]
    roots = []
    for k in range(n):
        if v[k]*v[k+1] < 0:
            roots.append(float(mp.findroot(lambda t: mp.re(Lam_l(lam, mp.mpf(0.5)+1j*t)), (ts[k], ts[k+1]), solver='anderson')))
    return roots
out["online_zeros_t_in_[12,19]"] = {str(l): online(l) for l in [0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]}
json.dump(out, open("r14f_quotient.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
