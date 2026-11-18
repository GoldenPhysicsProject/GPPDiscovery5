"""Round 14A checks: adelic Pasterski-Shao bulk (Tate / Weil representation).
Zeros are used only as test data, never as inputs to any argument.
Run: python3 r14a_adelic_ps.py  -> r14a_adelic_ps.json
"""
import json, mpmath as mp
mp.mp.dps = 30
out = {}

# 1. Tate integral of f0 = exp(-pi x^2) (x) 1_Zhat equals Lambda(s) = pi^{-s/2} Gamma(s/2) zeta(s)
#    Z(f0,s) = int_{C_Q} E-sum = int_0^inf (theta(t)-1)/2 * t^{s/2} dt/t,  theta(t)=sum_n exp(-pi n^2 t)
def tate_f0(s):
    def thf(t):  # full theta; modular relation theta(t)=t^{-1/2}theta(1/t) for small t (Poisson on Z)
        return mp.jtheta(3, 0, mp.exp(-mp.pi*t)) if t >= 0.2 else thf(1/t)/mp.sqrt(t)
    th = lambda t: thf(t) - 1
    return mp.quad(lambda t: th(t)/2 * t**(s/2-1), [0, 1, mp.inf])
for s in [mp.mpc(2.3, 1.7), mp.mpc(1.5, 6.0)]:
    a, b = tate_f0(s), mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    out[f"tate_vs_Lambda_s={complex(s)}"] = {"tate": str(mp.nstr(a, 15)), "Lambda": str(mp.nstr(b, 15)), "relerr": float(abs(a-b)/abs(b))}

# 2. Plancherel on C_Q (unramified, R_+ part): f = g - 2 g(2x), g = x^2 exp(-pi x^2), f_p = 1_Zp.
#    f(0)=0, fhat(0)=0 so E f in L^2(R_+, dt/t).  E f(t) = t^{1/2} sum_{n!=0} f(t n).
#    Mellin: M(s) = int E f(t) t^{s-1/2} dt/t = zeta(s)(1-2^{1-s}) pi^{-(s+2)/2} Gamma(s/2+1).
#    Check int |Ef|^2 dt/t  ==  (1/2pi) int |M(1/2+i tau)|^2 dtau.
g = lambda x: x**2*mp.exp(-mp.pi*x**2)
f = lambda x: g(x) - 2*g(2*x)
def Ef(t):
    N = int(12/t) + 5
    return mp.sqrt(t) * 2 * mp.fsum(f(t*n) for n in range(1, N))
lhs = mp.quad(lambda u: Ef(mp.exp(u))**2, [-4, -2, 0, 2, 4])  # dt/t = du
M = lambda s: mp.zeta(s)*(1-2**(1-s))*mp.pi**(-(s+2)/2)*mp.gamma(s/2+1)
rhs = mp.quad(lambda tau: abs(M(mp.mpf(0.5)+1j*tau))**2, mp.linspace(-80, 80, 33))/(2*mp.pi)
# Mellin of Ef itself vs formula at one point (checks E f <-> Tate integral)
s0 = mp.mpc(0.5, 3.0)
mel = mp.quad(lambda u: Ef(mp.exp(u))*mp.exp((s0-0.5)*u), [-4, -2, 0, 2, 4])
out["plancherel_CQ"] = {"int|Ef|^2": float(lhs), "(1/2pi)int|M|^2 on Re s=1/2": float(rhs),
                        "relerr": float(abs(lhs-rhs)/rhs),
                        "Mellin(Ef)(0.5+3i)": str(mp.nstr(mel, 12)), "formula": str(mp.nstr(M(s0), 12))}
# zeros of zeta on the line are zeros of M -- measure zero, invisible to the L^2 norm
out["M_at_first_crit_zero(test data)"] = float(abs(M(mp.zetazero(1))))

# 3. Weil-form arch weight vs Plancherel: w(tau) = (1/2pi)[Re psi(1/4+i tau/2) - log pi]
#    W restricted to a narrow packet at frequency T ~ w(T) * ||phi||^2 : unbounded relative to L^2.
w = lambda T: (mp.re(mp.digamma(mp.mpf(1)/4+1j*T/2)) - mp.log(mp.pi))/(2*mp.pi)
out["arch_weight"] = {str(T): float(w(T)) for T in [0, 5, 10, 14.13, 100, 1000, 1e4]}

# 4. Davenport-Heilbronn function: an exact GL1(A_Q) Tate integral with a NON-isotypic f_5.
#    chi mod 5 with chi(2)=i (odd).  kappa = (sqrt(10-2sqrt5)-2)/(sqrt5-1).
#    DH(s) = (1-i kappa)/2 L(s,chi) + (1+i kappa)/2 L(s,chibar)
#    = Z(f,|.|^s)/(arch) with f = x e^{-pi x^2} (x) f_5 (x) prod 1_Zp, f_5 = sum_a c(a) 1_{a+5Z_5}, c(a)=coeff of a.
chi = {1: 1, 2: 1j, 4: -1, 3: -1j, 0: 0}
kap = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
c = {a: (1-1j*kap)/2*chi[a] + (1+1j*kap)/2*mp.conj(chi[a]) for a in range(5)}
out["DH_coeffs_c(a)"] = {a: str(mp.nstr(c[a], 10)) for a in range(5)}
def DH(s):
    return mp.fsum(c[a]*mp.zeta(s, mp.mpf(a)/5) for a in range(1, 5)) * 5**(-s)
xi = lambda s: (mp.pi/5)**(-(s+1)/2)*mp.gamma((s+1)/2)*DH(s)
fe = {}
for s in [mp.mpc(0.3, 2.0), mp.mpc(0.8, 7.0)]:
    fe[str(complex(s))] = str(mp.nstr(xi(s)/xi(1-s), 15))
out["DH_FE_ratio xi(s)/xi(1-s)"] = fe
z = mp.findroot(DH, mp.mpc(0.8085, 85.6993))
out["DH_offline_zero(test data)"] = {"rho": str(mp.nstr(z, 12)), "|DH(rho)|": float(abs(DH(z)))}
# multiplicativity failure (no Euler product): c(2)c(3) vs c(6 mod 5 = 1)
out["DH_c2c3_vs_c1"] = [str(mp.nstr(c[2]*c[3], 10)), str(mp.nstr(c[1], 10))]

# 5. Epstein x^2+5y^2 : Z_Q = zeta L(chi_-20) + L(chi_-4) L(chi_5)  (13e decomposition); off-line zero test data
def Lmod(s, q, ch):
    return q**(-s)*mp.fsum(ch(a)*mp.zeta(s, mp.mpf(a)/q) for a in range(1, q) if ch(a) != 0)
def chi_m4(a): return [0, 1, 0, -1][a % 4]
def chi_5(a): return [0, 1, -1, -1, 1][a % 5]
def chi_m20(a):  # chi_-4 * chi_5
    return chi_m4(a)*chi_5(a)
ZQ = lambda s: mp.zeta(s)*Lmod(s, 20, chi_m20) + Lmod(s, 4, chi_m4)*Lmod(s, 5, chi_5)
# direct check of decomposition at s=2
direct = mp.nsum(lambda x, y: 0 if (x == 0 and y == 0) else (x**2+5*y**2)**(-2), [-mp.inf, mp.inf], [-mp.inf, mp.inf])
out["Epstein_decomp_s=2"] = {"direct": float(direct), "decomp": float(mp.re(ZQ(2)))}
zq = mp.findroot(ZQ, mp.mpc(0.93297, 15.66825))
out["Epstein_offline_zero(test data)"] = {"rho": str(mp.nstr(zq, 12)), "|Z_Q(rho)|": float(abs(ZQ(zq)))}

json.dump(out, open("r14a_adelic_ps.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
