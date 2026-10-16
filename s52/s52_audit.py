# Grok S52 H52a: audit of Codex 6225514212 items (1),(2),(4). Arb (python-flint) where CERTIFIED.
from flint import arb, acb, ctx
from fractions import Fraction as Fr
import mpmath as mp
ctx.prec = 300
def xi(s):  # completed xi(s)=1/2 s(s-1) pi^{-s/2} Gamma(s/2) zeta(s), acb
    s = acb(s); return s*(s-1)/2 * acb.pi()**(-s/2) * (s/2).gamma() * s.zeta()
x12 = xi(arb(1)/2).real; x52 = xi(arb(5)/2).real
print('xi(1/2) =', x12.str(25), ' > 1071/3200=0.3346875 :', x12 > arb(1071)/3200)
print('xi(5/2) =', x52.str(25), ' < 85/128=0.6640625 :', x52 < arb(85)/128)
d = 2*x12 - x52
print('2f(0)-f(4) = 2xi(1/2)-xi(5/2) =', d.str(25), ' > 17/3200=0.0053125 :', d > arb(17)/3200)
# f'(4) = F'(2)/(2*2) with F(z)=xi(1/2+z); xi'(5/2) via mpmath high precision (cross-check)
mp.mp.dps = 60
xim = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
fp4 = mp.diff(xim, mp.mpf(5)/2)/4
M_true = fp4/(2*xim(mp.mpf(1)/2)-xim(mp.mpf(5)/2))
print("f'(4) = xi'(5/2)/4 =", mp.nstr(fp4, 20))
print('M(true dens) = f\'(4)/(2f0-f4) =', mp.nstr(M_true, 15), '; with Codex lower denominator 17/3200: M <=', mp.nstr(fp4/(mp.mpf(17)/3200), 15))
print('trace bounds 16M/9, 4M/9 (true-denominator M):', mp.nstr(16*M_true/9, 10), mp.nstr(4*M_true/9, 10))
# h0 = f'(0)/f(0) = xi''(1/2)/(2 xi(1/2)) = sum_{gamma>0} gamma^{-2}  (zero-free)
h0 = mp.diff(xim, mp.mpf(1)/2, 2)/(2*xim(mp.mpf(1)/2))
g = [mp.mpf(l.split()[0]) for l in open('../s51/zeros_arb_3000.txt')]
mp.mp.dps = 40
def hn(n, G=g):
    return mp.fsum(gg**(-2*n-2) for gg in G)
print('h0 zero-free =', mp.nstr(h0, 20), ' ; sum over 3000 zeros =', mp.nstr(hn(0), 20), '(tail ~ log/H)')
# Actual trace norms (H0,H1 PSD numerically): tr H0 = sum_i h_{2i}, tr H1 = sum_i h_{2i+1}
trH0 = h0 + mp.fsum(hn(2*i) for i in range(1, 30)); trH1 = mp.fsum(hn(2*i+1) for i in range(30))
print('tr H0 ~', mp.nstr(trH0, 12), ' tr H1 ~', mp.nstr(trH1, 12), '(NUMERICAL; = trace norms if PSD)')
print('Hausdorff support under RH: max x = 1/gamma1^2 =', mp.nstr(1/g[0]**2, 12), ' (vs Codex 1/4)')
# |h_n| decay: radius gamma1^2 unconditionally (no zeros with 0<gamma<14.13)
for n in [1, 5, 10, 20]:
    print(' n=%d h_n=%s  Codex bound M 4^-n=%s  gamma1^{-2n-2}=%s' % (n, mp.nstr(hn(n), 6), mp.nstr(M_true*mp.mpf(4)**(-n), 6), mp.nstr(g[0]**(-2*n-2), 6)))
# (4) F_theta: positive pushforward measure check, and K_theta sign at tau=pi/(2 theta gamma1)
import numpy as np
from scipy.integrate import quad
def Phi(t):
    t = abs(t); s = 0.0
    for n in range(1, 8):
        s += (2*np.pi**2*n**4*np.exp(4.5*t) - 3*np.pi*n**2*np.exp(2.5*t))*np.exp(-np.pi*n**2*np.exp(2*t))
    return s
F = lambda z: quad(lambda t: Phi(t)*np.cosh(z*t), -6, 6, limit=200)[0]
print('F(0)=int Phi =', F(0.0), ' xi(1/2)=', float(x12.mid()))
th, z = 0.2, 0.37
gw = lambda w: quad(lambda t: Phi(t)*Phi(w-t)*np.cosh(th*(2*t-w)), -6, 6, limit=200)[0]
lhs = F(z-th)*F(z+th)
rhs = quad(lambda w: gw(w)*np.cosh(z*w), -10, 10, limit=200)[0]
print('F_theta(z)=F(z-th)F(z+th)=%.12g ; int cosh(zw) dmu_theta=%.12g ; rel diff %.2e' % (lhs, rhs, abs(lhs-rhs)/abs(lhs)))
print('min density of mu_theta on grid (>=0 expected):', min(gw(w) for w in np.linspace(-8, 8, 41)))
mp.mp.dps = 80
th = mp.mpf('0.2'); tau = mp.pi/(2*th*g[0])
Kth = mp.fsum(2*mp.exp((th**2-gg**2)*tau)*mp.cos(2*th*gg*tau) for gg in g[:200])
print('tau=pi/(2 theta gamma1)=', mp.nstr(tau, 15), ' K_theta(tau)=', mp.nstr(Kth, 8), ' first term=', mp.nstr(2*mp.exp((th**2-g[0]**2)*tau)*mp.cos(2*th*g[0]*tau), 8))
print('K_theta sign change predicted at pi/(4 theta gamma1)=', mp.nstr(mp.pi/(4*th*g[0]), 6))
