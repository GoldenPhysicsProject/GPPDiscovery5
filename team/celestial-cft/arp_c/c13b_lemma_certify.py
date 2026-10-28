"""Certified (Arb ball arithmetic, python-flint) constants for Lemma C.13.B and Conj C.13.H tables.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.
Omega(t) = Re digamma(1/4 + i t/2) - log(pi);  kappa(u) = e^{u/2}/(2 sinh u) = e^{3u/2}/(e^{2u}-1), u>0;
(Q kappa)(u) = -kappa''(u) + kappa(u)/4.  a = log 2, c2 = log(2)/sqrt(2).
Lemma C.13.B bracket(l) = delta*(1 - A l^3) - (c2 - Omega(0)) A l^3 - M l^3/pi^2,  A = (T/4 + T^3/3)/pi^3, T = 14,
delta = Omega(T) - c2, M = sup_{|v|<=lmax} |(Q kappa)(a+v)|.  bracket is decreasing in l on (0, lmax], so check at lmax.
"""
from flint import arb, acb, arb_series, ctx
ctx.prec = 128
pi = arb.pi(); a = arb(2).log(); c2 = a / arb(2).sqrt()
def Omega(t):
    return acb(arb(1)/4, arb(t)/2).digamma().real - pi.log()
def Qkappa_ball(lo, hi):
    u = arb((lo+hi)/2, (hi-lo)/2)            # ball covering [lo,hi]
    x = arb_series([u, 1], prec=3)
    k = (arb(3)/2*x).exp() / ((2*x).exp() - 1)
    k0, k2 = k[0], 2*k[2]
    return -k2 + k0/4
T = 14; lmax = arb("0.119")
O14, O0 = Omega(T), Omega(0)
delta = O14 - c2
A = (arb(T)/4 + arb(T)**3/3) / pi**3
N = 400; lo = a - lmax; hi = a + lmax; M = arb(0)
for i in range(N):
    b = Qkappa_ball(lo + (hi-lo)*i/N, lo + (hi-lo)*(i+1)/N)
    M = M.union(b) if i else b
Mup = abs(M).upper()
l3 = lmax**3
bracket = delta*(1 - A*l3) - (c2 - O0)*A*l3 - arb(Mup)*l3/pi**2
print("c2           =", c2)
print("Omega(0)     =", O0)
print("Omega(14)    =", O14)
print("delta        =", delta, " certified > 0:", delta > 0)
print("A            =", A, " A*l^3 < 1:", A*l3 < 1)
print("Qkappa range on [a-0.119, a+0.119] =", M, " => M_l <=", Mup)
print("bracket(0.119) =", bracket, " certified > 0:", bracket > 0)
# Omega levels used in C.13.H / C.14 (where Omega crosses c2 and 2c2): certified sign checks
for t, lvl, name in [(10.25, c2, "c2"), (10.30, c2, "c2"), (16.6, 2*c2, "2c2"), (16.8, 2*c2, "2c2")]:
    print(f"Omega({t}) - {name} =", Omega(arb(str(t))) - lvl)
# 2-adic multiplier extremes for C.14: Omega_2(t) = -2 log2 * sum_{m>=1} 2^{-m/2} cos(m a t)
r = 1/arb(2).sqrt()
print("min Omega_2 = -2 log2 * r/(1-r) =", -2*a*r/(1-r))
print("max Omega_2 =  2 log2 * r/(1+r) =", 2*a*r/(1+r))
