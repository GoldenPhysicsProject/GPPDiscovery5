# S49 (Grok): genus-0 Hadamard form of the Fisher margin in the shadow variable u=z^2,
#   -m'(u) = sum_{Im rho>0} 1/(u-(rho-1/2)^2)^2 ,
# (i) CONFIRM -m'_xi(0) = sum gamma^-4 (100k zeros + Trudgian tail bound) against the Arb value;
# (ii) table of -m'_xi(u) on u in [0,1e4]: Hadamard sum vs direct xi (NUMERICAL);
# (iii) per-pair sign rule and controls: c=5 two-mode fake (exact window), F_theta (theta=1/4),
#       DH (direct, completed), symmetric one-prime deformation at p=2.
import json, math, numpy as np, mpmath as mp
mp.mp.dps = 40
G = np.loadtxt('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt', usecols=1)
Gm = [mp.mpf(x) for x in open('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt').read().split()[1::2]]
assert len(Gm) == 100001 and len(G) == 100001
T = Gm[-1] + mp.mpf('1e-6')  # N(T) = 100001 exactly (list complete up to T)
NT = len(Gm)
def Nmain(t): return t/(2*mp.pi)*mp.log(t/(2*mp.pi*mp.e)) + mp.mpf(7)/8
def Qb(t): return mp.mpf('0.112')*mp.log(t) + mp.mpf('0.278')*mp.log(mp.log(t)) + mp.mpf('2.51')
# sum_{g>T} g^-4 = int_T^oo t^-4 dN = -N(T)/T^4 + 4 int_T^oo N(t) t^-5 dt   (Trudgian 2014: |N - Nmain| <= Qb, t >= e)
tail_up = -NT/T**4 + 4*mp.quad(lambda t: (Nmain(t) + Qb(t))/t**5, [T, mp.inf])
tail_lo = -NT/T**4 + 4*mp.quad(lambda t: (Nmain(t) - Qb(t))/t**5, [T, mp.inf])
tail_lo = max(tail_lo, mp.mpf(0))
S4 = mp.fsum([g**-4 for g in Gm])
cert = mp.mpf('3.717259928526968616486626247174057845365e-5')
print('sum_{100k} gamma^-4 =', mp.nstr(S4, 25))
print('tail in [', mp.nstr(tail_lo, 6), ',', mp.nstr(tail_up, 6), ']')
print('Arb -m\'(0) minus (S4 + tail) in [', mp.nstr(cert - S4 - tail_up, 6), ',', mp.nstr(cert - S4 - tail_lo, 6), ']')
res = dict(S4_100k=mp.nstr(S4, 30), tail=[mp.nstr(tail_lo, 8), mp.nstr(tail_up, 8)],
           residual=[mp.nstr(cert - S4 - tail_up, 8), mp.nstr(cert - S4 - tail_lo, 8)])
# (ii) table
def xi(s): return s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def mmp_direct(fun, U): return -mp.diff(lambda v: mp.log(fun(v)), U, 2)
G2 = G**2
def mmp_zeros(U):  # float, plus tail bound <= sum_{g>T} 1/(U+g^2)^2 <= tail_up
    return float(np.sum(1.0/(U + G2)**2))
tab = []
mp.mp.dps = 30
for U in [0, 1, 10, 100, 1000, 3000, 10000]:
    d = mmp_direct(lambda v: xi(mp.mpf(0.5) + mp.sqrt(v)), mp.mpf(U)) if U > 0 else cert
    z = mmp_zeros(U)
    tab.append(dict(u=U, direct=mp.nstr(d, 12), hadamard_100k=f'{z:.12e}', rel_diff=float((d - z)/d)))
    print(tab[-1])
res['table_u_0_1e4'] = tab
# (iii) sign rule constants
res['ratio_threshold_1_plus_sqrt2'] = mp.nstr(1 + mp.sqrt(2), 20)
H = mp.mpf(3)*10**12
res['u_floor_PlattTrudgian'] = mp.nstr(-(H**2 - H - mp.mpf(1)/4), 20)
g1 = Gm[0]
res['u_floor_no_computation_gamma1'] = mp.nstr(-(g1**2 - g1 - mp.mpf(1)/4), 20)
print('u floor (H=3e12):', res['u_floor_PlattTrudgian'], '  using only gamma>=gamma1:', res['u_floor_no_computation_gamma1'])
# c=5 two-mode fake
y = (5 + mp.sqrt(33))/4; dl = mp.acosh(y)
win = dl**2 - mp.pi**2 + 2*mp.pi*dl
print('c=5: delta=', mp.nstr(dl, 15), ' gamma=pi, gamma/delta=', mp.nstr(mp.pi/dl, 10), ' negative-term window u <', mp.nstr(win, 15))
def Fc(c): return lambda v: mp.cosh(2*mp.sqrt(v)) + c*mp.cosh(mp.sqrt(v))
def mmp_series(c, U):  # exact derivatives via Taylor in u: f(u)=sum (4^k + c)/(2k)! u^k
    a = [(mp.mpf(4)**k + c)/mp.factorial(2*k) for k in range(80)]
    f0 = mp.polyval(a[::-1], U); f1 = mp.polyval([k*a[k] for k in range(1, 80)][::-1], U)
    f2 = mp.polyval([k*(k - 1)*a[k] for k in range(2, 80)][::-1], U)
    return (f1**2 - f0*f2)/f0**2
fk = []
for U in [0, 0.5, 1, 2, 3, 3.1689, 5, 10, 50]:
    v = mmp_series(5, mp.mpf(U))
    u0 = dl**2 - mp.pi**2 + 2j*mp.pi*dl
    pair = 2*mp.re(1/(U - u0)**2)
    fk.append(dict(u=U, total=mp.nstr(v, 10), lowest_offaxis_pair=mp.nstr(pair, 10)))
print('c=5 fake -m\'(u):', fk)
cmin = []
for c in [1.5, 2, 5, 10, 30, 100, 1000]:
    vals = [mmp_series(c, mp.mpf(U)) for U in np.linspace(0, 40, 401)]
    cmin.append(dict(c=c, min_u_in_0_40=mp.nstr(min(vals), 8), at_0=mp.nstr(vals[0], 8)))
print('two-mode family min over u in [0,40]:', cmin)
res['c5_fake'] = dict(delta=mp.nstr(dl, 20), gamma_over_delta=mp.nstr(mp.pi/dl, 15), window_upper=mp.nstr(win, 20),
                      window_closed_form='acosh(y)^2 - pi^2 + 2 pi acosh(y), y=(5+sqrt33)/4', values=fk, family_min=cmin)
# F_theta, theta = 1/4: all zeros off line, delta = theta
th = 0.25
def mmp_Ftheta(U, th=th):
    w = U - (th + 1j*G)**2
    return float(np.sum(2*np.real(1/w**2)))
ft = [dict(u=U, val=f'{mmp_Ftheta(U):.6e}') for U in [0, 10, 1000]]
# direct check of F_theta at u=10
dft = mmp_direct(lambda v: xi(mp.mpf(0.5) + mp.sqrt(v) - th)*xi(mp.mpf(0.5) + mp.sqrt(v) + th), mp.mpf(10))
print('F_theta(1/4): ', ft, ' direct u=10:', mp.nstr(dft, 10))
# F_theta on u<0 near the lowest pair: own-term dip vs background
uu = np.linspace(-G[0]**2 - 30, -G[0]**2 + 30, 6001)
ftn = np.array([mmp_Ftheta(x) for x in uu])
print('F_theta min on window near -gamma1^2:', ftn.min(), 'at u=', uu[ftn.argmin()])
res['F_theta_quarter'] = dict(u_ge_0=ft, direct_u10=mp.nstr(dft, 12), min_near_minus_gamma1sq=float(ftn.min()), at=float(uu[ftn.argmin()]))
# DH completed function Lambda(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s)  (odd character mod 5)
mp.mp.dps = 30
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
chl = [0, 1, 1j, -1j, -1]; chb = [0, 1, -1j, 1j, -1]
def L(s, ch): return 5**(-s)*sum(ch[r]*mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5))
def fDH(s): return (1 - 1j*k)/2*L(s, chl) + (1 + 1j*k)/2*L(s, chb)
def LamDH(s): return (5/mp.pi)**((s + 1)/2)*mp.gamma((s + 1)/2)*fDH(s)
ev = abs(LamDH(mp.mpf(0.5) + 0.37 + 2j) - LamDH(mp.mpf(0.5) - 0.37 - 2j))/abs(LamDH(mp.mpf(0.5) + 0.37 + 2j))
print('DH evenness check |Lam(1/2+z)-Lam(1/2-z)|/|Lam| =', mp.nstr(ev, 5))
dh = []
for U in [0.001, 0.5, 1, 2, 5, 10, 25, 50, 100, 200, 400, 1000]:
    val = mmp_direct(lambda v: mp.re(LamDH(mp.mpf(0.5) + mp.sqrt(v))), mp.mpf(U))
    dh.append(dict(u=U, val=mp.nstr(val, 10), Lam=mp.nstr(mp.re(LamDH(mp.mpf(0.5) + mp.sqrt(U))), 8)))
print('DH -m\'(u):', dh)
# dense sign scan of DH on [0, 400]
dense = [float(mmp_direct(lambda v: mp.re(LamDH(mp.mpf(0.5) + mp.sqrt(v))), mp.mpf(U))) for U in np.linspace(0.01, 400, 400)]
print('DH dense min on (0,400]:', min(dense))
res['DH'] = dict(evenness_rel=mp.nstr(ev, 5), values=dh, dense_min_0_400=min(dense),
                 offline_zero=dict(delta=0.308517, gamma=85.699348, ratio=85.699348/0.308517))
# symmetric one-prime deformation at p=2: xi(s) [(1-2^{-s})(1-2^{s-1})]^m = xi(s) (3/2 - sqrt2 cosh(z log2))^m
for mdeg in (1, 2):
    vals = []
    for U in [0, 0.1, 0.25 + 1e-3, 1, 10, 100]:
        fun = lambda v: xi(mp.mpf(0.5) + mp.sqrt(v))*(mp.mpf(3)/2 - mp.sqrt(2)*mp.cosh(mp.sqrt(v)*mp.log(2)))**mdeg
        if U == 0:
            U = mp.mpf('1e-6')
        val = -mp.diff(lambda v: mp.log(abs(fun(v))), mp.mpf(U), 2)
        vals.append(dict(u=float(U), val=mp.nstr(val, 10)))
    res[f'one_prime_sym_m{mdeg}'] = vals
    print('one-prime sym m=', mdeg, vals)
json.dump(res, open('s49_hadamard.json', 'w'), indent=1, default=str)
