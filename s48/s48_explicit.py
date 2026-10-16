"""S48 H48a: exact constant and explicit formula for the log-mean numerator
G(X) = int_1^X (psi(t)-t) t^{-3/2} dt  (= -int_1^X d_psi dt/t with d_psi=(x-psi)/sqrt x).
Mellin: int_1^oo G(X) X^{-s-1} dX = (1/s) [ -zeta'/zeta(s+1/2)/(s+1/2) - 1/(s-1/2) ]   (PROVED, Re s>1/2, then mero.)
Residues: s=0 -> c_G = 2 - 2 zeta'/zeta(1/2) = 2 - gamma - pi/2 - log(8 pi);
 s=1/2 removable; s=rho-1/2 -> -X^{rho-1/2}/(rho(rho-1/2)); s=-1/2 -> +2 log(2pi) X^{-1/2};
 s=-2n-1/2 -> -X^{-2n-1/2}/(2n(2n+1/2)).
Check G(1)=0 <=> sum_rho 1/(rho(rho-1/2)) = c_G + 2log2pi - (4 - pi/2 - 3log2)  (= -(2+gamma-log 4pi)).
"""
import mpmath as mp, numpy as np, json, math
mp.mp.dps = 40
z = mp.zeta(0.5); zp = mp.zeta(0.5, derivative=1)
cG = 2 - 2*zp/z
cG_closed = 2 - mp.euler - mp.pi/2 - mp.log(8*mp.pi)
T_triv = 4 - mp.pi/2 - 3*mp.log(2)   # sum_{n>=1} 1/(2n(2n+1/2))
Srho_pred = cG + 2*mp.log(2*mp.pi) - T_triv
Srho_closed = -(2 + mp.euler - mp.log(4*mp.pi))
print('c_G (zeta\'/zeta) =', mp.nstr(cG, 30)); print('c_G (closed)      =', mp.nstr(cG_closed, 30))
print('sum_rho 1/(rho(rho-1/2)) predicted =', mp.nstr(Srho_pred, 25), ' closed -(2+g-log4pi)=', mp.nstr(Srho_closed, 25))
# zeros
gam = []
with open('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt') as f:
    for line in f:
        a = line.split()
        if len(a) == 2: gam.append(a[1])
gm = [mp.mpf(g) for g in gam]
# zero-side check (uses RH-located zeros from the list, which are verified on the line in this height range)
s = mp.fsum([-2/(g*g + mp.mpf(1)/4) for g in gm])
T = gm[-1]
# tail bound: sum_{g>T} 2/(g^2+1/4) <= 2 int_T^oo dN(t)/t^2 ; N(t) = t/2pi log(t/2pi e) + 7/8 + Q, |Q| <= 0.112 log t + 0.278 loglog t + 2.51 (Trudgian 2014)
def Nmain(t): return t/(2*mp.pi)*mp.log(t/(2*mp.pi*mp.e)) + mp.mpf(7)/8
def Qb(t): return mp.mpf('0.112')*mp.log(t) + mp.mpf('0.278')*mp.log(mp.log(t)) + mp.mpf('2.51')
# int_T^oo dN/t^2 = -N(T)/T^2 + 2 int_T^oo N/t^3 ; upper: N<=Nmain+Qb, N(T)>=Nmain-Qb
tail_up = 2*(-(Nmain(T) - Qb(T))/T**2 + 2*mp.quad(lambda t: (Nmain(t) + Qb(t))/t**3, [T, mp.inf]))
tail_lo = 2*(-(Nmain(T) + Qb(T))/T**2 + 2*mp.quad(lambda t: (Nmain(t) - Qb(t))/t**3, [T, mp.inf]))
tail_main = 2*mp.quad(lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)/t**2, [T, mp.inf])
print('zero sum (100k) =', mp.nstr(s, 20), ' + tail in [', mp.nstr(-tail_up, 6), ',', mp.nstr(-tail_lo, 6), '] main', mp.nstr(-tail_main, 8))
print('residual vs closed (with main tail) =', mp.nstr(s - tail_main - Srho_closed, 6))
env = mp.fsum([2/(g*mp.sqrt(g*g + mp.mpf(1)/4)) for g in gm]) + tail_up  # RH envelope sup|sum X^{ig}/(ig rho)|
print('RH envelope sum 1/(|g||rho|) <=', mp.nstr(env, 12))
res = dict(c_G=mp.nstr(cG, 30), c_G_closed=mp.nstr(cG_closed, 30), E1_limit=mp.nstr(cG/2, 25),
           Srho_closed=mp.nstr(Srho_closed, 25), Srho_100k=mp.nstr(s, 20), tail_interval=[mp.nstr(-tail_up, 8), mp.nstr(-tail_lo, 8)],
           RH_envelope=mp.nstr(env, 12), nzeros=len(gm), T=mp.nstr(T, 12))
# ---- numpy explicit formula predictions (RH assumed, zeros from list) ----
g = np.array([float(x) for x in gam]); rho = 0.5 + 1j*g
lp2 = math.log(2*math.pi); cGf = float(cG)
I_triv_const = float(mp.quad(lambda t: 0.5*mp.log(1 - t**-2)*t**-0.5, [1, mp.inf]))
def pred(x, K=None):
    gg = g if K is None else g[:K]; rr = rho if K is None else rho[:K]
    L = math.log(x); ph = np.exp(1j*gg*L)
    triv = sum(x**(-2*n - 0.5)/(2*n*(2*n + 0.5)) for n in range(1, 30))
    G = cGf - 2*np.sum((ph/(1j*gg*rr)).real) + 2*lp2/math.sqrt(x) - triv
    dpsi = (2*np.sum((ph/rr).real) + (lp2 + 0.5*math.log(1 - x**-2))/math.sqrt(x))
    # Lin/x = int_1^x d_psi dt / x
    Iq = I_triv_const - float(mp.quad(lambda t: 0.5*mp.log(1 - t**-2)*t**-0.5, [x, mp.inf]))
    lin = (2*np.sum(((x**0.5*ph*x**0.5 - 1)/(rr*(rr + 0.5))).real) + 2*lp2*(math.sqrt(x) - 1) + Iq)/x
    return dict(G=G, d_psi=dpsi, E1=G/2 - dpsi, logmean=-G/L, linmean=lin)
sieve = json.load(open('sieve_1e9.json')) if __import__('os').path.exists('sieve_1e9.json') else {}
try: sieve.update({k: v for k, v in json.load(open('sieve_10000000000.json')).items()})
except Exception: pass
rows = []
xs = sorted(set([int(k) for k in sieve] + [10**11, 5*10**10]))
for x in xs:
    p = pred(float(x)); p1 = pred(float(x), 10000)
    row = dict(x=x, pred=p, pred_10k=p1, sieve=sieve.get(str(x)))
    rows.append(row)
    sv = sieve.get(str(x))
    print(f"x={x:.0e} G_pred={p['G']:.6f} G_sieve={sv['G'] if sv else float('nan'):.6f} | d_psi pred={p['d_psi']:.5f} sieve={sv['d_psi'] if sv else float('nan'):.5f} | logmean pred={p['logmean']:.5f} | linmean pred={p['linmean']:.5f} sieve={sv['linmean'] if sv else float('nan'):.5f} | E1 pred={p['E1']:.5f}")
# dyadic quantities at 1e11
for x in [1e11, 1e10, 1e9]:
    a = pred(x); b = pred(x/2)
    print(f"x={x:.0e}: E1(x)-E1(x/2)={a['E1']-b['E1']:.5f}; dyadic log-mean of d_psi (G(x/2)-G(x))/log2 = {(b['G']-a['G'])/math.log(2):.5f}; sqrt(log x)={math.sqrt(math.log(x)):.3f}")
    res[f'dyadic_{x:.0e}'] = dict(dE1=a['E1']-b['E1'], dyad_logmean=(b['G']-a['G'])/math.log(2))
res['rows'] = rows
json.dump(res, open('s48_explicit.json', 'w'), indent=1, default=float)
