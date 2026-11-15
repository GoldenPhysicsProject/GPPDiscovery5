"""(sigma-rig) Rigorous lower envelope of sigma_x(t), t>=0, in Arb (python-flint).
sigma_x(t) = 2C0 + 6 Cin(t/2) + 6R(t) + B(t) + P(t)
  Cin(x) = gamma + log x - Ci(x) = int_0^x (1-cos v)/v dv
  R(t)  = int_0^{1/2} f(u)(1-cos ut) du,  f(u) = h(u) - 1/u >= 0  (h = e^{u/2}/sinh u)   [PROVED: u e^{u/2} >= sinh u]
  B(t)  = 4 int_0^{1/2} kappa(s) cos(st) ds
  P(t)  = 4 sum c_n cos(mu_n t)
Grid part: sigma at grid pts in balls (R dropped on [0,eps] only, R>=0), cell bound min(endpoints) - L2*d^2/8.
Tail part (t>=T_big): sigma >= L(t) = 2C0 + 6(gamma+log(t/2)) - 6(2/t+4/t^2) - (4/t)(|kappa(1/2)|+(1/2)sup|kappa'|) - 4 sum|c_n|  (R>=0 dropped),
  increasing in t."""
import json, sys
from flint import arb, acb, ctx
ctx.prec = 80
half = arb(1)/2
pi = arb.pi(); gam = arb.const_euler()
def h(u): return (u/2).exp()/u.sinh()
def hp(u): return (u/2).exp()*(u.sinh()/2 - u.cosh())/u.sinh()**2
def kap(s): return h(1+s) - h(2+s)/4 + h(1-s) - h(2-s)/4
def kapp(s): return hp(1+s) - hp(2+s)/4 - hp(1-s) + hp(2-s)/4
# constants
w0 = arb(0.25).digamma() - pi.log()
U = arb(200)
H5a = acb.integral(lambda u,an: (u/2).exp()/u.sinh(), half, U).real
tail_hi = 4*(-U/2).exp()/(1-(-2*U).exp())
H5 = H5a + arb(0, 1)*tail_hi/2 + tail_hi/2   # in [H5a, H5a+tail_hi]
C0 = 3*w0 + 3*H5
NSn = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
mus = []
for n,p in NSn:
    lam = arb(n).log(); ln = arb(p).log()/arb(n).sqrt()
    if abs(float(lam)-1) < 0.5: mus.append((abs(lam-1), 2*ln, n))
    if abs(float(lam)-2) < 0.5: mus.append((abs(lam-2), -ln/2, n))
S4 = 4*sum(abs(c) for _,c,_ in mus)
# kappa > 0 and sup |kappa'| on [0,1/2] by interval subdivision
NS = 400
kmin = None; kpmax = arb(0)
for i in range(NS):
    s = arb((2*i+1), 1)/(4*NS)  # ball [i/(2NS),(i+1)/(2NS)]
    kv = kap(s); kpv = abs(kapp(s))
    kmin = kv.lower() if kmin is None else min(kmin, kv.lower())
    kpmax = max(kpmax, kpv.upper())
kmin = float(kmin)
# second-derivative bound L2 >= sup|sigma''| = |6 int h u^2 cos - 4 int kappa s^2 cos - 4 sum c mu^2 cos|
I_hu2 = acb.integral(lambda u,an: u*u*(u/2).exp()/u.sinh() if u.real.contains(0)==False else u*(u/2).exp()*(u/u.sinh()), arb(10)**-12, half).real
I_hu2 = I_hu2 + arb(10)**-12   # [0,1e-12] piece: u^2 h <= 2u there
I_ks2 = acb.integral(lambda s,an: s*s*kap(s), 0, half).real   # kappa>0 so this is int s^2|kappa|
L2 = (6*I_hu2 + 4*I_ks2 + 4*sum(abs(c)*m*m for m,c,_ in mus)).upper()
eps = arb(10)**-8
def sigma(t):
    """ball enclosure of a LOWER-valid quantity: exact sigma minus the nonneg piece 6*int_0^eps f(1-cos)"""
    t = arb(t)
    if t == 0:
        cin = arb(0); Rv = arb(0)
    else:
        x = t/2; cin = gam + x.log() - x.ci()
        Rv = acb.integral(lambda u,an: (u*(u/2).exp()-u.sinh())/(u*u.sinh())*(1-(u*t).cos()), acb(eps), half).real
    Bv = 4*acb.integral(lambda s,an: kap(s)*(s*t).cos(), 0, half).real
    Pv = 4*sum(c*(m*t).cos() for m,c,_ in mus)
    return 2*C0 + 6*cin + 6*Rv + Bv + Pv
def Ltail(T):
    T = arb(T)
    return 2*C0 + 6*(gam + (T/2).log()) - 6*(2/T + 4/T**2) - (4/T)*(abs(kap(half)) + half*kpmax) - S4
if __name__ == '__main__':
    Tbig = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    from fractions import Fraction as F
    pts = [F(i,20) for i in range(0, 40*20)] + [F(i,2) for i in range(80, 2*Tbig+1)]
    vals = [sigma(arb(p.numerator)/p.denominator) for p in pts]
    cells = []
    for i in range(len(pts)-1):
        d = pts[i+1]-pts[i]
        lo = min(float(vals[i].lower()), float(vals[i+1].lower())) - float(L2)*float(d)**2/8 - 1e-12
        cells.append((float(pts[i]), float(pts[i+1]), lo))
    out = dict(prec=80, C0=C0.str(15), w0=w0.str(15), H5=H5.str(15), sum4abs_c=S4.str(15), kappa_min_lower=kmin,
               sup_abs_kappa_prime_upper=float(kpmax), kappa_half=kap(half).str(12), L2_upper=float(L2), Tbig=Tbig,
               Ltail_Tbig=Ltail(Tbig).str(12), sigma0=vals[0].str(12))
    # S1: last cell with nonpositive lower bound
    neg = [c for c in cells if c[2] <= 0]
    T0 = max(c[1] for c in neg)
    out['T0_sigma_pos_beyond'] = T0
    out['min_lower_on_0_T0'] = min(c[2] for c in cells if c[1] <= T0)
    for name, T in [('omega13', 13*4*float(pi)), ('omega73', 73*4*float(pi))]:
        sub = [c for c in cells if c[1] > T]
        out['min_lower_t_ge_'+name] = min(min(c[2] for c in sub), float(Ltail(Tbig).lower()))
        out[name] = T
    out['min_lower_on_T0_omega13'] = min(c[2] for c in cells if c[0] >= T0 and c[0] < 13*4*float(pi))
    out['grid_points'] = len(pts)
    print(json.dumps(out, indent=1)); json.dump(out, open('sigrig_results.json','w'), indent=1)
