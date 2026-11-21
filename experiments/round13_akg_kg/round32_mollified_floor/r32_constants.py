"""Round 32: certified (Arb ball arithmetic, python-flint) normalisation constants of the Moebius floor
  m1(L) = sum_{d<=e^L} mu_M(d) d^{-1/2}(L - log d),  Laplace transform 1/(b_M(s+1/2) s^2),
  zero-free main part (residue at s=0):  L / b(1/2) - b'(1/2)/b(1/2)^2,
and for ZETA also the m2 and m^(1) main parts and c1 = (digamma(1/4) - log pi)/2.  Also checks the DH / EPS
closed forms against their coefficient definitions at s = 2. No zeros anywhere."""
import json, math
from flint import acb, arb, acb_series, ctx
ctx.prec = 256
KAP = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
kap = (arb(10) - 2*arb(5).sqrt()).sqrt() - 2; kap = kap/(arb(5).sqrt() - 1)

def ser(s0, n=3): return acb_series([acb(s0), 1], prec=n)
def hz(s, a): return s.zeta(acb(a))            # Hurwitz zeta series in s
def Lchar(s, q, chi):                            # L(s, chi) = q^-s sum chi(a) zeta(s, a/q)
    tot = None
    for a in range(1, q):
        if chi[a % q] == 0: continue
        t = hz(s, arb(a)/q)*chi[a % q]
        tot = t if tot is None else tot + t
    return tot*(acb(q)**(-s))  if False else tot*((-s)*acb(q).log()).exp()

chi_m4 = {0: 0, 1: 1, 2: 0, 3: -1}
chi_5 = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
def chi_m20(a):
    a %= 20
    if math.gcd(a, 20) > 1: return 0
    return chi_m4[a % 4]*chi_5[a % 5]
chi_m20d = {a: chi_m20(a) for a in range(20)}

def b_series(name, s0, th=None):
    s = ser(s0)
    if name == 'ZETA': return s.zeta()
    if name == 'FTH': return (s + th).zeta()*(s - th).zeta()
    if name == 'DH':
        f = hz(s, arb(1)/5) + hz(s, arb(2)/5)*kap - hz(s, arb(3)/5)*kap - hz(s, arb(4)/5)
        return f*((-s)*acb(5).log()).exp()
    if name == 'EPS':
        return (s.zeta()*Lchar(s, 20, chi_m20d) + Lchar(s, 4, chi_m4)*Lchar(s, 5, chi_5))/2

def coef(name, n, th=None):
    if name == 'DH': return [0, 1, KAP, -KAP, -1][n % 5]
    if name == 'EPS':
        M = int(n**0.5) + 2
        return sum(1 for x in range(-M, M + 1) for y in range(-M, M + 1) if x*x + 5*y*y == n)/2

out = {}
for name, th in [('ZETA', None), ('FTH', 0.0), ('FTH', 0.1), ('FTH', 0.25), ('DH', None), ('EPS', None)]:
    f = b_series(name, 0.5, th)
    b0, b1, b2 = f[0], f[1], 2*f[2]
    key = name if th is None else f'FTH{th}'
    out[key] = dict(b_half=str(b0.real), db_half=str(b1.real), d2b_half=str(b2.real),
                    slope=float((1/b0).real.mid()), intercept=float((-b1/b0**2).real.mid()),
                    slope_ball=str((1/b0).real), intercept_ball=str((-b1/b0**2).real))
    if name in ('DH', 'EPS'):     # closed form vs coefficients at s = 2
        direct = sum(coef(name, n)/n**2 for n in range(1, 4000))
        out[key]['check_s2'] = [float(b_series(name, 2.0)[0].real.mid()), direct]
z0, z1, z2 = (lambda f: (f[0], f[1], 2*f[2]))(b_series('ZETA', 0.5))
inv2 = 2*z1**2/z0**3 - z2/z0**2                  # (1/zeta)''(1/2)
c1 = (acb(0.25).digamma() - acb(math.pi).log())/2 if False else (arb(0.25).digamma() - arb.pi().log())/2
out['ZETA'].update(c1=str(c1), zeta_log_deriv_half=str((z1/z0).real), check_c1_plus_zlogderiv=str((c1 + z1/z0).real),
                   inv_zeta_dd_half=str(inv2.real),
                   m2_main='(1/2)[L^2/z0 - 2 L z1/z0^2 + (1/zeta)\'\'(1/2)]')
out['ZETA'].update(z0=float(z0.real.mid()), z1=float(z1.real.mid()), z2=float(z2.real.mid()), inv2=float(inv2.real.mid()),
                   c1f=float(c1.mid()))
json.dump(out, open('r32_constants.json', 'w'), indent=1)
for k, v in out.items(): print(k, {a: v[a] for a in ('b_half', 'db_half', 'slope', 'intercept')}, v.get('check_s2', ''))
print('c1 + zeta\'/zeta(1/2) =', out['ZETA']['check_c1_plus_zlogderiv'])
