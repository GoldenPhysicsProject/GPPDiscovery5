"""S84: zero-side (spectral) decomposition of fixed rational Ritz vectors of the window Weil form.
Even basis c_m(u)=cos(a_m u), a_m=(2m+1)pi/(2L); odd s_m(u)=sin(b_m u), b_m=m pi/L (m>=1), on [-L,L] (s62/s75 dhritz convention).
Closed forms (derived here, checked numerically by quadrature in selftest):
  even: fhat(t)=int f e^{itu}du = 2 cos(tL) sum_m (-1)^m v_m a_m/(a_m^2-t^2)
  odd : |fhat(t)| = 2|sin(tL) sum_m (-1)^m v_m b_m/(b_m^2-t^2)|
Explicit formula (S14 normalisation, pole term = 2 fhat(i/2)^2): Q(f) = sum over ALL zeros rho=1/2+i gamma of fhat(gamma) conj(fhat(conj gamma))
 = 2 sum_{gamma>0} |fhat(gamma)|^2 for zeta under RH (the zeros used here are on the line, verified tables).
NUMERICAL (mpmath); the Rayleigh quotients compared against are the CERTIFIED Arb values of the cosvec files."""
import json, sys, mpmath as mp
from fractions import Fraction
import os; mp.mp.dps = int(os.environ.get("DPS","120"))
Z3 = [mp.mpf(l.split()[0]) for l in open('/workspace/GPPDiscovery5/s51/zeros_arb_3000.txt')]
def load_tail(nmax):
    out = []
    for l in open('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt'):
        i, g = l.split()
        if 3000 < int(i) <= nmax: out.append(mp.mpf(g))
    return out
def fhat_factory(vec, L, sector):
    Lm = mp.mpf(Fraction(L).numerator) / Fraction(L).denominator
    v = [mp.mpf(Fraction(q).numerator) / Fraction(q).denominator for q in vec]; K = len(v)
    if sector == 'even': fr = [(2*m+1)*mp.pi/(2*Lm) for m in range(K)]
    else: fr = [(m+1)*mp.pi/Lm for m in range(K)]
    c = [(-1)**m * v[m] * fr[m] for m in range(K)]
    if sector == 'even':
        f = lambda t: 2*mp.cos(t*Lm)*mp.fsum(c[m]/(fr[m]**2 - t**2) for m in range(K))
    else:
        f = lambda t: 2*mp.sin(t*Lm)*mp.fsum(c[m]/(fr[m]**2 - t**2) for m in range(K))
    norm2 = Lm*mp.fsum(x**2 for x in v)
    moms = [mp.fsum(c[m]*fr[m]**(2*j) for m in range(K)) for j in range(4)]
    return f, norm2, Lm, v, fr, moms
def analyse(fn, sector, ntail=20000, tag=''):
    d = json.load(open(fn)); vec = d['vec']; L = d['L'].replace('_', '/')
    f, norm2, Lm, v, fr, moms = fhat_factory(vec, L, sector)
    rq = mp.mpf(d['RQ'].strip('[').split('+/-')[0])
    target = rq*norm2
    terms = [2*f(g)**2 for g in Z3]
    tz = load_tail(3000 + ntail); terms += [2*f(g)**2 for g in tz]
    gam = Z3 + tz
    S = mp.fsum(terms)
    T = gam[-1]; M1 = moms[0]
    # tail beyond T: |fhat|^2 ~ 4 cos^2(tL) M1^2/t^4 -> avg 2 M1^2/t^4, density log(t/2pi)/(2pi), times 2 (both signs)
    tail = 2*mp.quad(lambda t: 2*M1**2/t**4*mp.log(t/(2*mp.pi))/(2*mp.pi), [T, mp.inf])
    order = sorted(range(len(terms)), key=lambda i: -terms[i])
    cum = lambda X: mp.fsum(terms[i] for i in range(len(terms)) if gam[i] <= X)
    Tstar = 2*mp.pi*mp.e**(2*Lm)   # zero density log(T/2pi)/2pi = L/pi (Paley-Wiener sampling density of type L)
    res = {"file": fn.split('/')[-1], "tag": tag, "L": L, "sector": sector, "K": len(v),
           "RQ_certified": d['RQ'][:40], "Q_geom=RQ*|f|^2": mp.nstr(target, 12),
           "Q_spec_partial(N=%d zeros, gamma<=%.0f)" % (len(gam), float(T)): mp.nstr(S, 12),
           "tail_est_beyond_T": mp.nstr(tail, 4), "spec/geom": mp.nstr((S+tail)/target, 10) if target > 0 else None,
           "M1=sum(-1)^m v_m a_m (1/t^2 coefficient)/|v|": mp.nstr(M1/mp.sqrt(norm2/Lm), 4),
           "T*=2pi e^{2L} (density=L/pi)": mp.nstr(Tstar, 5),
           "frac_gamma1": mp.nstr(terms[0]/S, 6), "frac_gamma2": mp.nstr(terms[1]/S, 6), "frac_gamma3": mp.nstr(terms[2]/S, 6),
           "frac_gamma<=50": mp.nstr(cum(50)/S, 6), "frac_gamma<=T*": mp.nstr(cum(Tstar)/S, 6),
           "frac_gamma<=2T*": mp.nstr(cum(2*Tstar)/S, 6), "frac_gamma<=max_freq": mp.nstr(cum(fr[-1])/S, 6),
           "max_basis_freq": mp.nstr(fr[-1], 5),
           "top5_zeros(index,gamma,frac)": [(i+1, mp.nstr(gam[i], 8), mp.nstr(terms[i]/S, 4)) for i in order[:5]],
           "|fhat(gamma_k)|^2/|f|^2 k=1..5": [mp.nstr(terms[k]/2/norm2, 4) for k in range(5)]}
    if sector == 'even':
        pole = 2*f(mp.mpc(0, 0.5))**2
        res["pole_term_2fhat(i/2)^2/|f|^2"] = mp.nstr(mp.re(pole)/norm2, 6)
    return res, terms, gam
if __name__ == '__main__':
    fn, sec = sys.argv[1], sys.argv[2]; nt = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
    r, _, _ = analyse(fn, sec, nt); print(json.dumps(r, indent=1))
