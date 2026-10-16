"""Grok S35: Arb re-evaluation of the float (libm) tail bookkeeping of a schurboost certificate JSON
(s33_onestroke_schurboost.py lines ~265-284: sumD, pD2, pA, eps_B, eps_CDD, lamD, certified_lower).
Same formulas, every step in Arb (python-flint) with explicit geometric tail bounds (no libm, no truncation by fudge factor):
  term(d) = 2L(2d+1) (LT)^(2d) / ((2d+1)!!)^2,  ratio term(d+2)/term(d) = (LT)^4 (2d+5)/((2d+1)(2d+3)^2(2d+5)^2) decreasing;
  sumD = sum_{d>=dK, step 2} term(d)  <= partial + term_next/(1-r_next);
  p_d <= sqrt(2L(2d+1)) (L/2)^d/(2d+1)!! * exp((L/2)^2/2)   (series bound for i_d);
  pD2 = sum_{d>=dK, step 2} p_d^2 (geometric tail), pA^2 = sum_{d in degs} p_d^2 <= sum over all d of parity (upper bound);
  eps_B = (T/pi) M0 sqrt(sumA sumD) + 2 pA sqrt(pD2);  eps_CDD = (T/pi) M0 sumD;  lamD = beta - eps_CDD - [odd] 2 pD2;
  certified_lower = min(lam0 - rE, lamD) - eps_B with rE = stored normE_bound * (1 + 2^-50) (K*max rad, float product).
Audits ARITHMETIC only (not the analytic derivation of the formulas). Usage: python s35_tail_audit.py cert.json   RH OPEN."""
import sys, json
from flint import arb, ctx
ctx.prec = 256
d = json.load(open(sys.argv[1]))
L = arb(d['L']); T = arb(d['T']); K = d['K']; M0 = arb(d['M0']).upper() if '+/-' in d['M0'] else arb(d['M0'])
M0 = arb(M0) if not isinstance(M0, arb) else M0
beta = arb(d['beta']); out = {}
def ldf(n):  # log (2n+1)!!
    return arb(2 * n + 2).lgamma() - n * arb(2).log() - arb(n + 1).lgamma()
def term(dd): return (2 * L * (2 * dd + 1)) * ((L * T).log() * (2 * dd) - 2 * ldf(dd)).exp()
def ratio(dd): return (L * T) ** 4 * (2 * dd + 5) / ((2 * dd + 1) * (2 * dd + 3) ** 2 * (2 * dd + 5) ** 2)
def p2(dd): return 2 * L * (2 * dd + 1) * (((L / 2).log() * dd - ldf(dd)) * 2 + (L / 2) ** 2).exp()
for par in ('even', 'odd'):
    e = d[par]; pr = 0 if par == 'even' else 1
    dK = e['first_omitted_degree']; degs = list(range(pr, dK, 2)); assert len(degs) == K
    s = arb(0); dd = dK
    while True:
        r = ratio(dd)
        if r < arb('0.5'): s += term(dd) / (1 - r); break
        s += term(dd); dd += 2
    sumD = s
    # pD2: p2 ratio also decreasing; (L/2)^4/((2d+3)^2(2d+5)^2)*(2d+5)/(2d+1) << 1
    q = (L / 2) ** 4 * (2 * dK + 5) / ((2 * dK + 1) * (2 * dK + 3) ** 2 * (2 * dK + 5) ** 2)
    pD2 = p2(dK) / (1 - q)
    pA2 = arb(0)
    for dd in range(pr, 60, 2): pA2 += p2(dd)
    pA2 += p2(60 + pr) * 2  # ratio at d>=60 is < 1e-8
    pA = pA2.sqrt()
    sumA = sum((2 * L * (2 * x + 1) for x in degs), arb(0))
    eps_B = T / arb.pi() * M0 * (sumA * sumD).sqrt() + 2 * pA * pD2.sqrt()
    eps_CDD = T / arb.pi() * M0 * sumD
    lamD = beta - eps_CDD - (2 * pD2 if pr == 1 else 0)
    lam0 = arb(e['lam0_schur']); rE = arb(e['normE_bound']) * (1 + arb(2) ** -50)
    a1 = lam0 - rE
    low = (a1 if a1.upper() < lamD.lower() else lamD) - eps_B
    out[par] = dict(sumD_upper=float(sumD.upper()), pD2_upper=float(pD2.upper()), pA_upper=float(pA.upper()),
                    eps_B_arb_upper=float(eps_B.upper()), eps_B_float=e['eps_B'], eps_CDD_arb_upper=float(eps_CDD.upper()),
                    eps_CDD_float=e['eps_CDD'], lamD_arb_lower=float(lamD.lower()), certified_lower_arb=low.lower().str(10, radius=False),
                    certified_lower_float=e['certified_lower'],
                    agrees=bool(low.lower() > 0 and float(eps_B.upper()) <= e['eps_B'] * (1 + 1e-6) + 1e-300))
print(json.dumps(out, indent=1))
