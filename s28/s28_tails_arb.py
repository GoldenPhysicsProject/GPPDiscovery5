"""Grok S28: outward-rounded (Arb) recomputation of the final one-stroke error budget, answering Codex 6204288080.
Inputs: a cert JSON from s28_onestroke_fast.py (or s27).  Replaces every float auxiliary bound by Arb:
 * Legendre-block tail (Codex's exact ratio bound): a_d^2 = 2L(2d+1)(LT)^{2d}/((2d+1)!!)^2 (from |j_d(x)| <= x^d/(2d+1)!!),
   r_d = a_{d+2}^2/a_d^2 = (LT)^4/[(2d+1)(2d+3)^2(2d+5)] decreasing in d, so sum_{d'>=d, same parity} a_{d'}^2 <= a_d^2/(1-r_d) if r_d < 1.
 * pole tail: p_d^2 <= 2L(2d+1)(L/2)^{2d} e^{L^2/4}/((2d+1)!!)^2 (i_d(x) <= x^d e^{x^2/2}/(2d+1)!!), same ratio with L/2 in place of LT.
 * pole head |p|^2 = sum_{d in block} p_d^2 with p_d = sqrt(2L(2d+1)) i_d(L/2) in Arb.
 * ||A - Atil|| <= K * max radius: the float radius is inflated by (1 + 2^-40) (float(mag) is round-to-nearest).
 * certified = min(lam0 - rE, lamD) - eps_B, lamD = beta - eps_CDD - [2 |p_tail|^2 if odd]; all as Arb lower/upper endpoints.
Status: computer-assisted (Arb), NOT Lean.  RH OPEN."""
import json, sys, os
from flint import arb, ctx, fmpq
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's27'))
from s27_onestroke_arb import sph_i
ctx.prec = 256

def dfact(d):  # (2d+1)!! exact
    r = arb(1)
    for j in range(1, 2 * d + 2, 2): r *= j
    return r

def tail(L, X, d, extra=arb(1)):
    a2 = 2 * L * (2 * d + 1) * X ** (2 * d) / dfact(d) ** 2 * extra
    r = X ** 4 / ((2 * d + 1) * (2 * d + 3) ** 2 * (2 * d + 5))
    assert r < 1, 'ratio test fails'
    return a2 / (1 - r)

def audit(fn):
    d = json.load(open(fn))
    L = arb(fmpq(*map(int, d['L'].split('/')))) if '/' in d['L'] else arb(d['L'])
    T = arb(d['T']); K = d['K']; beta = arb(d['beta']); M0 = arb(d['M0']) if 'M0' in d else None
    if M0 is None: raise SystemExit('need M0 in JSON')
    out = {}
    for par, key in ((0, 'even'), (1, 'odd')):
        s = d.get(key)
        if s is None or s.get('lam0_schur', s.get('lam0_cholesky')) is None: continue
        lam0 = arb(s.get('lam0_schur') or s.get('lam0_cholesky'))
        degs = [par + 2 * k for k in range(K)]; dK = par + 2 * K
        sumD = tail(L, L * T, dK)
        sumA = sum((2 * L * (2 * dd + 1) for dd in degs), arb(0))
        pA2 = sum(((2 * L * (2 * dd + 1)) * sph_i(dd, L / 2) ** 2 for dd in degs), arb(0))
        pD2 = tail(L, L / 2, dK, extra=(L * L / 4).exp())
        epsB = T / arb.pi() * M0 * (sumA * sumD).sqrt() + 2 * pA2.sqrt() * pD2.sqrt()
        epsC = T / arb.pi() * M0 * sumD
        lamD = beta - epsC - (2 * pD2 if par == 1 else 0)
        rE = arb(s['normE_bound']) * (1 + arb(2) ** -40)
        # quadrature: redo the Bernstein bound with a rational rho = 347/100 (< true rho 3.4868 for w=1/2, b=2/5;
        # its ellipse overhang 0.2198 < the swept 0.222), all in Arb; rE_new <= rE + K (Eq_new - Eq_old)
        w = arb(fmpq(*map(int, d['w'].split('/')))) if '/' in d['w'] else arb(d['w'])
        b = arb(fmpq(*map(int, d['b'].split('/')))) if '/' in d['b'] else arb(d['b'])
        assert d['w'] == '1/2' and d['b'] == '2/5'
        rho = arb(fmpq(347, 100)); npan = int(round(float(T / w)))
        Mint = arb(d['M_S']).upper() / arb.pi() * 2 * L * (2 * degs[-1] + 1) * (2 * L * b).exp()
        EqN = npan * w / 2 * arb(64) / 15 * Mint / ((rho ** 2 - 1) * rho ** (2 * (d['ng'] - 1)))
        rE = rE + K * (EqN - arb(s['quad_err_entry'])).upper() if EqN > arb(s['quad_err_entry']) else rE
        cert_lo = (lam0 - rE).lower(); lD = lamD.lower()
        cert = (arb(min(cert_lo, lD, key=float)) - epsB).lower()
        out[key] = dict(certified_lower_arb=float(cert), eps_B_upper=float(epsB.upper()), lamD_lower=float(lD),
                        tail_sumD_upper=float(sumD.upper()), pole_tail_upper=float(pD2.upper()), lam0=str(lam0), rE=float(rE.upper()),
                        certified_lower_float_pipeline=s.get('certified_lower'))
        print(fn, key, json.dumps(out[key]))
    return out

if __name__ == '__main__':
    res = {fn: audit(fn) for fn in sys.argv[1:]}
    json.dump(res, open('tails_arb_audit.json', 'w'), indent=1)
