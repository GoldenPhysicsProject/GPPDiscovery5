"""S84 (2): channel analysis of the CERTIFIED DH negative witnesses (s75 cosvecs). NUMERICAL (mpmath).
 - zeta zero-side decomposition of Q_zeta(f_DH) (certified value at 137/80: 1.50064301068282e-8 * |f|^2)
 - pole term 2 fhat(i/2)^2 (the term zeta has and DH lacks)
 - fhat profile; overlaps with zeta's near-null vector at the same L (same basis -> coefficient overlap)."""
import json, sys, mpmath as mp
from spec_decomp import fhat_factory, Z3, load_tail
mp.mp.dps = 60
def zeta_side(fn, ntail=5000):
    d = json.load(open(fn)); L = d['L'].replace('_','/')
    f, n2, Lm, v, fr, moms = fhat_factory(d['vec'], L, 'even')
    gam = Z3 + load_tail(3000+ntail); terms = [2*f(g)**2 for g in gam]; S = mp.fsum(terms)
    M1 = moms[0]; T = gam[-1]
    tail = 2*mp.quad(lambda t: 2*M1**2/t**4*mp.log(t/(2*mp.pi))/(2*mp.pi), [T, mp.inf])
    order = sorted(range(len(terms)), key=lambda i: -terms[i])
    pole = 2*f(mp.mpc(0,0.5))**2
    prof = {str(t): mp.nstr(abs(f(t))/mp.sqrt(n2), 3) for t in [0, 5, 10, 14.1347, 20, 30, 50, 85.7, 114.2, 150, 200]}
    return {"file": fn.split('/')[-1], "L": L, "K": len(v), "RQ_DH_certified": d['RQ'][:40],
            "Q_zeta_spec/|f|^2 (N zeros + tail)": mp.nstr((S+tail)/n2, 10),
            "frac_gamma1,2,3": [mp.nstr(terms[k]/S, 5) for k in range(3)],
            "frac_gamma<=50": mp.nstr(mp.fsum(terms[i] for i in range(len(gam)) if gam[i] <= 50)/S, 6),
            "top5": [(i+1, mp.nstr(gam[i], 8), mp.nstr(terms[i]/S, 4)) for i in order[:5]],
            "pole_term/|f|^2": mp.nstr(mp.re(pole)/n2, 8),
            "|fhat(t)|/|f|": prof}, (f, n2, v, Lm)
if __name__ == '__main__':
    out = []
    for fn in ["../s75/cosvec_DH_even_L137_80_K129.json", "../s75/cosvec_DH_even_L3421_2000_K300.json"]:
        r, _ = zeta_side(fn); out.append(r); print(json.dumps(r, indent=1))
    json.dump(out, open("dh_channels_results.json", "w"), indent=1)
