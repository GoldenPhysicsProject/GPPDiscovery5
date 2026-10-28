# N.15 NUMERICAL: is the S1-only template (R - c2 ||xi||^2_{H(E2)} - c3 ||xi||^2_{H(E3)} >= 0) viable on {2,3} windows?
# E_p = S-^{(p)} u S+^{(p)}, strips of length 2L - log p.  Ritz upper bounds of inf (.)/||xi||_H^2 per parity.
import sys, numpy as np
src = open(__file__.replace('n15_S1_obstruction.py', 'n15_P2_margins.py')).read().split('c2, c3, c4 =')[0]
exec(src)
def HE(x, h, L, ell):
    # Gram of int_{E}(u'v' + uv/4), E = [-L,-L+ell] u [L-ell,L]
    g, w = np.polynomial.legendre.leggauss(4); out = 0
    for lo, hi in [(-L, -L+ell), (L-ell, L)]:
        bp = np.unique(np.r_[lo, x[(x > lo) & (x < hi)], hi])
        y = np.concatenate([(a+b)/2 + (b-a)/2*g for a, b in zip(bp[:-1], bp[1:])])
        wt = np.concatenate([(b-a)/2*w for a, b in zip(bp[:-1], bp[1:])])
        v, s = hatval(x, h, y); out = out + (s*wt) @ s.T + 0.25*(v*wt) @ v.T
    return out
c2, c3 = np.log(2)/np.sqrt(2), np.log(3)/np.sqrt(3)
for lab, twoL, use3 in [("log2+0.028", np.log(2)+0.028, False), ("log2+0.1", np.log(2)+0.1, False), ("log3", np.log(3), False), ("log3+0.028", np.log(3)+0.028, True)]:
    L = twoL/2; N = 80
    x, h, GR, GH = grams(L, N); n = x.size
    A = GR - c2*HE(x, h, L, twoL - np.log(2)) - (c3*HE(x, h, L, twoL - np.log(3)) if use3 else 0)
    Rf = np.eye(n)[::-1]; out = []
    for sgn in (+1, -1):
        P = (np.eye(n) + sgn*Rf)/2; U, s, _ = np.linalg.svd(P); B = U[:, s > 0.5]
        Ci = np.linalg.inv(np.linalg.cholesky(B.T @ GH @ B)); out.append(np.linalg.eigvalsh(Ci @ (B.T @ A @ B) @ Ci.T)[0])
    print(f"2L={lab:11s} S1-only bound form: even <= {out[0]: .4f}   odd <= {out[1]: .4f}", flush=True)
