"""S45: Arb-CERTIFIED outer witness that F_θ prime-weight deformation W is indefinite.
Find Ritz λ_min < 0 eigenvector of W, round to rationals, enclose x^T W x < 0 in Arb.
Same style as s44_cert.py. CERTIFIED for the prime-weight deformation only (ζ arch).
"""
import sys, json
sys.path.insert(0, '.')
from s45_ftheta_scan import *

def cert(M, h, theta, prec):
    A, T, ns = build_all(M, h, prec)
    # Build Arb W for ftheta_def
    W = arb_mat(M, M)
    for i in range(M):
        for j in range(M):
            s = A[i, j]
            for n in ns:
                p = vonmangoldt(n)
                # w = 2 log(p)/sqrt(n) * (n^θ + n^{-θ})
                wn = (2 * arb(p).log() / arb(n).sqrt()
                      * (arb(n) ** arb(theta) + arb(n) ** (-arb(theta))))
                s = s - wn * T[n][i, j]
            W[i, j] = s
    row = {"M": M, "h": h, "theta": str(theta), "prime_powers": ns,
           "kind": "ftheta_prime_def", "arch": "zeta"}
    for nm, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        w = tomp(W, I)
        ev, V = mp.eigsy(w)
        idx = min(range(len(I)), key=lambda i: ev[i])
        lmin = ev[idx]
        y = V[:, idx]
        sc = max(abs(v) for v in y)
        y = [v / sc for v in y]
        # round to rationals with 10^70 denominator
        xs = [fmpq(int(mp.nint(v * mp.mpf(10) ** 70)), 10 ** 70) for v in y]
        xa = [arb(q.p) / q.q for q in xs]
        wv = sum(xa[a] * xa[b] * W[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
        # also enclose Rayleigh vs ||x||^2 mid (not Gram-normalized rigorously)
        nn = sum(xa[a] * xa[a] for a in range(len(I)))
        row[nm] = {
            "ritz_lam_min": mp.nstr(lmin, 14),
            "xTx_enclosure": wv.str(20),
            "xTx_negative_certified": bool(wv < 0),
            "x_norm2_mid": nn.mid().str(12, radius=False),
            "x_rationals_den": 10 ** 70,
            "x_num_preview": [str(q.p)[:20] + "..." for q in xs[:4]],
        }
    return row

if __name__ == "__main__":
    # argv: M hs_csv theta [dps] [prec]
    M = int(sys.argv[1])
    hs = [float(x) for x in sys.argv[2].split(',')]
    theta = mp.mpf(sys.argv[3])
    dps = int(sys.argv[4]) if len(sys.argv) > 4 else 110
    prec = int(sys.argv[5]) if len(sys.argv) > 5 else 768
    mp.mp.dps = dps
    for h in hs:
        print(json.dumps(cert(M, h, theta, prec)), flush=True)
