"""Grok S54 H54b: knockout test of the S53 'odd floor is prime 2' finding.
S14/S44 exact Arb assembly (Legendre*(1-u^2)^4 basis, M functions; reuses s53_odd_floor.build_all), then per sector
lambda_min of W_zeta with selected prime-power terms -w_n T_n removed. Ritz (mpmath eigsy). NUMERICAL.
Ritz lambda_min >= true inf of the form on L^2[-L,L], so a NEGATIVE Ritz value is a genuine witness (up to numerics)
that the knocked-out form is indefinite; a positive one is not a proof of positivity.
usage: s54_knockout.py M h1,h2,... [dps] [prec]   RH OPEN."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's53'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s53_odd_floor import build_all, tomp, weight
import mpmath as mp
VARIANTS = {"full": [], "no2": [2], "no2pow": [2, 4, 8], "no3": [3], "no3pow": [3, 9], "no5": [5], "no7": [7], "no11": [11], "no13": [13], "noprimes": None}
def run(M, h, prec):
    t0 = time.time(); A, P, T, ns = build_all(M, h, prec)
    row = {"M": M, "h": h, "prime_powers": ns}
    for nm, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        a = tomp(A, I); p = tomp(P, I); Tm = {n: tomp(T[n], I) for n in ns}
        row[nm] = {}
        for vn, drop in VARIANTS.items():
            W = a + p
            for n in ns:
                if drop is None or n in drop: continue
                W -= weight(n) * Tm[n]
            ev = sorted(mp.eigsy(W, eigvals_only=True))
            row[nm][vn] = [mp.nstr(ev[0], 10), mp.nstr(ev[1], 6)]
    row["wall_s"] = round(time.time() - t0, 1); return row
if __name__ == "__main__":
    M = int(sys.argv[1]); hs = sys.argv[2].split(",")
    mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 130
    prec = int(sys.argv[4]) if len(sys.argv) > 4 else 1024
    for h in hs: print(json.dumps(run(M, h, prec)), flush=True)
