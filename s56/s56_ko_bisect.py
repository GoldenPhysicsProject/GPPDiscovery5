"""Grok S56 H56b: resolution test of the knockout margin delta_p = L_p* - (1/2)log p.
Ritz space = Legendre*(1-u^2)^EDGEK (EDGEK env, default 4 = S54/S55 basis), M functions (s14_exact copy in this folder).
For a sector (odd/even) and prime p, bisect L in [(1/2)log p, (1/2)log p + hi] for the first L where W_zeta + w_p T_p
(prime p removed, other terms kept) is NOT positive definite (mpmath Cholesky; sign is basis-independent by Sylvester).
A failure is a genuine witness (Ritz space is a subspace of L^2[-L,L]); so each reported upper end is an UPPER bound on the
true delta_p. If delta_p(M) keeps shrinking as M grows, the S55 'c/log p' law is a resolution artifact. NUMERICAL. RH OPEN.
usage: (cd s56 && python s56_ko_bisect.py SECTOR p M [hi] [steps])"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's53'))
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from s53_odd_floor import build_all, tomp, weight, psd
import mpmath as mp
sec, p, M = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
hi = float(sys.argv[4]) if len(sys.argv) > 4 else 0.04; steps = int(sys.argv[5]) if len(sys.argv) > 5 else 6
mp.mp.dps = 130; I = list(range(0 if sec == 'even' else 1, M, 2))
def fails(L):
    A, P, T, ns = build_all(M, '%.6f' % L, 1024)
    W = tomp(A, I) + tomp(P, I)
    for n in ns:
        if n != p: W -= weight(n) * tomp(T[n], I)
    return not psd(W)
L0 = 0.5 * float(mp.log(p)); lo, up = 0.0, hi; t0 = time.time()
if not fails(L0 + up):
    print(json.dumps(dict(sec=sec, p=p, M=M, EDGEK=os.environ.get('EDGEK', '4'), delta_gt=hi, note='no failure up to hi')), flush=True); sys.exit()
for _ in range(steps):
    mid = (lo + up) / 2
    if fails(L0 + mid): up = mid
    else: lo = mid
print(json.dumps(dict(sec=sec, p=p, M=M, EDGEK=os.environ.get('EDGEK', '4'), delta_in=[round(lo, 5), round(up, 5)],
                      dlogp_in=[round(lo * 2 * L0, 5), round(up * 2 * L0, 5)], wall_s=round(time.time() - t0))), flush=True)
