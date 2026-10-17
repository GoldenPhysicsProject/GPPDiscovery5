"""Grok S57 H57b(ii): knockout margins delta_p and lambda_min(Q_zeta) at support log p, Ritz basis Legendre*(1-u^2)^EDGEK
(EDGEK=0 allowed: f in L^2[-L,L] need not vanish at the edge; the Weil form is finite there since g=f*f~ is Lipschitz).
modes:  lam SECTOR M L1,L2,...      -> generalized Ritz lambda_min/lambda_2 of Q_zeta (relative to the Gram matrix)
        bis SECTOR p M lo hi steps  -> bisect delta in [lo,hi] for first failure of Q_zeta + sum_{n=p^k} w_n T_n (Cholesky)
Failures are genuine (Ritz subspace of L^2[-L,L]) => upper ends are UPPER bounds on the true delta_p. NUMERICAL. RH OPEN."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's56')); sys.path.insert(1, os.path.join(HERE, '..', 's53'))
os.chdir(os.path.join(HERE, '..', 's56'))
from s53_odd_floor import build_all, tomp, weight, vonmangoldt
from s14_exact import G as Gpoly
import mpmath as mp
mp.mp.dps = 130; K = os.environ.get('EDGEK', '4')
def mats(sec, M, L):
    I = list(range(0 if sec == 'even' else 1, M, 2))
    A, P, T, ns = build_all(M, mp.nstr(mp.mpf(L), 40), 1024)
    W = tomp(A, I) + tomp(P, I); Tp = {}
    for n in ns:
        Tn = tomp(T[n], I); W -= weight(n) * Tn; Tp[n] = Tn
    Gm = mp.matrix(len(I), len(I))
    for a, i in enumerate(I):
        for b, j in enumerate(I):
            c0 = ((Gpoly(i, j) + Gpoly(j, i)) / 2).coeffs()[0]; Gm[a, b] = mp.mpf(L) * mp.mpf(int(c0.p)) / int(c0.q)
    return W, Tp, Gm, ns
def genlam(W, Gm):
    Li = mp.inverse(mp.cholesky(Gm)); S = Li * W * Li.T; S = (S + S.T) / 2
    ev = sorted(mp.eigsy(S, eigvals_only=True)); return ev[0], ev[1]
def psd(X):
    try: mp.cholesky(X); return True
    except Exception: return False
mode, sec = sys.argv[1], sys.argv[2]; t0 = time.time()
if mode == 'lam':
    M = int(sys.argv[3])
    for Ls in sys.argv[4].split(','):
        L = mp.mpf(eval(Ls, {'log': mp.log})); W, Tp, Gm, ns = mats(sec, M, L); l1, l2 = genlam(W, Gm)
        print(json.dumps(dict(mode='lam', sec=sec, M=M, EDGEK=K, L=Ls, Lnum=mp.nstr(L, 12), lam_min=mp.nstr(l1, 10), lam_2=mp.nstr(l2, 10),
                              wall_s=round(time.time() - t0))), flush=True)
else:
    p, M, lo, hi, steps = int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])
    L0 = mp.log(p) / 2
    def fails(d):
        W, Tp, Gm, ns = mats(sec, M, L0 + d)
        for n in ns:
            if vonmangoldt(n) == p: W += weight(n) * Tp[n]
        return not psd(W)
    note = ''
    if fails(lo): note = 'fails already at lo'
    elif not fails(hi): note = 'no failure up to hi'
    else:
        for _ in range(steps):
            mid = (lo + hi) / 2
            if fails(mid): hi = mid
            else: lo = mid
    print(json.dumps(dict(mode='bis', sec=sec, p=p, M=M, EDGEK=K, delta_in=[round(lo, 7), round(hi, 7)], note=note,
                          wall_s=round(time.time() - t0))), flush=True)
# Grok S57 (2026-10-17). RH OPEN. GRH OPEN.
