"""Grok S68 H68b (v3): lambda_min(T) and its eigenvector for a dumped (indefinite) T, NUMERICAL at 640 bits.
 1. Bisection on l0 (log scale, negative) with the s36 approximate blocked Cholesky of T - l0 I (success => T - l0 I > 0
    numerically; failure pivot recorded): brackets lambda_min to relative width 1e-3 (numerically; the s36 rigorous residual
    makes the 'success' side CERTIFIED for the dump if wanted).
 2. Shifted inverse iteration at sigma = success end of the bracket (T - sigma I > 0, so the nearest eigenvalue is lambda_min,
    contraction <= (lam_min - sigma)/(lam_2 - sigma)); Rayleigh quotient + residual reported.
 3. Mode profile of the eigenvector by Legendre degree d = 2i + parity.
usage: h68b_bisect_vec.py dump.json parity lo_exp hi_exp out.json   (lo_exp: Cholesky succeeds at -10^lo_exp; hi_exp fails)"""
import sys, json, time, math
from flint import arb, arb_mat, ctx
ctx.prec = 640; t0 = time.time(); bs = 64
fn, par, lo, hi, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
S = json.load(open(fn)); K = len(S)
Amid = [[arb(s).mid() for s in row] for row in S]
def chol_ok(l0):
    M = [[(Amid[i][j] - (l0 if i == j else 0)).mid() for j in range(K)] for i in range(K)]
    for j0 in range(0, K, bs):
        b = min(bs, K - j0); Rb = [[arb(0)] * b for _ in range(b)]
        for k in range(b):
            s = M[j0 + k][j0 + k]
            for t in range(k): s = (s - Rb[t][k] * Rb[t][k]).mid()
            if not (s > 0): return False, k + j0
            Rb[k][k] = s.sqrt().mid()
            for i in range(k + 1, b):
                s = M[j0 + k][j0 + i]
                for t in range(k): s = (s - Rb[t][k] * Rb[t][i]).mid()
                Rb[k][i] = (s / Rb[k][k]).mid()
        n2 = K - j0 - b
        if n2 == 0: break
        A12 = arb_mat([[M[j0 + k][j0 + b + c] for c in range(n2)] for k in range(b)])
        R11T = arb_mat([[Rb[c][r] for c in range(b)] for r in range(b)])
        R12 = R11T.solve(A12, algorithm='approx').mid(); S22 = (R12.transpose() * R12).mid()
        for r in range(n2):
            for c in range(n2): M[j0 + b + r][j0 + b + c] = (M[j0 + b + r][j0 + b + c] - S22[r, c]).mid()
    return True, None
log = []
while not chol_ok(-arb(10) ** arb(lo))[0]:
    log.append((lo, False, None)); hi, lo = lo, lo + 10; print('widen lo ->', lo, flush=True)
for _ in range(40):   # lo: success exponent (more negative l0), hi: failure exponent
    if lo - hi < 4e-4: break
    mid = (lo + hi) / 2; ok, piv = chol_ok(-arb(10) ** arb(mid))
    log.append((mid, ok, piv)); print('l0=-1e%.4f ok=%s piv=%s [%.0fs]' % (mid, ok, piv, time.time() - t0), flush=True)
    if ok: lo = mid
    else: hi = mid
sig = -arb(10) ** arb(lo)
I = arb_mat([[arb(1) if i == j else arb(0) for j in range(K)] for i in range(K)])
T = arb_mat([[Amid[i][j] for j in range(K)] for i in range(K)])
Mi = (T - I * sig).solve(I, algorithm='approx').mid()
v = arb_mat([[arb(1)] for i in range(K)])
for _ in range(30):
    u = (Mi * v).mid(); nu = sum((u[i, 0] ** 2 for i in range(K)), arb(0)).sqrt(); v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(K)])
Tv = (T * v).mid(); lam = sum((v[i, 0] * Tv[i, 0] for i in range(K)), arb(0))
res = sum(((Tv[i, 0] - lam * v[i, 0]) ** 2 for i in range(K)), arb(0)).sqrt()
w = [float((v[i, 0] ** 2).mid()) for i in range(K)]; deg = [2 * i + par for i in range(K)]
top = sorted(range(K), key=lambda i: -w[i])[:15]
d = dict(dump=fn, parity=par, K=K, bracket=[-10 ** lo, -10 ** hi], bisect=log, lam_min=float(lam.mid()), resid=float(res.mid()),
         weight_first={k: sum(w[:k]) for k in [5, 10, 20, 40, 80, 160, 320]}, weight_last={k: sum(w[-k:]) for k in [20, 50, 100, 200, 316]},
         top_modes=[(deg[i], w[i]) for i in top], centroid_deg=sum(deg[i] * w[i] for i in range(K)),
         degree_at_mass={q: next(deg[i] for i in range(K) if sum(w[:i + 1]) >= q) for q in [0.5, 0.9, 0.99]}, seconds=time.time() - t0)
print(json.dumps(d, indent=1)); json.dump(d, open(out, 'w'), indent=1)
