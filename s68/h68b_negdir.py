"""Grok S68 H68b: locate the negative direction of the L=71/50 Mori-route T (S67 dumps, m=1/2, Omega=536, s=763/200, N=1272).
NUMERICAL (640-bit midpoint arithmetic on the dumped T; the sign itself is CERTIFIED in S67 by exact Rayleigh witnesses).
 1. Ritz: eigenvalues of leading k x k blocks (acb_mat.eig approx) -> negative count and min vs k (low-mode vs truncation edge).
 2. Full T: shifted inverse iteration at sigma = 1.5 * (most negative Ritz value) -> eigenvector of lambda_min(T);
    residual ||Tv - lam v|| reported; mode-weight profile |v_i|^2 by Legendre degree d = 2i + parity.
 3. Unshifted inverse iteration (S66 lam_near0) converges to the eigenvalue of smallest |lambda|: shown explicitly.
usage: h68b_negdir.py dump.json parity out.json"""
import sys, json, time
from flint import arb, arb_mat, acb_mat, ctx
ctx.prec = 640; t0 = time.time()
fn, par, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
S = json.load(open(fn)); n = len(S)
T = arb_mat([[arb(x).mid() for x in r] for r in S])
def block(k): return arb_mat([[T[i, j] for j in range(k)] for i in range(k)])
res = dict(dump=fn, parity=par, n=n, ritz={})
for k in [40, 60, 70, 80, 90, 100, 120, 160, 240, 320]:
    E = acb_mat(block(k)).eig(algorithm='approx')
    ev = sorted(float(e.real.mid()) for e in E)
    neg = [e for e in ev if e < 0]
    res['ritz'][k] = dict(n_neg=len(neg), min=ev[0], second=ev[1], neg_list=neg[:6])
    print(k, len(neg), '%.4e' % ev[0], '%.4e' % ev[1], '[%.0fs]' % (time.time() - t0), flush=True)
def invit(sig, its=40):
    I = arb_mat([[arb(1) if i == j else arb(0) for j in range(n)] for i in range(n)])
    M = T - I * sig
    Mi = M.solve(I, algorithm='approx').mid()
    v = arb_mat([[arb(1) / (i + 1)] for i in range(n)])
    for _ in range(its):
        u = (Mi * v).mid(); nu = sum((u[i, 0] ** 2 for i in range(n)), arb(0)).sqrt()
        v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(n)])
    Tv = (T * v).mid(); lam = sum((v[i, 0] * Tv[i, 0] for i in range(n)), arb(0))
    r = sum(((Tv[i, 0] - lam * v[i, 0]) ** 2 for i in range(n)), arb(0)).sqrt()
    return lam, r, v
sig = arb(1.5 * min(res['ritz'][k]['min'] for k in res['ritz']))
lam, r, v = invit(sig)
w = [float((v[i, 0] ** 2).mid()) for i in range(n)]
deg = [2 * i + par for i in range(n)]
cum = lambda K: sum(w[:K])
top = sorted(range(n), key=lambda i: -w[i])[:12]
res['full'] = dict(sigma=float(sig.mid()), lam_min=float(lam.mid()), resid=float(r.mid()),
                   weight_first=dict((K, cum(K)) for K in [10, 20, 40, 60, 80, 100, 120, 200]),
                   weight_last100=sum(w[-100:]), top_modes=[(deg[i], w[i]) for i in top],
                   centroid_deg=sum(deg[i] * w[i] for i in range(n)))
print(json.dumps(res['full']), flush=True)
lam0, r0, _ = invit(arb(0), 60)
res['nearest_zero'] = dict(lam=float(lam0.mid()), resid=float(r0.mid()))
print('nearest-zero (S66 lam_near0):', res['nearest_zero'], '[%.0fs]' % (time.time() - t0), flush=True)
json.dump(res, open(out, 'w'), indent=1)
