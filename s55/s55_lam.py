"""Grok S55: FIX for the S36-S54 `lam_estimate` bug (inverse iteration returns the eigenvalue of smallest MODULUS,
so it printed +1.26e-70 for the indefinite S54 L=1.354 even matrix).
lam_estimate_fixed(Atil) = (mu, n_neg):
  n_neg = number of negative pivots of a blocked LDL^T of the midpoint matrix (Sylvester inertia; non-rigorous
          float-ball midpoints at the working precision, used for REPORTING only),
  mu    = smallest-modulus eigenvalue (inverse iteration), which equals the smallest ALGEBRAIC eigenvalue iff n_neg == 0.
  If n_neg > 0 the function returns lam = -inf-flag (None) and the caller must not read mu as lam_min.
PD certificates are unaffected: they were always decided by s36_cholres + s37_ldlt_exact with a failing negative control.
usage (on a dump): s55_lam.py atil.json [prec] [bs]"""
import sys, json, time
from flint import arb, arb_mat, ctx

def inertia_neg(rows, bs=48):
    M = arb_mat(rows).mid(); n = M.nrows(); neg = 0; zero = 0; off = 0
    while n > 0:
        b = min(bs, n)
        A11 = [[M[i, j] for j in range(b)] for i in range(b)]
        # scalar LDL^T of the small block (no pivoting; midpoints)
        Lm = [[arb(0)] * b for _ in range(b)]; d = [arb(0)] * b
        for k in range(b):
            s = A11[k][k]
            for j in range(k): s -= Lm[k][j] * Lm[k][j] * d[j]
            d[k] = s.mid()
            if d[k] == 0: zero += 1; d[k] = arb(10) ** (-ctx.prec // 4)
            elif d[k] < 0: neg += 1
            Lm[k][k] = arb(1)
            for i in range(k + 1, b):
                s = A11[i][k]
                for j in range(k): s -= Lm[i][j] * Lm[k][j] * d[j]
                Lm[i][k] = (s / d[k]).mid()
        if n == b: break
        A11m = arb_mat([[M[i, j] for j in range(b)] for i in range(b)])
        A12 = arb_mat([[M[i, j] for j in range(b, n)] for i in range(b)])
        A21 = A12.transpose()
        A22 = arb_mat([[M[i, j] for j in range(b, n)] for i in range(b, n)])
        M = (A22 - A21 * A11m.solve(A12).mid()).mid(); n -= b
    return neg, zero

def smallest_modulus(rows, iters=60):
    K = len(rows); Ai = arb_mat(rows).inv().mid()
    v = arb_mat([[arb(1) / (i + 1)] for i in range(K)]); lam = None
    for _ in range(iters):
        u = Ai * v
        nu = sum((u[i, 0] ** 2 for i in range(K)), arb(0)).sqrt()
        rq = sum((u[i, 0] * v[i, 0] for i in range(K)), arb(0)) / sum((v[i, 0] ** 2 for i in range(K)), arb(0))
        lam = 1 / rq; v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(K)])
    return float(lam.mid())

def lam_estimate_fixed(rows, bs=48):
    neg, zero = inertia_neg(rows, bs)
    mu = smallest_modulus(rows)
    return dict(mu_smallest_modulus=mu, n_neg_pivots=neg, n_zero_pivots=zero,
                lam_min_estimate=(mu if neg == 0 and zero == 0 else None),
                verdict=('PD-looking (mu = lam_min)' if neg == 0 and zero == 0 else f'INDEFINITE: {neg} negative pivots; mu is NOT lam_min'))

if __name__ == '__main__':
    f = sys.argv[1]; ctx.prec = int(sys.argv[2]) if len(sys.argv) > 2 else 640; bs = int(sys.argv[3]) if len(sys.argv) > 3 else 48
    t0 = time.time(); rows = [[arb(x) for x in r] for r in json.load(open(f))]
    out = lam_estimate_fixed(rows, bs); out.update(file=f, K=len(rows), prec=ctx.prec, secs=round(time.time() - t0))
    print(json.dumps(out), flush=True)
