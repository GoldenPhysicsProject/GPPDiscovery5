"""Grok S66: recompute the Mori finite condition T from a dump (A, B, J) with a chosen X.
X = 'opt': X = (B - mA)^{-1}(A - mI), the Loewner minimiser of W(X) = m^{-1}I + X*CX - DX - X*D
(C = m^{-1}(B - mA), D = m^{-1}A - I; completing the square: W(X) - W(X_opt) = (X - X_opt)*C(X - X_opt) >= 0).
X = 'ls': X = B^{-1}A.  NUMERICAL (midpoint arithmetic)."""
import sys, pickle, json, argparse
from flint import arb, arb_mat, ctx
sys.path.insert(0, '.')
from mori_route import lam_near0
from mori_consts import consts
ap = argparse.ArgumentParser(); ap.add_argument('dump'); ap.add_argument('--m', default='1/2'); ap.add_argument('--X', default='opt')
ap.add_argument('--prec', type=int, default=448); ap.add_argument('--alpha', default='0'); ap.add_argument('--beta', default='0')
ap.add_argument('--dumpT', default='')
a = ap.parse_args(); ctx.prec = a.prec
D = pickle.load(open(a.dump, 'rb')); A = arb_mat([[arb(x) for x in r] for r in D['A']]); B = arb_mat([[arb(x) for x in r] for r in D['B']])
J = arb_mat([[arb(x) for x in r] for r in D['J']]); n = A.nrows(); m = arb(a.m)
I = arb_mat([[arb(1) if i == j else arb(0) for j in range(n)] for i in range(n)])
X = ((B - m * A).solve(A - m * I) if a.X == 'opt' else B.solve(A)).mid(); XT = X.transpose()
W = X + XT - XT * A * X + (I - A * X - XT * A + XT * B * X) / m; W = (W + W.transpose()) / 2
T = W.inv() + J - arb(a.alpha) * I - arb(a.beta) * (B - A * A); T = ((T + T.transpose()) / 2).mid()
lam, v = lam_near0(T.tolist())
print(json.dumps(dict(dump=a.dump, X=a.X, lam_T=lam.str(12))))
if a.dumpT:
    pickle.dump([[T[i, j].str(n=int(a.prec * 0.302) + 5, radius=False) for j in range(n)] for i in range(n)], open(a.dumpT, 'wb'))
