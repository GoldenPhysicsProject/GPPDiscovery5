"""Grok S67-resumed: negative-eigenvalue witness for a dumped T (proves T NOT PSD, so the Mori finite condition fails
at those parameters).  v = bottom eigenvector of the leading k x k block (mpmath, any v works), rounded to exact
rationals; Rayleigh q = v^T T_dump v / v^T v computed EXACTLY (fractions) from the decimal strings.
Then lambda_min(T_true) <= q + err, err = ||T_true - T_dump||_2 (radT row sum + inverse error + quad bound + string
error, passed in).  usage: neg_witness.py dump.json k err_upper"""
import sys, json, mpmath as mp
from fractions import Fraction as Fr
fn, k, err = sys.argv[1], int(sys.argv[2]), Fr(sys.argv[3])
mp.mp.dps = 200
S = json.load(open(fn))
A = mp.matrix([[mp.mpf(x) for x in r[:k]] for r in S[:k]])
E, Q = mp.eigsy(A); i0 = min(range(k), key=lambda i: E[i])
v = [Fr(int(mp.nint(Q[i, i0] * mp.mpf(2) ** 600)), 2 ** 600) for i in range(k)]
T = [[Fr(x) for x in r[:k]] for r in S[:k]]
num = sum(v[i] * sum(T[i][j] * v[j] for j in range(k)) for i in range(k)); den = sum(x * x for x in v)
q = num / den; up = q + err
print(json.dumps(dict(dump=fn, k=k, eig_block_min=mp.nstr(E[i0], 8), rayleigh_exact=float(q), err=float(err),
                      lam_min_upper=float(up), NOT_PSD=bool(up < 0))))
