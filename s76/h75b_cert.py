"""Grok S76 H75b: rigorous two-level PSD certificate for T - c I (one parity of a Mori-route dump).
T = [[A,B],[B^T,C]], A = T[:n,:n].  T - cI is PD  <=>  A - cI PD  and  S_c = C - cI - B^T (A - cI)^-1 B PD.
 Step 1 (exact, s71/s37_ldlt_exact.py on the sub-dump A, sha-pinned): lambda_min(A) >= a_cert; need a_cert > c.
 Step 2 (Arb balls, prec P, rigorous arb_mat.solve, every dump entry widened by 10^-157 as in s37):
        S_c as a ball matrix; Gershgorin lower bound g = min_i (lo(S_ii) - sum_{j!=i} up|S_ij|) > 0  =>  S_c PD.
 Then T - cI PD, so lambda_min(T) > c.  CERTIFIED (computer-assisted, not Lean).  RH OPEN.
usage: h75b_cert.py dump.json sha256 n c a_cert [prec]"""
import sys, json, hashlib, time
from fractions import Fraction
from flint import arb, arb_mat, ctx
fn, sha_exp, n, c, acert = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
ctx.prec = int(sys.argv[6]) if len(sys.argv) > 6 else 1024
assert Fraction(acert) > Fraction(c), 'need a_cert > c'
t0 = time.time(); raw = open(fn, 'rb').read(); assert hashlib.sha256(raw).hexdigest() == sha_exp
rows = json.loads(raw); del raw; K = len(rows); m = K - n
RAD = arb(0, arb('1e-157')); cc = arb(c)
def e(i, j): return arb(rows[i][j]) + RAD - (cc if i == j else 0)
A = arb_mat([[e(i, j) for j in range(n)] for i in range(n)])
B = arb_mat([[e(i, j) for j in range(n, K)] for i in range(n)])
C = arb_mat([[e(i, j) for j in range(n, K)] for i in range(n, K)])
del rows
X = A.solve(B)  # rigorous ball solve (raises if A - cI cannot be certified invertible at this precision)
S = C - B.transpose() * X
maxrad = max(float(S[i, j].rad()) for i in range(m) for j in range(m))
g = min(float((S[i, i] - sum((abs(S[i, j]) for j in range(m) if j != i), arb(0))).lower()) for i in range(m))
ok = g > 0
print(json.dumps(dict(dump=fn, n=n, m=m, c=c, a_cert=acert, prec=ctx.prec, S_maxrad=maxrad, gersh_S_lower=g,
                      S_PD=ok, verdict=('T - cI PD (given step 1)' if ok else 'FAIL: Gershgorin on S_c not positive'),
                      secs=round(time.time() - t0))), flush=True)
