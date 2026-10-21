"""Grok S76 H75b' (refined two-level): rigorous PD certificate for T - cI with a SMALL exact low block.
 Step 1 (exact, s71/s37_ldlt_exact.py on sub-dump A = T[:n,:n], sha-pinned): lambda_min(A) >= a_cert > c.
 Step 2 (Arb balls, prec P, rigorous arb_mat.solve; dump entries widened by 10^-157 as in s37):
        S = C - cI - B^T (A - cI)^-1 B as a ball matrix (m x m, m = K - n).
 Step 3 (S is well conditioned, NUMERICAL lambda_min(S) ~ 0.5): Sf = float64(mid S); float Cholesky Sf - tI = L L^T;
        exact-dyadic residual R = Sf - tI - L L^T in Arb (balls, prec 400); then rigorously
        lambda_min(S) >= t - ||R||_F - ||S - Sf||_F  (Weyl; ||S - Sf||_F <= sqrt(sum (|mid - Sf| + rad)^2)).
 If > 0: S PD and A - cI PD => T - cI PD.  CERTIFIED (computer-assisted, not Lean). RH OPEN.
usage: h75b_cert2.py dump.json sha256 n c a_cert t [prec]"""
import sys, json, hashlib, time
from fractions import Fraction
import numpy as np
from flint import arb, arb_mat, ctx
fn, sha_exp, n, c, acert, t = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], float(sys.argv[6])
ctx.prec = int(sys.argv[7]) if len(sys.argv) > 7 else 1024
t0 = time.time(); raw = open(fn, 'rb').read(); assert hashlib.sha256(raw).hexdigest() == sha_exp
rows = json.loads(raw); del raw; K = len(rows); m = K - n
RAD = arb(0, arb('1e-157')); cc = arb(c)
def e(i, j): return arb(rows[i][j]) + RAD - (cc if i == j else 0)
A = arb_mat([[e(i, j) for j in range(n)] for i in range(n)])
B = arb_mat([[e(i, j) for j in range(n, K)] for i in range(n)])
C = arb_mat([[e(i, j) for j in range(n, K)] for i in range(n, K)])
del rows
X = A.solve(B); S = C - B.transpose() * X; del X, A, B, C
t1 = time.time()
Sf = np.array([[float(S[i, j].mid()) for j in range(m)] for i in range(m)]); Sf = (Sf + Sf.T) / 2
dev2 = arb(0)
for i in range(m):
    for j in range(m):
        d = abs(S[i, j].mid() - arb(Sf[i, j])) + S[i, j].rad(); dev2 += d * d
dev = dev2.sqrt()
lam_f = float(np.linalg.eigvalsh(Sf).min())
res = dict(dump=fn, n=n, m=m, c=c, a_cert=acert, t=t, prec=ctx.prec, S_dev_F=float(dev.upper()), lamS_float=lam_f)
try:
    L = np.linalg.cholesky(Sf - t * np.eye(m))
    ctx.prec = 400
    La = arb_mat([[arb(float(L[i, j])) for j in range(m)] for i in range(m)])
    M = arb_mat([[arb(float(Sf[i, j])) - (arb(t) if i == j else 0) for j in range(m)] for i in range(m)])
    R = M - La * La.transpose()
    RF = sum((R[i, j] ** 2 for i in range(m) for j in range(m)), arb(0)).sqrt()
    bound = arb(t) - RF - dev
    res.update(R_F=float(RF.upper()), lamS_lower=float(bound.lower()), S_PD=bool(bound.lower() > 0))
except np.linalg.LinAlgError:
    res.update(S_PD=False, note='float Cholesky of Sf - tI failed')
res['verdict'] = 'T - cI PD given step 1 (lambda_min(A) >= a_cert > c)' if res['S_PD'] else 'FAIL at step 3'
res['secs'] = round(time.time() - t0); print(json.dumps(res), flush=True)
