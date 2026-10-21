"""Grok S76 H75b: two-level (low-block LDL^T + Schur/Gershgorin high block) PSD test on a dumped Mori-route
matrix T (one parity, K x K).  Split T = [[A,B],[B^T,C]], A = T[:n,:n].
Sufficient test (H75b):  lambda_min(T) >= 0 if  a := lambda_min(A) > 0  and  g_C - ||B||_F^2 / a > 0,
  g_C = min_i (C_ii - sum_{j!=i} |C_ij|)  (Gershgorin lower bound on lambda_min(C)),
  since S = C - B^T A^-1 B >= lambda_min(C) - ||B||_2^2/a  and  T PSD <=> A PD and S PSD.
CERTIFIED parts: g_C (upper/lower ends) and ||B||_F (lower/upper ends) in Arb balls; each entry gets an extra radius
 10^-157 (dump strings are printed to e_last <= -158, see s37).  a: certified LOWER bound via Cauchy interlacing,
 a >= lambda_min(T) >= l_cert (s37 exact LDL^T); NUMERICAL estimate of lambda_min(A) by inverse iteration.
 Since a <= lambda_min(A) the test with the TRUE lambda_min(A) is the most favourable case; we report both.
NUMERICAL diagnostics: Schur complement S (approx solve, prec 640): min diag, Gershgorin, lambda_min(S);
 scaled variant D T D, D = diag(T)^-1/2: Gershgorin of scaled C, and rho(|offdiag scaled C|) (generalized diagonal
 dominance with the best positive diagonal scaling holds iff this Perron root < 1 -- float, NUMERICAL).
usage: h75b_two_level.py dump.json sha256 l_cert n1,n2,... [shift]   (shift adds shift*I: sanity positive control)
RH OPEN."""
import sys, json, hashlib, time
import numpy as np
from flint import arb, arb_mat, ctx
ctx.prec = 640
fn, sha_exp, lcert, ns = sys.argv[1], sys.argv[2], float(sys.argv[3]), [int(x) for x in sys.argv[4].split(',')]
shift = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
t0 = time.time(); raw = open(fn, 'rb').read(); sha = hashlib.sha256(raw).hexdigest(); assert sha == sha_exp, sha
rows = json.loads(raw); del raw; K = len(rows)
RAD = arb(0, arb('1e-157'))
T = arb_mat([[arb(x) + RAD + (arb(shift) if i == j else 0) for j, x in enumerate(r)] for i, r in enumerate(rows)])
Tf = np.array([[float(x) for x in r] for r in rows]) + shift * np.eye(K); del rows
print(f'K={K} sha OK shift={shift} [{time.time()-t0:.0f}s]', flush=True)
print(f'NUMERICAL float: diag range [{Tf.diagonal().min():.3e},{Tf.diagonal().max():.3e}], max|T_ij| {abs(Tf).max():.3e}, ||T||_F {np.linalg.norm(Tf):.3e}', flush=True)

def sub(M, r0, r1, c0, c1):
    return arb_mat([[M[i, j] for j in range(c0, c1)] for i in range(r0, r1)])
def lam_inv_iter(Mi, m, its=40):  # Mi approx inverse (mid), NUMERICAL lambda_min
    v = arb_mat([[arb(1) / (i + 1)] for i in range(m)]); lam = None
    for _ in range(its):
        u = (Mi * v).mid(); nu = sum((u[i, 0] ** 2 for i in range(m)), arb(0)).sqrt()
        lam = sum((v[i, 0] ** 2 for i in range(m)), arb(0)) / sum((u[i, 0] * v[i, 0] for i in range(m)), arb(0))
        v = arb_mat([[(u[i, 0] / nu).mid()] for i in range(m)])
    return lam.mid(), v
out = []
for n in ns:
    t1 = time.time(); m = K - n
    # ---- certified pieces
    gC = None
    for i in range(n, K):
        s = T[i, i] - sum((abs(T[i, j]) for j in range(n, K) if j != i), arb(0))
        gC = s if gC is None else (s if s.upper() < gC.upper() else gC)  # min (by upper end; we report the ball)
    BF2 = sum((T[i, j] ** 2 for i in range(n) for j in range(n, K)), arb(0))
    # scaled Gershgorin of C (certified ball)
    gCs = None
    for i in range(n, K):
        s = 1 - sum((abs(T[i, j]) / (T[i, i] * T[j, j]).sqrt() for j in range(n, K) if j != i), arb(0))
        gCs = s if gCs is None else (s if s.upper() < gCs.upper() else gCs)
    # ---- numerical pieces
    A = sub(T, 0, n, 0, n); B = sub(T, 0, n, n, K); C = sub(T, n, K, n, K)
    I_n = arb_mat([[arb(1) if i == j else arb(0) for j in range(n)] for i in range(n)])
    X = A.solve(arb_mat([[I_n[i, j] for j in range(n)] + [B[i, j] for j in range(m)] for i in range(n)]), algorithm='approx').mid()
    Ai = arb_mat([[X[i, j] for j in range(n)] for i in range(n)])
    AiB = arb_mat([[X[i, n + j] for j in range(m)] for i in range(n)])
    lamA, _ = lam_inv_iter(Ai, n)
    BAB = (B.transpose() * AiB).mid(); trBAB = float(sum((BAB[i, i] for i in range(m)), arb(0)))
    S = (C - BAB).mid()
    Sf = np.array([[float(S[i, j]) for j in range(m)] for i in range(m)])
    gS = min(Sf[i, i] - (np.abs(Sf[i]).sum() - abs(Sf[i, i])) for i in range(m))
    if m > 1:
        Si = S.solve(arb_mat([[arb(1) if i == j else arb(0) for j in range(m)] for i in range(m)]), algorithm='approx').mid()
        lamS, _ = lam_inv_iter(Si, m)
    else:
        lamS = S[0, 0]
    Cf = Tf[n:, n:]; d = 1 / np.sqrt(Cf.diagonal()); Cs = d[:, None] * Cf * d[None, :]
    O = np.abs(Cs - np.diag(Cs.diagonal())); rho = max(abs(np.linalg.eigvals(O))) if m > 1 else 0.0
    lamC = np.linalg.eigvalsh(Cf).min()
    BF2u = float(BF2.upper()); aU = float(lamA)
    rec = dict(n=n, m=m,
               gC=gC.str(6), gC_upper=float(gC.upper()), BF=float(BF2.upper()) ** 0.5,
               a_cert=lcert, a_num=float(lamA),
               test_cert=float(gC.upper()) - float(BF2.lower()) / lcert,
               test_best=float(gC.upper()) - float(BF2.lower()) / aU,
               BF2_over_a=float(BF2.lower()) / aU,
               scaled_gC=gCs.str(6), scaled_gC_upper=float(gCs.upper()), rho_scaled_offdiag=float(rho),
               lamC_float=float(lamC), S_mindiag=float(Sf.diagonal().min()), S_gersh=float(gS), lamS_num=float(lamS),
               tr_BtAinvB=trBAB, test_refined_num=float(gC.upper()) - trBAB,
               secs=round(time.time() - t1))
    out.append(rec); print(json.dumps(rec), flush=True)
json.dump(dict(dump=fn, sha256=sha, K=K, shift=shift, l_cert=lcert, results=out), open(sys.argv[6] if len(sys.argv) > 6 else '/dev/null', 'w'), indent=1)
print('DONE', round(time.time() - t0))
