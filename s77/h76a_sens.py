"""Grok S77 H76a-0: precision budget of the two-level (H75b') certificate, using Astra's weighted Schur bound
(PR #11 comment 6328829674, bridge 6fdaeb8).  Simulates low-precision assembly on an existing sha-pinned dump:
  B_hat = B rounded to the grid eps_B (so |B - B_hat| <= eps_B + 1e-157 entrywise, exactly known),
  C_hat = C rounded to the grid eps_C (|dC| <= eps_C + 1e-157 entrywise),  A exact (dump balls).
Rigorous chain (A exact, lambda_min(A) >= a_cert from the exact LDL^T step 1):
  S_hat = C_hat - B_hat^T A^-1 B_hat  (Arb balls)  ->  mu = certified lower bound of lambda_min(S_hat) (float Cholesky + Arb residual + Weyl)
  b^2 >= trace(B_hat^T A^-1 B_hat) (Arb upper),  e = sqrt(n m) eps_B' / sqrt(a_cert),  c = m eps_C'
  lambda_min(S_true) >= mu - c - (2 b e + e^2)            [Astra: ||B^T A^-1 B - Bh^T A^-1 Bh|| <= 2be + e^2]
If > 0: S_true PD, and with A PD, T PD.  CERTIFIED (computer-assisted). Simulation of the precision budget only:
the dump itself was assembled at 448 bits.  RH OPEN.
usage: h76a_sens.py dump sha n a_cert t epsB_list epsC_list"""
import sys, json, hashlib, time
import numpy as np
from flint import arb, arb_mat, ctx
fn, sha_exp, n, acert, t = sys.argv[1], sys.argv[2], int(sys.argv[3]), arb(sys.argv[4]), float(sys.argv[5])
epsBs = sys.argv[6].split(','); epsCs = sys.argv[7].split(',')
ctx.prec = 1024
raw = open(fn, 'rb').read(); assert hashlib.sha256(raw).hexdigest() == sha_exp
rows = json.loads(raw); del raw; K = len(rows); m = K - n
RAD = arb(0, arb('1e-157'))
A = arb_mat([[arb(rows[i][j]) + RAD for j in range(n)] for i in range(n)])
Bm = [[arb(rows[i][j]) for j in range(n, K)] for i in range(n)]
Cf = np.array([[float(arb(rows[i][j]).mid()) for j in range(n, K)] for i in range(n, K)]); assert np.abs(Cf).max() < 10, np.abs(Cf).max()
del rows
def rnd(x, eps):  # x rounded down to the grid eps (|x - rnd| < eps)
    qi = (x.mid() / eps).mid().floor().unique_fmpz()
    return arb(qi) * eps
def lower_lam(Sf_arb_mid, t):
    Sf = (Sf_arb_mid + Sf_arb_mid.T) / 2
    try: L = np.linalg.cholesky(Sf - t * np.eye(m))
    except np.linalg.LinAlgError: return None, float(np.linalg.eigvalsh(Sf).min()), Sf
    ctx.prec = 200
    La = arb_mat([[arb(float(L[i, j])) for j in range(m)] for i in range(m)])
    M = arb_mat([[arb(float(Sf[i, j])) - (arb(t) if i == j else 0) for j in range(m)] for i in range(m)])
    R = M - La * La.transpose(); RF = sum((R[i, j] ** 2 for i in range(m) for j in range(m)), arb(0)).sqrt()
    ctx.prec = 1024
    return arb(t) - RF, float(np.linalg.eigvalsh(Sf).min()), Sf
out = []
for eB in epsBs:
    t0 = time.time(); epsB = arb(eB)
    Bh = arb_mat([[rnd(x, epsB) for x in r] for r in Bm])
    X = A.solve(Bh); Q = Bh.transpose() * X; del X
    tr = sum((Q[i, i] for i in range(m)), arb(0)); b = arb(tr.upper()).sqrt()
    e = (arb(n * m).sqrt() * (epsB + arb('1e-157'))) / acert.sqrt()
    Qmid = np.array([[float(Q[i, j].mid()) for j in range(m)] for i in range(m)])
    # ||Q - Qmid||_F (balls + float rounding of midpoints)
    d2 = arb(0)
    for i in range(m):
        for j in range(m):
            d = abs(Q[i, j].mid() - arb(Qmid[i, j])) + Q[i, j].rad(); d2 += d * d
    dQ = d2.sqrt(); del Q
    for eC in epsCs:
        epsC = float(eC); Ch = np.round(Cf / epsC) * epsC if epsC > 0 else Cf
        # Ch - Cf deviation is the simulated assembly error; the certificate only uses Ch and the bound eps_C (+ float rounding of C, < 1e-15 rel)
        Sf = Ch - Qmid
        mu, lamf, Sfs = lower_lam(Sf, t)
        cC = arb(m) * (arb(eC) + arb('1e-14'))   # ||dC||_2 <= ||dC||_F <= m * max entry error (float-rounding of C absorbed by 1e-14, |C_ij| < 10)
        rec = dict(dump=fn.split('/')[-1], n=n, m=m, epsB=eB, epsC=eC, b=float(b.upper()), e=float(e.upper()), dQ_F=float(dQ.upper()),
                   lamS_float=lamf)
        if mu is None: rec.update(cert=False, note='float Cholesky of S_hat - tI failed')
        else:
            lb = mu - dQ - cC - (2 * b * e + e * e)
            rec.update(mu=float(mu.lower()), loss_B=float((2 * b * e + e * e).upper()), loss_C=float(cC.upper()), lamS_lower=float(lb.lower()), cert=bool(lb.lower() > 0))
        rec['secs'] = round(time.time() - t0); print(json.dumps(rec), flush=True); out.append(rec)
