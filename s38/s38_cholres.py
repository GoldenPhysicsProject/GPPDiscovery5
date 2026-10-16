"""Grok S36: rigorous PD certificate for a dumped midpoint Gram matrix WITHOUT ball-elimination radius growth.
Method (Rump-style 'approximate Cholesky + rigorous residual'):
  1. A_s = matrix parsed (Arb, prec P) from the dumped decimal strings (each a tiny ball); M = mid(A_s) - l0 I.
  2. Approximate (non-rigorous) blocked Cholesky M ~ R^T R in point arithmetic at prec P (any R works).
  3. RIGOROUS: E = (A_s - l0 I) - R^T R in ball arithmetic (R exact point matrix); ||E||_2 <= max_i sum_j |E_ij|.
     Since R^T R >= 0, every symmetric matrix in A_s satisfies lambda_min >= l0 - ||E||_2.
  4. String error: the dumped strings differ from the 384-bit midpoints Atil by <= 10^(e-n+1) per entry
     (n significant digits printed); delta_str = K * max that bound covers ||Atil - A_s||_2.
Output: certified lambda_min(Atil) >= l0 - ||E|| - delta_str (then subtract rE and eps_B as in the cert).
Status: computer-assisted (Arb), single implementation, NOT Lean.  RH OPEN.
usage: s36_cholres.py dump.json P fac1,fac2,... lam_est [bs]"""
import sys, json, time, math
from flint import arb, arb_mat, ctx
fn = sys.argv[1]; P = int(sys.argv[2]); facs = [float(x) for x in sys.argv[3].split(',')]; lam_est = float(sys.argv[4])
bs = int(sys.argv[5]) if len(sys.argv) > 5 else 64
ctx.prec = P; t0 = time.time()
S = json.load(open(fn)); K = len(S)
def strerr(s):
    m = s.lstrip('-'); mant, _, ex = m.partition('e'); ex = int(ex) if ex else 0
    digs = mant.replace('.', '').lstrip('0'); nd = len(digs)
    if nd == 0: return 0.0
    ip = mant.split('.')[0].lstrip('0'); 
    lead = (len(ip) - 1) if ip else -(len(mant.split('.')[1]) - len(mant.split('.')[1].lstrip('0')) + 1)
    return 10.0 ** (lead + ex - nd + 1)
dstr = K * max(strerr(s) for row in S for s in row)
A = [[arb(s) for s in row] for row in S]
Amid = [[x.mid() for x in row] for row in A]
print(f'K={K} prec={P} delta_str<={dstr:.3e} parse [{time.time()-t0:.0f}s]', flush=True)
def chol_approx(l0):
    M = [[Amid[i][j] - (l0 if i == j else 0) for j in range(K)] for i in range(K)]
    M = [[x.mid() for x in row] for row in M]
    R = [[arb(0)] * K for _ in range(K)]; minpiv = None
    for j0 in range(0, K, bs):
        b = min(bs, K - j0); idx = range(j0, j0 + b)
        # small dense Cholesky of the leading block (point arithmetic: mid after every op)
        Rb = [[arb(0)] * b for _ in range(b)]
        for k in range(b):
            s = M[j0 + k][j0 + k]
            for t in range(k): s = (s - Rb[t][k] * Rb[t][k]).mid()
            if not (s > 0): return None, k + j0
            minpiv = float(s) if minpiv is None or float(s) < minpiv else minpiv
            Rb[k][k] = s.sqrt().mid()
            for i in range(k + 1, b):
                s = M[j0 + k][j0 + i]
                for t in range(k): s = (s - Rb[t][k] * Rb[t][i]).mid()
                Rb[k][i] = (s / Rb[k][k]).mid()
        for k in range(b):
            for i in range(k, b): R[j0 + k][j0 + i] = Rb[k][i]
        n2 = K - j0 - b
        if n2 == 0: break
        A12 = arb_mat([[M[j0 + k][j0 + b + c] for c in range(n2)] for k in range(b)])
        R11T = arb_mat([[Rb[c][r] for c in range(b)] for r in range(b)])
        R12 = R11T.solve(A12, algorithm='approx').mid()
        for k in range(b):
            for c in range(n2): R[j0 + k][j0 + b + c] = R12[k, c]
        S22 = (R12.transpose() * R12).mid()
        for r in range(n2):
            Mr = M[j0 + b + r]
            for c in range(n2): Mr[j0 + b + c] = (Mr[j0 + b + c] - S22[r, c]).mid()
    return R, minpiv
res = {'dump': fn, 'K': K, 'prec': P, 'delta_str': dstr, 'lam_est': lam_est, 'tries': []}
for fac in facs:
    c = float('%.4g' % (lam_est * fac)); l0 = arb(repr(c)); t1 = time.time()
    R, info = chol_approx(l0)
    if R is None:
        print(f' fac={fac} l0={c:.4e}: approx Cholesky FAILED at pivot {info} (A - l0 I numerically indefinite) [{time.time()-t1:.0f}s]', flush=True)
        res['tries'].append({'fac': fac, 'l0': c, 'chol': 'fail', 'pivot': info}); continue
    Rm = arb_mat(R); Am = arb_mat(A)
    for i in range(K): Am[i, i] = Am[i, i] - l0
    E = Am - Rm.transpose() * Rm
    nE = max(sum((abs(E[i, j]).upper() for j in range(K)), arb(0)) for i in range(K))
    lam_c = l0 - arb(nE.upper()) - arb(dstr)
    ok = lam_c > 0
    print(f' fac={fac} l0={c:.4e}: minpiv={info:.3e} ||E||<={float(nE):.3e} => lambda_min(Atil) >= {float(lam_c.lower()):.6e} {"OK" if ok else "NO"} [{time.time()-t1:.0f}s]', flush=True)
    res['tries'].append({'fac': fac, 'l0': c, 'minpiv': info, 'normE': float(nE), 'lam_lower': float(lam_c.lower()), 'ok': bool(ok)})
    if ok: res['lam_lower'] = float(lam_c.lower()); res['l0'] = c; break
json.dump(res, open(fn.replace('.json', f'_cholres_p{P}.json'), 'w'), indent=1)
print('DONE', time.time() - t0, flush=True)
