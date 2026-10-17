"""Grok S62 H62a (simplicity): second even min-max value of the zeta Weil form on L2[-27/20,27/20], lower bound.
Input: S37 even certificate dump Atil (780 orthonormal-Legendre even modes, reduced one-stroke form R_T, T=700),
  sha256 61bb468c...; R <= Q as forms on all of even L2 (Zhu Thm 1.1 + invisible lift), so lambda_2(Q) >= lambda_2(R).
Courant-Fischer: for ANY vector w in the leading block, lambda_2(R) >= min_{x perp (w,0)} R(x)/|x|^2
  >= min( lambda_min(Atil restricted to w-perp) - rE - delta_str, lamD ) - eps_B   (same two-block bound as the cert).
Restriction parametrised EXACTLY: w rational, j0=argmax|w|, columns p_i=e_i - v_i e_j0 (v=w/w_j0), i!=j0;
  P^T Atil P - s2 P^T P >= 0  <=>  lambda_min(M') >= s2 with M' = P^T X P - s2 v v^T (P^T P = I + v v^T).
This script: parses X (dump strings) in Arb, finds w and lambda_1, lambda_2 estimates (block inverse iteration +
Rayleigh-Ritz, non-rigorous; any w is valid), builds M' in Arb, dumps mid(M') strings + its max radius.
The PD check of M' - s2 I is then done by the two independent tools s36_cholres.py and s37_ldlt_exact.py.
usage: lam2_dump.py dump.json prec out.json"""
import sys, json, time, hashlib, random
from flint import arb, arb_mat, fmpq, ctx
import mpmath as mp
fn, P, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
ctx.prec = P; mp.mp.dps = int(P * 0.30); t0 = time.time()
raw = open(fn, 'rb').read(); sha = hashlib.sha256(raw).hexdigest(); S = json.loads(raw); del raw
K = len(S); X = arb_mat(K, K)
for i in range(K):
    row = S[i]
    for j in range(K): X[i, j] = arb(row[j])
    S[i] = None
del S
print(f'K={K} sha={sha[:12]} parsed [{time.time()-t0:.0f}s]', flush=True)
prad = 0.0
for i in range(K):            # symmetrise in place; record the parse radius (added to the final error budget)
    for j in range(i, K):
        prad = max(prad, float(X[i, j].rad()), float(X[j, i].rad()))
        m = ((X[i, j] + X[j, i]) / 2).mid(); X[i, j] = m; X[j, i] = m
random.seed(62); nb = 4
B = arb_mat([[random.uniform(-1, 1) for _ in range(nb)] for _ in range(K)])
Y = B
for it in range(3):
    Y = X.solve(Y, algorithm='approx').mid()
    mx = max(abs(float(Y[i, j])) for i in range(K) for j in range(nb)); Y = (Y * arb(1 / mx)).mid()
    print(f' inverse iteration {it} [{time.time()-t0:.0f}s]', flush=True)
G = (Y.transpose() * Y).mid(); H = (Y.transpose() * X * Y).mid()
tm = lambda M: mp.matrix([[mp.mpf(M[i, j].mid().str(mp.mp.dps + 10, radius=False)) for j in range(M.ncols())] for i in range(M.nrows())])
Gm, Hm = tm(G), tm(H); Lc = mp.cholesky(Gm); Li = mp.inverse(Lc); C = Li * Hm * Li.T; C = (C + C.T) / 2
ev, U = mp.eigsy(C); o = sorted(range(nb), key=lambda t: ev[t]); c1 = Li.T * U[:, o[0]]
lam1, lam2 = ev[o[0]], ev[o[1]]
print(f' Ritz estimates lam1={mp.nstr(lam1,8)} lam2={mp.nstr(lam2,8)} lam3={mp.nstr(ev[o[2]],8)}', flush=True)
Ym = tm(Y); w = Ym * c1; sc = max(abs(t) for t in w)
wq = [fmpq(int(mp.nint(t / sc * mp.mpf(10) ** 80)), 10 ** 80) for t in w]
j0 = max(range(K), key=lambda i: abs(wq[i])); v = [arb(q.p) / q.q / (arb(wq[j0].p) / wq[j0].q) for q in wq]
s2 = float(mp.nstr(lam2 * mp.mpf('0.98'), 4))
idx = [i for i in range(K) if i != j0]; n = len(idx); s2a = arb(repr(s2))
a = [X[i, j0] for i in range(K)]; xjj = X[j0, j0]
Mp = arb_mat(n, n); mr = 0.0
for r, i in enumerate(idx):
    for c, k in enumerate(idx):
        val = X[i, k] - v[k] * a[i] - v[i] * a[k] + (xjj - s2a) * v[i] * v[k]
        Mp[r, c] = val; mr = max(mr, float(val.rad()))
nd = int(P * 0.30) - 5
rows = [[Mp[r, c].mid().str(nd, radius=False) for c in range(n)] for r in range(n)]
txt = json.dumps(rows).encode(); open(out, 'wb').write(txt); osha = hashlib.sha256(txt).hexdigest()
meta = {'dump': fn, 'dump_sha256': sha, 'K': K, 'prec': P, 'lam1_est': mp.nstr(lam1, 12), 'lam2_est': mp.nstr(lam2, 12),
        'lam3_est': mp.nstr(ev[o[2]], 12), 'j0': j0, 's2': s2, 'Mprime': out, 'Mprime_sha256': osha, 'Mprime_n': n,
        'Mprime_maxrad': mr, 'parse_rad_max': prad, 'parse_rad_spectral_bound': K * prad * (1 + 2 * 3) , 'Mprime_rad_spectral_bound': n * mr, 'w': [str(q) for q in wq], 'sec': round(time.time() - t0, 1)}
json.dump(meta, open(out.replace('.json', '_meta.json'), 'w'), indent=1); print(json.dumps({k: meta[k] for k in meta if k != 'w'}), flush=True)
