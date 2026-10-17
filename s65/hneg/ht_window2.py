"""Grok S64 H64c stage 3b: same zero-side window test as ht_window.py but all tiers in Arb (prec 256), so the noise floor
drops from ~1e-21 to ~1e-60 (NUMERICAL: zero list + proxy caveat as before). usage: ht_window2.py t L K parity"""
import sys, json
from flint import arb, acb, arb_mat, acb_mat, ctx
import mpmath as mp
tt, Ls, K, par = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
import os
ctx.prec = int(os.environ.get('PREC','256')); mp.mp.dps = int(ctx.prec*0.29)
PSHIFT = os.environ.get('PSHIFT','0')
def parse(s):
    if '/' in s: p, q = s.split('/'); return arb(int(p))/int(q)
    return arb(s)
L = parse(Ls); PI = arb.pi()
zz = [l.split()[1] for l in open('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt')]
if tt == '0': real = [arb(g) for g in zz if float(g) <= 150]; cplx = []
else:
    real = [arb(repr(x))/2 for x in json.load(open(f'scan_t{tt}.json'))['real_zeros_x']]
    cplx = [acb(arb(repr(a)), arb(repr(b))) for a, b in json.load(open(f'cplx_t{tt}.json'))['gamma_form']]
prox = [arb(g)+arb(PSHIFT) for g in zz if float(g) > 150]
fr = [(2*m+1)*PI/(2*L) for m in range(K)] if par == 'even' else [m*PI/L for m in range(1, K+1)]
sgn = [(-1)**m for m in range(K)]
def vec_real(g):   # real zero: even -> real; odd -> i*real (|.|^2 same) -> use real part magnitude
    if par == 'even': c = (g*L).cos(); return [sgn[m]*c*2*a/(a*a-g*g) for m, a in enumerate(fr)]
    s = (g*L).sin(); return [-sgn[m]*s*2*b/(b*b-g*g) for m, b in enumerate(fr)]
def vec_c(z):
    if par == 'even': c = (z*L).cos(); return [sgn[m]*c*2*a/(a*a-z*z) for m, a in enumerate(fr)]
    s = (z*L).sin(); return [acb(0, 1)*(-sgn[m])*s*2*b/(b*b-z*z) for m, b in enumerate(fr)]
allz = real + prox; Gm = arb_mat(K, K)
for c0 in range(0, len(allz), 5000):           # chunked to bound memory
    V = arb_mat([vec_real(g) for g in allz[c0:c0+5000]]); Gm += V.transpose()*V*2   # +-gamma
    del V
G = [[Gm[i, j] for j in range(K)] for i in range(K)]
for z in cplx:
    for zc in (z, z.conjugate()):
        u = vec_c(zc); w = vec_c(zc.conjugate())
        for i in range(K):
            for j in range(K): G[i][j] += 2*(u[i]*w[j].conjugate()).real
M = mp.matrix(K, K)
for i in range(K):
    for j in range(K): M[i, j] = mp.mpf(((G[i][j]+G[j][i])/2/L).mid().str(mp.mp.dps-3, radius=False))
E, Q = mp.eigsy(M); lam = min(E); k = list(E).index(lam)
v = [arb(mp.nstr(Q[i, k], mp.mp.dps-3)) for i in range(K)]
def form(zs):
    s = arb(0)
    for z in zs:
        for zc in (z, z.conjugate()):
            u = vec_c(zc); w = vec_c(zc.conjugate())
            s += 2*(sum((u[i]*v[i] for i in range(K)), acb(0))*sum((w[i]*v[i] for i in range(K)), acb(0)).conjugate()).real
    return s/L
pc = form(cplx) if cplx else arb(0)
r = {"t": tt, "L": Ls, "K": K, "parity": par, "lam_min": mp.nstr(lam, 8), "pair_part_on_minvec": pc.mid().str(8, radius=False),
     "n_real": len(real), "n_cplx_pairs": len(cplx), "engine": "arb%d" % ctx.prec, "pshift": PSHIFT}
print(json.dumps(r)); open('window2.jsonl', 'a').write(json.dumps(r)+'\n')
