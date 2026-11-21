"""Round 34 controls: brute-force ratio-class layers (r34_pairs.c) for ZETA (cross-check of r34_layers.c), ZETA^2, F_0.1,
EPS, DH (each with its own Dirichlet-inverse mu; centering c = 1/b(1/2) from the R32 Arb constants) and Beurling
(R32 DMV-type primes, c = 1/zeta_B(1/2)). X = 3e4. Layer k = Omega(u)+Omega(v) of the reduced ratio. Property (a)
(reduced-ratio sign preservation) is tested directly on m <= 150th support element, n <= 5000. No zeros."""
import json, math, sys, subprocess, numpy as np
sys.path.insert(0, '../round32_mollified_floor'); from r22_models import KAP
X = 30000; L = math.log(X)
K = json.load(open('../round32_mollified_floor/r32_constants.json')); BZ = json.load(open('../round32_mollified_floor/r32_beur_zeta.json'))
spf = np.zeros(X + 1, int)
for p in range(2, X + 1):
    if spf[p] == 0: spf[p::p] = np.where(spf[p::p] == 0, p, spf[p::p])
def factor(n):
    out = []
    while n > 1: out.append(int(spf[n])); n //= spf[n]
    return out
facs = [None] + [factor(k) for k in range(1, X + 1)]
def inverse_fast(b):
    binv = np.zeros(X + 1); binv[1] = 1.0; acc = np.zeros(X + 1); nzb = np.nonzero(b[2:])[0] + 2
    for d in range(1, X + 1):
        if d > 1: binv[d] = -acc[d]
        if binv[d] != 0 and 2*d <= X:
            k = nzb[:np.searchsorted(nzb, X//d, side='right')]; acc[d*k] += binv[d]*b[k]
    return binv
def mult_mu(fp, sqfree=False):  # multiplicative: mu(p) = fp(p), mu(p^2) = 1 (cube-free) or 0 (squarefree)
    mu = np.ones(X + 1); mu[0] = 0
    for n in range(2, X + 1):
        fs = facs[n]; v = 1.0
        for p in set(fs):
            e = fs.count(p); v *= fp(p) if e == 1 else ((0.0 if sqfree else 1.0) if e == 2 else 0.0)
        mu[n] = v
    return mu
def records(vals, fl, logs):
    dt = np.dtype([('x', 'f8'), ('a', 'f8'), ('nf', 'i4'), ('p', 'i4', 24)], align=True); rec = np.zeros(len(vals), dt)
    rec['x'] = logs; rec['a'] = vals; rec['nf'] = [len(f) for f in fl]
    for i, fs in enumerate(fl): rec['p'][i, :len(fs)] = sorted(fs)
    o = np.argsort(rec['x'], kind='stable'); return rec[o]
def run(name, rec):
    rec.tofile(f'rec_{name}.bin'); subprocess.run(['./r34_pairs', f'rec_{name}.bin', repr(L), f'pairs_{name}.json'], check=True)
    return json.load(open(f'pairs_{name}.json'))['Ek']
def sign_check(mu):
    bad = tot = 0; ex = None; idx = np.nonzero(mu[1:5001])[0] + 1
    for m in idx[:150]:
        for nn in idx[idx > m]:
            g = math.gcd(int(m), int(nn)); u, v = m // g, nn // g; tot += 1
            if np.sign(mu[m]*mu[nn]) != np.sign(mu[u]*mu[v]):
                bad += 1
                if ex is None: ex = dict(m=int(m), n=int(nn), u=int(u), v=int(v), mu_m_mu_n=float(mu[m]*mu[nn]), mu_u_mu_v=float(mu[u]*mu[v]))
    return dict(pairs_tested=tot, sign_mismatch=bad, first_example=ex)
n = np.arange(1, X + 1); out = {}
models = {'ZETA': (mult_mu(lambda p: -1.0, True), K['ZETA']['slope']),
          'ZETA2': (mult_mu(lambda p: -2.0), K['FTH0.0']['slope']),
          'FTH0.1': (mult_mu(lambda p: -(p**0.1 + p**-0.1)), K['FTH0.1']['slope'])}
bDH = np.array([[0, 1, KAP, -KAP, -1][k % 5] for k in range(X + 1)], float); models['DH'] = (inverse_fast(bDH), K['DH']['slope'])
bE = np.zeros(X + 1); M = int(math.isqrt(X)) + 1
for x in range(-M, M + 1):
    y = np.arange(-M, M + 1); v = x*x + 5*y*y; v = v[(v > 0) & (v <= X)]; np.add.at(bE, v, 1)
models['EPS'] = (inverse_fast(bE/2), K['EPS']['slope'])
for name, (mu, c) in models.items():
    a = mu[1:]/np.sqrt(n); a[0] -= c; keep = np.nonzero(a)[0]
    E = run(name, records(a[keep], [facs[i + 1] for i in keep], np.log(n[keep].astype(float))))
    out[name] = dict(c=c, Ek=E, support=int(len(keep)), sign_check=sign_check(mu))
    print(name, 'E0=%.2f E1=%.2f R=%.2f G=%.4f' % (E[0], E[1], sum(E[1:]), sum(E)), out[name]['sign_check'], flush=True)
P = np.fromfile('../round32_mollified_floor/beur_primes.bin'); P = P[P <= X]
gl, gs, gf = [], [], []
def dfs(start, g, sgn, fs):
    gl.append(math.log(g)); gs.append(sgn); gf.append(list(fs))
    for i in range(start, len(P)):
        if g*P[i] > X: break
        dfs(i + 1, g*P[i], -sgn, fs + [i])
dfs(0, 1.0, 1, [])
cB = BZ['slope']; a = np.array(gs, float)*np.exp(-np.array(gl)/2); a[0] -= cB
E = run('BEUR', records(a, gf, np.array(gl)))
out['BEUR'] = dict(c=cB, Ek=E, support=len(gl), sign_check='holds by construction (multiplicative, squarefree support)')
print('BEUR', len(gl), 'E0=%.2f E1=%.2f R=%.2f G=%.4f' % (E[0], E[1], sum(E[1:]), sum(E)))
json.dump(dict(X=X, L=L, models=out), open('r34_controls.json', 'w'), indent=1)
