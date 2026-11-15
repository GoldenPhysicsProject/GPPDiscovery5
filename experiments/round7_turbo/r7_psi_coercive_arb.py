"""Round 7 (T1) finite section: certify Gram of psi_j = phi_j - w_j delta_l (j = 13..72) >= c1 * diag(D_j), D_j = 3 w(omega_j)/omega_j^4,
from the certified entries (exact linear transform in ball arithmetic). NUMERICAL (interval-certified; not Lean)."""
from flint import arb, acb, ctx
import json, glob
import mpmath as mp
ctx.prec = 160; mp.mp.dps = 50
NT = 72; N = 12
E = {}
for f in sorted(glob.glob('shard_*.json')): E.update(json.load(open(f)))
l = arb(1)/4; pi = arb.pi(); c2 = (2/l).sqrt()
om = lambda j: j*pi/l
w = lambda j: c2*(-1)**(j+1)/om(j)               # 1-indexed
def ent(a, b):
    if a == 'd' and b == 'd': return arb(E['d,d'])
    if b == 'd': return arb(E[f"{a-1},d"])
    if a == 'd': return arb(E[f"{b-1},d"])
    a, b = min(a, b), max(a, b); return arb(E[f"{a-1},{b-1}"])
J = list(range(N+1, NT+1)); n = len(J)
g = ent('d','d')
G = [[ent(j,k) - w(k)*ent(j,'d') - w(j)*ent(k,'d') + w(j)*w(k)*g for k in J] for j in J]
def wsym(x): return ((acb(1)/4 + acb(0,1)*x/2).digamma()).real - pi.log()
D = [3*wsym(om(j))/om(j)**4 for j in J]
Gn = mp.matrix([[mp.mpf((G[a][b]/(D[a]*D[b]).sqrt()).mid().str(40, radius=False)) for b in range(n)] for a in range(n)])
lam = sorted(mp.eigsy(Gn, eigvals_only=True))[0]
def ldl(c):
    A = [[G[a][b] - (c*D[a] if a == b else 0) for b in range(n)] for a in range(n)]
    L = [[arb(0)]*n for _ in range(n)]; P = [arb(0)]*n
    for j in range(n):
        P[j] = A[j][j] - sum((L[j][k]**2*P[k] for k in range(j)), arb(0))
        if not P[j] > 0: return None
        for i in range(j+1, n): L[i][j] = (A[i][j] - sum((L[i][k]*L[j][k]*P[k] for k in range(j)), arb(0)))/P[j]
    return P
res = dict(label="NUMERICAL (interval-certified; not Lean)", range="psi_13..psi_72", lam_min_Dnorm_float=mp.nstr(lam, 10))
for fac in ['0.99', '0.95', '0.9', '0.8']:
    c = arb(mp.nstr(lam*mp.mpf(fac), 5)); P = ldl(c)
    if P: res.update(c1_certified=c.str(8), min_pivot_lower=min(p.lower() for p in P).str(5)); break
print(json.dumps(res, indent=1)); json.dump(res, open('r7_psi_coercive_results.json','w'), indent=1)
