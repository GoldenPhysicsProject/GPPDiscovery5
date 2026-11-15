"""Round 7 GAP (a) certificate, well-conditioned form. NUMERICAL (interval-certified via arb; not Lean).
The direct ball LDL^T in the phi-basis loses certification to radius growth (entries ~1e-9, pivots ~1e-11; see r7_certify_results.json).
Fix (same span, exact change of basis): span{phi_13..phi_72, delta_l} = span{psi_13..psi_72, delta_l}, psi_j = phi_j - w_j delta_l.
Scaled basis: psi_j / sqrt(D_j) (D_j = 3 w(omega_j)/omega_j^4), delta_l / sqrt(g_mid), phi_i / sqrt(1e-8).
Certificate: A = S Gram S (ball matrix) = Am + E, |E_ab| <= R_ab. Let r_a = sum_b R_ab. Then x^T E x >= -sum_a r_a x_a^2, so
  A - beta S P S >= Am - beta S P S - diag(r),  and the latter is certified > 0 by LDL^T of an exact-midpoint matrix at 600 bits.
P = projector onto the E_12 coordinates. Positive pivots => Schur complement onto E_12 of Gram(E_12 ⊕ Y ⊕ delta) >= beta I (L^2 on [0,1/4]),
and Gram on Y ⊕ delta is positive definite."""
from flint import arb, acb, ctx
import json, glob
import mpmath as mp
ctx.prec = 160; mp.mp.dps = 50
NT = 72; N = 12
E = {}
for f in sorted(glob.glob('shard_*.json')): E.update(json.load(open(f)))
l = arb(1)/4; pi = arb.pi(); c2 = (2/l).sqrt()
om = lambda j: j*pi/l
w = lambda j: c2*(-1)**(j+1)/om(j)
def ent(a, b):
    if a == 'd' and b == 'd': return arb(E['d,d'])
    if b == 'd': return arb(E[f"{a-1},d"])
    if a == 'd': return arb(E[f"{b-1},d"])
    a, b = min(a, b), max(a, b); return arb(E[f"{a-1},{b-1}"])
def wsym(x): return ((acb(1)/4 + acb(0,1)*x/2).digamma()).real - pi.log()
Y = list(range(N+1, NT+1)); X = list(range(1, N+1))
g = ent('d','d')
# basis vectors as (kind, index): coefficients in phi/delta
basis = [('psi', j) for j in Y] + [('d', None)] + [('phi', i) for i in X]
def comps(b):
    k, j = b
    if k == 'psi': return [(j, arb(1)), ('d', -w(j))]
    if k == 'd': return [('d', arb(1))]
    return [(j, arb(1))]
scale = []
for (k, j) in basis:
    if k == 'psi': scale.append(1/(3*wsym(om(j))/om(j)**4).sqrt())
    elif k == 'd': scale.append(1/arb(g.mid()).sqrt())
    else: scale.append(1/arb('1e-8').sqrt())
n = len(basis); m = len(Y) + 1
A = [[None]*n for _ in range(n)]
for p in range(n):
    for q in range(p, n):
        v = sum((cp*cq*ent(a, b) for a, cp in comps(basis[p]) for b, cq in comps(basis[q])), arb(0))
        A[p][q] = A[q][p] = v*scale[p]*scale[q]
rowr = [sum(float(A[p][q].rad()) for q in range(n))*1.000001 + 1e-30 for p in range(n)]
rho = max(rowr)
import numpy as _np; _rr=_np.array([[float(A[p][q].rad()) for q in range(n)] for p in range(n)]); _i=_np.unravel_index(_rr.argmax(), _rr.shape); print('max rad at', _i, basis[_i[0]], basis[_i[1]], _rr.max(), 'rowsums delta,first psi,last psi,phi1', rowr[m-1], rowr[0], rowr[m-2], rowr[m])
Am = mp.matrix([[mp.mpf(A[p][q].mid().str(45, radius=False)) for q in range(n)] for p in range(n)])
# float Schur onto E_12 in scaled units -> beta in L^2 units = lam * 1e-8
A22 = Am[0:m, 0:m]; A12 = Am[m:n, 0:m]; A11 = Am[m:n, m:n]
lam = sorted(mp.eigsy(A11 - A12*mp.inverse(A22)*A12.T, eigvals_only=True))[0] * mp.mpf('1e-8')
lamD = sorted(mp.eigsy(A22, eigvals_only=True))[0]
ctx.prec = 600
mids = [[arb(A[p][q].mid()) for q in range(n)] for p in range(n)]
def ldl(beta):
    bs = beta/arb('1e-8')
    B = [[mids[p][q] - ((bs if p >= m else 0) + arb(rowr[p]) if p == q else 0) for q in range(n)] for p in range(n)]
    L = [[arb(0)]*n for _ in range(n)]; D = [arb(0)]*n
    for j in range(n):
        D[j] = B[j][j] - sum((L[j][k]**2*D[k] for k in range(j)), arb(0))
        if not D[j] > 0: return None, j
        for i in range(j+1, n): L[i][j] = (B[i][j] - sum((L[i][k]*L[j][k]*D[k] for k in range(j)), arb(0)))/D[j]
    return D, None
out = dict(label="NUMERICAL (interval-certified; not Lean)", NT=NT, N=N, basis="psi_13..psi_72 (D-scaled), delta_l, phi_1..phi_12",
           rho_perturbation_bound_max_rowsum=rho, perturbation='E <= diag(row sums of R) (|x_a x_b| <= (x_a^2+x_b^2)/2)', max_scaled_entry_radius=max(float(A[p][q].rad()) for p in range(n) for q in range(n)),
           lam_min_schur_float=mp.nstr(lam, 12), lam_min_scaled_Y_delta_block_float=mp.nstr(lamD, 10), tried=[])
for fac in ['0.995', '0.99', '0.98', '0.95', '0.9']:
    beta = arb(mp.nstr(lam*mp.mpf(fac), 6)); D, bad = ldl(beta)
    out['tried'].append(dict(factor=fac, beta=beta.str(8), certified=D is not None, fail_at=bad))
    if D is not None:
        out.update(beta1_certified=beta.str(10), pivots_Y_delta_min_lower=min(D[j].lower() for j in range(m)).str(6),
                   pivot_delta_lower=D[m-1].lower().str(6), pivots_E12_lower_scaled=[D[j].lower().str(6) for j in range(m, n)],
                   pivots_E12_lower_L2units=[(D[j].lower()*arb('1e-8')).str(6) for j in range(m, n)])
        break
print(json.dumps(out, indent=1, default=str)); json.dump(out, open('r7_certify2_results.json','w'), indent=1, default=str)
