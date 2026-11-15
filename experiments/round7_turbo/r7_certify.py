"""Round 7 GAP (a) certificate: ball LDL^T of the bordered Gram on {phi_1..phi_72, delta_l} minus diag(beta1 I_12, 0).
Elimination order: Y = phi_13..phi_72, then delta_l, then E_12. Positive pivots on Y+delta certify the Gram on Y ⊕ {delta} is positive
definite (corrected P2 there); positive pivots on E_12 certify Schur complement onto E_12 minus beta1 I is positive definite.
NUMERICAL (interval-certified via arb, 160-bit); not Lean."""
from flint import arb, ctx
import json, glob
import mpmath as mp
ctx.prec = 160; mp.mp.dps = 50
NT = 72; N = 12
E = {}
for f in sorted(glob.glob('shard_*.json')): E.update(json.load(open(f)))
assert len(E) == NT*(NT+1)//2 + NT + 1, len(E)
def ent(a, b):
    if a == 'd' and b == 'd': return arb(E['d,d'])
    if b == 'd': return arb(E[f"{a},d"])
    if a == 'd': return arb(E[f"{b},d"])
    a, b = min(a, b), max(a, b); return arb(E[f"{a},{b}"])
order = list(range(N, NT)) + ['d'] + list(range(N))
n = len(order)
M = [[ent(order[p], order[q]) for q in range(n)] for p in range(n)]
maxrad = max(float(M[p][q].rad()) for p in range(n) for q in range(n))
# consistency vs Round 5 block and Round 6 g,u
R5 = json.load(open('../round5_turbo/p1_oddblock_L05_N12_results.json')); R6 = json.load(open('../round6_turbo/p2_bordered_arb_results.json'))
d5 = max(abs(float(ent(j,k).mid()) - float(R5['M'][j][k])) for j in range(N) for k in range(N))
d6g = abs(float(ent('d','d').mid()) - float(arb(R6['g']).mid())); d6u = max(abs(float(ent(i,'d').mid()) - float(arb(R6['u'][i]).mid())) for i in range(N))
# float Schur from midpoints
Mm = mp.matrix([[mp.mpf(M[p][q].mid().str(45, radius=False)) for q in range(n)] for p in range(n)])
m = NT - N + 1
A22 = Mm[0:m, 0:m]; A12 = Mm[m:n, 0:m]; A11 = Mm[m:n, m:n]
Sc = A11 - A12 * mp.inverse(A22) * A12.T
lam = sorted(mp.eigsy(Sc, eigvals_only=True))[0]
lamD = sorted(mp.eigsy(A22, eigvals_only=True))[0]
def ldl(beta):
    A = [[M[p][q] - (beta if (p == q and p >= m) else 0) for q in range(n)] for p in range(n)]
    L = [[arb(0)]*n for _ in range(n)]; D = [arb(0)]*n
    for j in range(n):
        D[j] = A[j][j] - sum((L[j][k]**2*D[k] for k in range(j)), arb(0))
        if not D[j] > 0: return None, j, D[j]
        for i in range(j+1, n):
            L[i][j] = (A[i][j] - sum((L[i][k]*L[j][k]*D[k] for k in range(j)), arb(0)))/D[j]
    return D, None, None
best = None; tried = []
for fac in ['0.99', '0.98', '0.95', '0.9']:
    beta = arb(mp.nstr(lam*mp.mpf(fac), 6))
    D, bad, val = ldl(beta)
    tried.append(dict(factor=fac, beta=beta.str(8), certified=D is not None, fail_at=bad, fail_val=None if val is None else val.str(5)))
    print(tried[-1], flush=True)
    if D is not None: best = (beta, D); break
out = dict(label="NUMERICAL (interval-certified; not Lean)", NT=NT, N=N, order="Y(phi_13..phi_72), delta_l, E_12", prec_bits=160,
           max_entry_radius=maxrad, check_vs_round5_M11_maxabs=d5, check_vs_round6_g=d6g, check_vs_round6_u_maxabs=d6u,
           lam_min_schur_float=mp.nstr(lam, 12), lam_min_Gram_Y_delta_float=mp.nstr(lamD, 12), tried=tried)
if best:
    beta, D = best
    out.update(beta1_certified=beta.str(10), pivots_Y_delta_lower_min=min(D[j].lower() for j in range(m)).str(6),
               pivot_delta=D[m-1].str(8), pivots_E12_lower=[D[j].lower().str(6) for j in range(m, n)],
               pivots_Y_lower_first_last=[D[0].lower().str(6), D[m-2].lower().str(6)],
               max_pivot_radius=max(float(D[j].rad()) for j in range(n)))
print(json.dumps(out, indent=1, default=str)); json.dump(out, open('r7_certify_results.json','w'), indent=1, default=str)
