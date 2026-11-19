"""Round 16 Part 1 ablation: remove one ingredient at a time from Bombieri's auxiliary function and see the bound collapse.
 A (no additivity): u^p -> u^e with e = p+1 (multiplicative, NOT additive). F still vanishes at rational points (S4 intact),
   but F is no longer an e-th power, so only multiplicity >= 1 is guaranteed: N <= 1 + pole(F)  (trivial: > q).
 B (Frobenius does not fix the points): v^q -> v^p (the F_p-Frobenius, additive AND multiplicative, but it does not fix
   F_q-points): F(P) != delta(F)(P), F vanishes only where v^p(P) = v(P) on the monomials, i.e. at F_p-points.
 C (Riemann-Roch slack): take m one smaller than needed: the kernel of delta is 0, no auxiliary function exists."""
import json, numpy as np
from ffq import Fq, Curve
from stepanov_ff import nullspace, coeffvec, choose_lm, count_Fp

def mpow(cur, i, j, e):
    F = cur.F; X = F.shiftx(F.P([1]), i*e); z = X[:, :0]; k = j*e
    return (F.mul(X, cur.fpow(k//2)), z) if k % 2 == 0 else (z, F.mul(X, cur.fpow(k//2)))

def run(cur, eu, ev, l, m):
    F = cur.F; p, q = F.p, F.q; Bl, Bm = cur.basis(l), cur.basis(m)
    top = max(2*i + j*cur.d for i, j in Bl)*eu + m
    cols, pairs = [], []
    for (i, j) in Bl:
        up = mpow(cur, i, j, eu)
        for (k, jj) in Bm: cols.append(coeffvec(cur, cur.mul(up, cur.mono(k, jj)), top)); pairs.append(((i, j), (k, jj)))
    ker, rank = nullspace(F, np.array(cols).T)
    if not ker: return dict(kernel_dim=0)
    c = ker[0]; Fn = (F.P([]), F.P([]))
    for cc, ((i, j), (k, jj)) in zip(c, pairs):
        if cc: Fn = cur.add(Fn, cur.smul(int(cc), cur.mul(mpow(cur, i, j, eu), mpow(cur, k, jj, ev))))
    pts = cur.points(); nz = sum(1 for x, y in pts if cur.ev(Fn, x, y) == 0)
    return dict(kernel_dim=len(ker), pole_F=cur.pole(Fn), rational_points=len(pts) + 1, F_vanishes_at=nz)

out = []
for kind, p, f in [('E', 7, [3, 1, 0, 1]), ('E', 11, [1, 1, 0, 1]), ('G2', 11, [1, 1, 0, 2, 0, 1])]:
    F = Fq(p); C = Curve(F, f); l, m, b = choose_lm(C); q = F.q; N = len(C.points()) + 1
    full = run(C, p, q, l, m); A = run(C, p + 1, q, l, m); B = run(C, p, p, l, m); Cc = run(C, p, q, l, m - 2)
    rec = dict(kind=kind, p=p, q=q, N=N, NFp=count_Fp(C), l=l, m=m,
               full=dict(full, bound=1 + full['pole_F']//p),
               A_no_additivity=dict(A, bound=1 + A['pole_F'] if A.get('pole_F') else None),
               B_frob_not_fixing=B, C_RR_short=Cc)
    out.append(rec); print(rec, flush=True)
json.dump(out, open('ablation_ff.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
