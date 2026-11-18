"""Heart arm: Woodbury (D2 Result D) and parity split (D2 sec.10) checked symbolically/exactly.
(1) generic 3x3 symbolic G (symmetric), symbolic U columns: B^{-1} = G^{-1} - G^{-1}U Q^{-1} U^T G^{-1}, Q=J^{-1}+U^T G^{-1} U.
(2) 200 random exact-rational instances in dim 6 incl. moment Gram M = Y^T B^{-1} Y = M_G - V^T Q^{-1} V.
(3) parity: G commuting with reflection R, Rc=c, Rs=-s => <c,G^{-1}s>=0, Q diagonal, B>0 <=> 2 beta<1 (G>0), rank-one formula.
"""
import sympy as s, random
random.seed(13)
g = s.symbols('g0:6'); u = s.symbols('u0:6')
G = s.Matrix([[g[0],g1:=g[1],g[2]],[g[1],g[3],g[4]],[g[2],g[4],g[5]]])
U = s.Matrix([[u[0],u[3]],[u[1],u[4]],[u[2],u[5]]]); J = s.diag(2,-2)
B = G+U*J*U.T; Gi = G.inv(); Q = J.inv()+U.T*Gi*U
D = (Gi - Gi*U*Q.inv()*U.T*Gi)*B - s.eye(3)
assert all(s.simplify(x) == 0 for x in D); print('PASS symbolic 3x3 Woodbury')
def rr(): return s.Rational(random.randint(-9,9), random.randint(1,5))
n_ok = 0
for _ in range(200):
    X = s.Matrix(6,6,lambda i,j: rr()); G = X*X.T + s.eye(6)
    U = s.Matrix(6,2,lambda i,j: rr()); B = G+U*J*U.T
    if B.det() == 0: continue
    Gi = G.inv(); Q = J.inv()+U.T*Gi*U
    assert Gi - Gi*U*Q.inv()*U.T*Gi == B.inv()
    Y = s.Matrix(6,3,lambda i,j: rr())
    assert Y.T*B.inv()*Y == Y.T*Gi*Y - (U.T*Gi*Y).T*Q.inv()*(U.T*Gi*Y)
    n_ok += 1
print('PASS exact-rational Woodbury + moment Gram,', n_ok, 'instances')
# parity: basis (x_-2..x_2) dim 5 with R reversing order
R = s.Matrix(5,5,lambda i,j: 1 if i+j==4 else 0); n_ok = 0; n_pos = 0
for _ in range(200):
    X = s.Matrix(5,5,lambda i,j: rr()); G0 = X*X.T + s.eye(5); G = (G0 + R*G0*R)/2
    a = [rr(),rr(),rr()]; c = s.Matrix([a[0],a[1],a[2],a[1],a[0]])
    b = [rr(),rr()]; sv = s.Matrix([b[0],b[1],0,-b[1],-b[0]])
    Gi = G.inv(); assert (c.T*Gi*sv)[0] == 0
    al = (c.T*Gi*c)[0]; be = (sv.T*Gi*sv)[0]
    B = G + 2*c*c.T - 2*sv*sv.T
    posdef = all(B[:k,:k].det() > 0 for k in range(1,6))
    assert posdef == (2*be < 1)
    if 2*be != 1:
        assert B.inv() == Gi - 2/(1+2*al)*(Gi*c)*(Gi*c).T + 2/(1-2*be)*(Gi*sv)*(Gi*sv).T
    n_ok += 1; n_pos += posdef
print('PASS parity split: cross term 0, B>0 <=> 2beta<1, rank-one inverse;', n_ok, 'instances,', n_pos, 'positive')
