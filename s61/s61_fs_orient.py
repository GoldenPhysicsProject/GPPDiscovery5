"""Grok S61 H61a(1): CP1 antipode R z = -1/conj(z) on the affine chart, symbolic (sympy).
Checks det DR, DR J = -J DR, R^* g_FS = g_FS, R^* dmu_FS = dmu_FS (|det|), R^* omega_FS = -omega_FS (signed det).
FS metric on the chart: g = rho(r) (dx^2+dy^2), rho = (1+r)^-2 (constant factor irrelevant), omega = rho dx^dy. CLASSICAL."""
import sympy as sp, json
x, y = sp.symbols('x y', real=True); r = x**2 + y**2
X, Y = -x/r, -y/r                              # -1/conj(z) = -z/|z|^2
DR = sp.Matrix([[sp.diff(X, x), sp.diff(X, y)], [sp.diff(Y, x), sp.diff(Y, y)]])
J = sp.Matrix([[0, -1], [1, 0]])
rho = lambda rr: 1/(1+rr)**2
rR = sp.simplify(X**2 + Y**2)
det = sp.simplify(DR.det())
out = dict(
  R_r=str(rR), detDR=str(det), detDR_eq_minus_r_m2=sp.simplify(det + r**-2) == 0,
  DRJ_plus_JDR_zero=sp.simplify(DR*J + J*DR) == sp.zeros(2, 2),
  metric_pullback_ok=sp.simplify(DR.T*(rho(rR)*sp.eye(2))*DR - rho(r)*sp.eye(2)) == sp.zeros(2, 2),
  measure_pullback_abs_ok=sp.simplify(rho(rR)*(-det) - rho(r)) == 0,      # |det| = -det since det<0
  omega_pullback_signed=str(sp.simplify(rho(rR)*det/rho(r))),             # -1 => R^*omega = -omega
  involution_ok=sp.simplify((-X/(X**2+Y**2)) - x) == 0 and sp.simplify((-Y/(X**2+Y**2)) - y) == 0,
  fixed_points="none on CP1 (|z|^2 = -1 impossible): free antiholomorphic involution, CLASSICAL")
# Kahler compatibility g(u,v) = omega(u, J v) is preserved: both sides flip sign under R
print(json.dumps(out, indent=1)); json.dump(out, open('fs_orient.json', 'w'), indent=1)
