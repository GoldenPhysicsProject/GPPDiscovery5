"""Round 14E2: t-grid check for mpgap.build: trapezoid h = 0.1 vs h = 0.2 (every other node) on the 90-digit line data.
Raw family members do not satisfy the pole constraint, so F(1/2+it)H^(t) has F's pole at t = -i/2 and the trapezoid error is
~exp(-pi/h) (1e-7 at h = 0.2).  After projecting onto the pole-constrained combinations (B, as mpgap uses) the integrand is entire,
and the h = 0.1 / 0.2 difference must collapse; that is the number reported."""
import sys, mpmath as mp
import mpgap as G
model, L, F = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
a = G.build(model, L, F, 1); b = G.build(model, L, F, 2); B = G.nullspace_rows(a[4])
def mx(X): return max(abs(X[r, c]) for r in range(X.rows) for c in range(X.cols))
for name, i, f in [('Gram ||EH||^2', 2, lambda X: B.T*X*B), ('window coeffs <phi,EH>', 3, lambda X: X*B)]:
    print(model, L, F, name, 'raw members: max|h=0.1 - h=0.2| =', mp.nstr(mx(a[i]-b[i]), 4),
          '| pole-constrained: ', mp.nstr(mx(f(a[i])-f(b[i])), 4), ' max entry', mp.nstr(mx(f(a[i])), 4), flush=True)
