"""Grok S35: EMPIRICAL (not a proof) ulp error of numpy float64 cos/sin on this box vs Arb (python-flint, 200 bits),
on random phases in the range used by the torus B&B (|phase| <= |K|_1 * pi ~ 70). Informational for the libm assumption."""
import numpy as np, math
from flint import arb, ctx
ctx.prec = 200
rs = np.random.default_rng(1); x = np.concatenate([rs.uniform(-70, 70, 200000), rs.uniform(-1e-3, 1e-3, 20000) + np.pi / 2 * rs.integers(-40, 40, 20000)])
for name, f, af in (('cos', np.cos, lambda a: a.cos()), ('sin', np.sin, lambda a: a.sin())):
    y = f(x); worst = 0.0; worst_abs = 0.0
    for xi, yi in zip(x, y):
        t = af(arb(float(xi))); err = abs(arb(float(yi)) - t); e = float(err.upper())
        ulp = math.ulp(abs(float(t.mid()))) if float(t.mid()) != 0 else 5e-324
        worst = max(worst, e / ulp); worst_abs = max(worst_abs, e)
    print(f'numpy {np.__version__} {name}: max error {worst:.3f} ulp (max abs {worst_abs:.3e}) over {len(x)} samples')
