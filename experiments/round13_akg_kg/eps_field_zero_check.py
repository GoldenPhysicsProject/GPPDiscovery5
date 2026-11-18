"""EPS field vs zero sum with the complete off-line list (epstein_zeros_lowband.json, argument-principle count 7 in
[0.501,1.6]x[0.5,60.5]) and on-line zeros t<60 from epstein_zeros.json. Test data only."""
import numpy as np, json, math
import harness as h
from akg_kg import make_g
on = json.load(open('epstein_zeros.json'))['online']; off = json.load(open('epstein_zeros_lowband.json'))['offline_right_half']
E = h.epstein_model(); out = {}
def G(g, L, z):
    y = np.linspace(-L/2, L/2, 4001); return np.trapezoid(g(y)*np.exp(z*y), y)
for L in (2.0, 3.0, 4.0):
    g = make_g(L, 'even'); us = np.arange(max(2.0, L/2+0.1), 13.0, 0.004); F = h.field(E, g, L, us)
    zs = [s*1j*t for t in on for s in (1, -1)] + [complex(sb*(b-0.5), st*t) for b, t in off for sb in (1, -1) for st in (1, -1)]
    Fz = sum((G(g, L, z)*np.exp(z*us)).real for z in zs)
    Foff = sum((G(g, L, complex(b-0.5, st*t))*np.exp(complex(b-0.5, st*t)*us)).real for b, t in off for st in (1, -1))
    sel = us > 9
    # u at which the off-line part would reach the size of the on-line part (extrapolating the dominant zero)
    b, t = max(off, key=lambda z: abs(G(g, L, complex(z[0]-0.5, z[1])))*math.exp((z[0]-0.5)*13))
    out[str(L)] = dict(rel_err=float(np.sqrt(np.mean((F-Fz)[sel]**2)/np.mean(F[sel]**2))),
                       offline_rms_share={str(a): float(np.sqrt(np.mean(Foff[(us>=a)&(us<a+1)]**2)/np.mean(F[(us>=a)&(us<a+1)]**2))) for a in (3.0, 6.0, 9.0, 12.0)},
                       dominant_offline_zero=[b, t])
    print(L, out[str(L)], flush=True)
json.dump(out, open('eps_field_zero_check.json', 'w'), indent=1)
