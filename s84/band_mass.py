"""S84: where does ||f||^2 sit in frequency? (1/2pi) int |fhat|^2 dt = ||f||^2 (Plancherel). Fraction in |t|<X; mode-mass profile;
nodal count of f on [0,L] at high precision (S82-OS float count may be rounding-polluted). NUMERICAL (mpmath)."""
import json, sys, mpmath as mp
from spec_decomp import fhat_factory
import os; mp.mp.dps = int(os.environ.get("DPS","120"))
def run(fn, sec, Xs):
    d = json.load(open(fn)); L = d['L'].replace('_','/')
    f, n2, Lm, v, fr, moms = fhat_factory(d['vec'], L, sec)
    Tstar = 2*mp.pi*mp.e**(2*Lm)
    out = {"file": fn.split('/')[-1], "L": L, "T*": float(Tstar)}
    pts = [0] + [mp.mpf(x) for x in range(10, 410, 10)]
    seg = [mp.quad(lambda t: f(t)**2, [pts[i], pts[i+1]]) for i in range(len(pts)-1)]
    tot = mp.pi*n2
    for X in Xs + [float(Tstar)]:
        out["frac_|t|<%.0f" % X] = mp.nstr(mp.fsum(s for i, s in enumerate(seg) if pts[i+1] <= X)/tot, 6)
    out["frac_|t|<400"] = mp.nstr(mp.fsum(seg)/tot, 8)
    w = [x**2 for x in v]; W = mp.fsum(w)
    out["mode_mass_top10"] = mp.nstr(mp.fsum(w[-10:])/W, 4); out["mode_mass_top25pct"] = mp.nstr(mp.fsum(w[3*len(w)//4:])/W, 4)
    out["mean_mode_freq"] = mp.nstr(mp.fsum(w[m]*fr[m] for m in range(len(v)))/W, 5)
    # nodal count on [0,L]
    N = 3000; xs = [Lm*i/N for i in range(N)]  # exclude u=L (all modes vanish there in even case)
    if sec == 'even': g = lambda u: mp.fsum(v[m]*mp.cos(fr[m]*u) for m in range(len(v)))
    else: g = lambda u: mp.fsum(v[m]*mp.sin(fr[m]*u) for m in range(len(v)))
    vals = [g(x) for x in xs[1:]]
    out["sign_changes_(0,L)_mp"] = sum(1 for i in range(len(vals)-1) if vals[i]*vals[i+1] < 0)
    out["max|f|/rms"] = mp.nstr(max(abs(x) for x in vals)/mp.sqrt(n2/(2*Lm)), 4)
    return out
if __name__ == '__main__':
    print(json.dumps(run(sys.argv[1], sys.argv[2], [50, 100, 150, 236, 300, 500]), indent=1))
