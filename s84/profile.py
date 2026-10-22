"""S84: fhat profile + low-band mass fractions for the near-null / witness vectors. NUMERICAL (mpmath)."""
import json, mpmath as mp, os
os.environ.setdefault("DPS", "220")
from spec_decomp import fhat_factory
mp.mp.dps = 220
out = []
for fn in ["hp/cosvec_ZETA_even_L71_50_K128.json", "dh/cosvec_DH_even_L71_50_K96.json", "hp2/cosvec_ZETA_even_L137_80_K160.json", "../s75/cosvec_DH_even_L137_80_K129.json"]:
    d = json.load(open(fn)); f, n2, Lm, v, fr, m = fhat_factory(d['vec'], d['L'].replace('_','/'), 'even')
    tot = mp.pi*n2; r = {"file": fn, "L": d['L']}
    for X in [2, 5, 10, 20]: r["frac|t|<%d" % X] = mp.nstr(mp.quad(lambda t: f(t)**2, mp.linspace(0, X, 9))/tot, 6)
    r["|fhat(t)|/|f|"] = {str(t): mp.nstr(abs(f(mp.mpf(t)))/mp.sqrt(n2), 3) for t in [0, 2, 5, 10, "14.134725141734693790", 30, 50, 100]}
    r["fhat(i/2)/|f|"] = mp.nstr(mp.re(f(mp.mpc(0, 0.5)))/mp.sqrt(n2), 6)
    out.append(r); print(json.dumps(r))
json.dump(out, open("profile_results.json", "w"), indent=1)
