"""Clue probe C: Lagarias-type positivity Re xi'/xi(s) > 0 for Re s > 1/2 (equivalent to RH for zeta; Lagarias,
Acta Arith. 89 (1999)). Sample min_t Re(Lam'/Lam)(sigma+it) for completed zeta, EPS (x^2+5y^2) and DH, sigma in
{0.55,...,1.5}, t in [2, 120]. Where is positivity lost, and for which model? Zeros are not inputs."""
import mpmath as mp, json, numpy as np
mp.mp.dps = 20
chi_m4 = [0, 1, 0, -1]; chi_5 = [0, 1, -1, -1, 1]
chi_m20 = [chi_m4[n % 4]*chi_5[n % 5] for n in range(20)]
chi_c = [0, 1, 1j, -1j, -1]; chi_cb = [0, 1, -1j, 1j, -1]
kap = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1); al, be = (1-1j*kap)/2, (1+1j*kap)/2
def lam_zeta(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def lam_eps(s): return s*(s-1)*(mp.sqrt(20)/(2*mp.pi))**s*mp.gamma(s)*(mp.zeta(s)*mp.dirichlet(s, chi_m20)+mp.dirichlet(s, chi_m4)*mp.dirichlet(s, chi_5))
def lam_dh(s): return (5/mp.pi)**(s/2)*mp.gamma((s+1)/2)*(al*mp.dirichlet(s, chi_c)+be*mp.dirichlet(s, chi_cb))
def logder(f, s, h=mp.mpf('1e-6')): return (f(s+h)-f(s-h))/(2*h)/f(s)
sig = [0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5]
ts = np.arange(2.0, 120.0, 0.25)
out = {}
for name, f in [('ZETA', lam_zeta), ('EPS', lam_eps), ('DH', lam_dh)]:
    row = {}
    for s0 in sig:
        v = [float(mp.re(logder(f, mp.mpc(s0, t)))) for t in ts]
        i = int(np.argmin(v)); neg = [float(ts[k]) for k in range(len(v)) if v[k] < 0]
        row[str(s0)] = dict(min=v[i], t_at_min=float(ts[i]), n_negative=len(neg), neg_t_ranges=[neg[0], neg[-1]] if neg else None)
        print(name, s0, row[str(s0)], flush=True)
    out[name] = row
json.dump(out, open('lagarias_pos.json', 'w'), indent=1)
