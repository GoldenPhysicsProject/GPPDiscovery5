"""Spliced DMV-type proxy B'': the real prime powers p^k with p < 55, plus discrete Beurling primes >= 55 tracking
dPsi_C = 1_{x>=55}(1 - sum_q x^{rho_q-1})dx (quartet rho_q = 1/2 +- 0.1 +- 28i). The Euler factors (1-p^-s)^-1, p<55,
have no zeros or poles in Re s > 0, so zeta_{B''} has exactly the quartet as zeros in Re s > 0 (Pi_B - Pi_C bounded).
Since e^4 < 55, every self-block operator at width L <= 4 (arch, pole, prime, signed prime-pole, W) built on B'' with
zeta's arch is IDENTICAL to zeta's. The field at large u is where B'' differs."""
import numpy as np, math, json
import harness as h
from akg_kg import make_g, hann_mean, vonmangoldt_table
Z = h.zeta_model(); B = h.beurling_model(x0=55.0)
m = Z['n'] < 55
S = dict(Z, name='SPLICE', n=np.concatenate([Z['n'][m], B['n']]), c=np.concatenate([Z['c'][m], B['c']]))
o = np.argsort(S['n']); S['n'] = S['n'][o]; S['c'] = S['c'][o]
out = dict(nprimes_beurling=B['nprimes'])
def G(g, L, z):
    y = np.linspace(-L/2, L/2, 4001); return np.trapezoid(g(y)*np.exp(z*y), y)
for L in (0.5, 1.0, 2.0, 3.0, 4.0):
    M = int(max(24, 10*L)); oz = h.operators(Z, L, M); os_ = h.operators(S, L, M)
    diff = max(float(np.max(np.abs(a-b))) for a, b in zip(oz, os_))
    g = make_g(L, 'even'); u0 = math.log(55)+L/2+0.05; us = np.arange(u0, 13.0, 0.004)
    Fs = h.field(S, g, L, us); Fz = h.field(Z, g, L, us); W = h.W_scalar(Z, g, L)
    Fq = np.zeros(len(us))
    for r in B['quartet']: Fq += (G(g, L, complex(r)-0.5)*np.exp((complex(r)-0.5)*us)).real
    row = dict(max_abs_operator_difference_vs_zeta=diff,
               hann_over_W_zeta={str(U): hann_mean(us, Fz**2, u0, U)/W for U in (9.0, 11.0, 13.0)},
               hann_over_W_splice={str(U): hann_mean(us, Fs**2, u0, U)/W for U in (9.0, 11.0, 13.0)},
               rms_splice_by_u={str(a): float(np.sqrt(np.mean(Fs[(us>=a)&(us<a+1)]**2))) for a in (6.0, 9.0, 12.0)},
               rms_quartet_by_u={str(a): float(np.sqrt(np.mean(Fq[(us>=a)&(us<a+1)]**2))) for a in (6.0, 9.0, 12.0)},
               rms_splice_minus_quartet_by_u={str(a): float(np.sqrt(np.mean((Fs-Fq)[(us>=a)&(us<a+1)]**2))) for a in (6.0, 9.0, 12.0)})
    out[str(L)] = row; print(L, json.dumps(row), flush=True)
json.dump(out, open('beur_spliced.json', 'w'), indent=1)
