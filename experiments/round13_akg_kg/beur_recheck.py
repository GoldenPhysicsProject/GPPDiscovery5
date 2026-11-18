"""Re-run BEUR rows (harness.operators/field) after the prime-placement fix; also report psi_B - psi_C and field vs quartet."""
import numpy as np, math, json
import harness as h
from akg_kg import make_g
B = h.beurling_model(); out = dict(nprimes=B['nprimes'])
q = B['quartet']
def psiC(X, x0=20.0):
    return (X - x0) - sum((X**complex(r) - x0**complex(r))/complex(r) for r in q).real
for X in (1e3, 1e4, 1e5, 1e6):
    out[f'psiB_minus_psiC_{X:.0e}'] = float(np.sum(B['c'][B['n'] < X]) - psiC(X))
def G(g, L, z):
    y = np.linspace(-L/2, L/2, 4001); return np.trapezoid(g(y)*np.exp(z*y), y)
rows = []
for L in (0.5, 1.0, 2.0, 3.0, 4.0):
    M = int(max(24, 10*L)); arch, pole, prime, absprime = h.operators(B, L, M)
    row = dict(L=L, lam_min_W_no_arch=float(np.linalg.eigvalsh(pole-prime)[0]),
               lam_max_prime=float(np.linalg.eigvalsh(prime)[-1]) if np.any(prime) else 0.0,
               lam_max_prime_minus_pole=float(np.linalg.eigvalsh(prime-pole)[-1]),
               lam_min_pole=float(np.linalg.eigvalsh(pole)[0]))
    # same atoms with zeta's arch borrowed (B has no FE, so this arch is not its own): does arch rescue it?
    Z = h.zeta_model(); az, _, _, _ = h.operators(dict(Z, n=np.array([2.0]), c=np.array([0.0])), L, M)
    row['lam_min_with_zeta_arch'] = float(np.linalg.eigvalsh(az+pole-prime)[0])
    g = make_g(L, 'even'); u0 = math.log(20)+L/2+0.05; us = np.arange(u0, 13.0, 0.004); F = h.field(B, g, L, us)
    Fq = np.zeros(len(us))
    for r in q: Fq -= (G(g, L, complex(r)-0.5)*np.exp((complex(r)-0.5)*us)).real
    row['rms_field_by_u'] = {str(a): float(np.sqrt(np.mean(F[(us>=a)&(us<a+1)]**2))) for a in (6.0, 9.0, 12.0)}
    row['rms_quartet_by_u'] = {str(a): float(np.sqrt(np.mean(Fq[(us>=a)&(us<a+1)]**2))) for a in (6.0, 9.0, 12.0)}
    row['rms_resid_by_u'] = {str(a): float(np.sqrt(np.mean((F-Fq)[(us>=a)&(us<a+1)]**2))) for a in (6.0, 9.0, 12.0)}
    rows.append(row); print(json.dumps(row), flush=True)
out['rows'] = rows; print(json.dumps({k: v for k, v in out.items() if k != 'rows'}))
json.dump(out, open('beur_recheck.json', 'w'), indent=1)
