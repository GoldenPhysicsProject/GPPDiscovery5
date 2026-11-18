"""Round 13e checks: (1) L=4 field rows with u0 = L/2+0.1 (fixes w(0) NaN); (2) EPS field vs zero sum from
epstein_zeros.json (test data only), incl. off-line quartets rho, conj, 1-rho, 1-conj; (3) BEUR field vs quartet
prediction; (4) Galerkin M-convergence of lam_min(W_EPS) at L=4 and frequency content of the negative eigenvector."""
import numpy as np, math, json
import harness as h
from akg_kg import make_g, hann_mean, integ
from defect_kappa import basis
Z = json.load(open('epstein_zeros.json'))
out = {}
mods = {m['name']: m for m in (h.zeta_model(), h.epstein_model(), h.beurling_model())}

def G(g, L, z):  # int g(y) e^{z y} dy
    y = np.linspace(-L/2, L/2, 4001); return np.trapezoid(g(y)*np.exp(z*y), y)

for L in (2.0, 3.0, 4.0):
    g = make_g(L, 'even'); u0 = L/2+0.1; us = np.arange(max(u0, 2.0), 13.0, 0.004)
    # EPS zero-sum field
    F = h.field(mods['EPS'], g, L, us)
    zs = []
    for t in Z['online']: zs += [complex(0, t), complex(0, -t)]
    for b, t in Z['offline']:
        for zz in (complex(b-0.5, t), complex(b-0.5, -t), complex(0.5-b, t), complex(0.5-b, -t)): zs.append(zz)
    Fz = np.zeros(len(us))
    for zz in zs: Fz += (G(g, L, zz)*np.exp(zz*us)).real
    Foff = np.zeros(len(us))
    for b, t in Z['offline']:
        for zz in (complex(b-0.5, t), complex(b-0.5, -t)): Foff += (G(g, L, zz)*np.exp(zz*us)).real
    sel = us > 9
    out[f'EPS_L{L}'] = dict(rel_err_field_vs_zeros_u_gt_9=float(np.sqrt(np.mean((F-Fz)[sel]**2)/np.mean(F[sel]**2))),
        offline_share_rms_by_u={str(a): float(np.sqrt(np.mean(Foff[(us>=a)&(us<a+1)]**2)/np.mean(F[(us>=a)&(us<a+1)]**2))) for a in (3.0, 6.0, 9.0, 12.0)},
        hann_over_W={str(U): hann_mean(us, F**2, us[0], U)/h.W_scalar(mods['EPS'], g, L) for U in (9.0, 11.0, 13.0)})
    # ZETA L=4 row fix
    Fzeta = h.field(mods['ZETA'], g, L, us)
    out[f'ZETA_L{L}_hann_over_W'] = {str(U): hann_mean(us, Fzeta**2, us[0], U)/h.W_scalar(mods['ZETA'], g, L) for U in (9.0, 11.0, 13.0)}
    # BEUR vs quartet
    u0b = math.log(20)+L/2+0.05; ub = np.arange(u0b, 13.0, 0.004); Fb = h.field(mods['BEUR'], g, L, ub)
    Fq = np.zeros(len(ub))
    for r in mods['BEUR']['quartet']: Fq -= (G(g, L, complex(r)-0.5)*np.exp((complex(r)-0.5)*ub)).real
    out[f'BEUR_L{L}'] = dict(rms_field_by_u={str(a): float(np.sqrt(np.mean(Fb[(ub>=a)&(ub<a+1)]**2))) for a in (6.0, 9.0, 12.0)},
        rms_quartet_pred_by_u={str(a): float(np.sqrt(np.mean(Fq[(ub>=a)&(ub<a+1)]**2))) for a in (6.0, 9.0, 12.0)},
        rms_residual_by_u={str(a): float(np.sqrt(np.mean((Fb-Fq)[(ub>=a)&(ub<a+1)]**2))) for a in (6.0, 9.0, 12.0)})
    print(L, json.dumps({k: v for k, v in out.items() if k.endswith(f'L{L}') or f'L{L}_' in k}), flush=True)

# M-convergence and eigenvector frequency at L=4
for M in (40, 56, 72):
    arch, pole, prime, absprime = h.operators(mods['EPS'], 4.0, M); W = arch+pole-prime
    ev, V = np.linalg.eigh(W); v = V[:, 0]
    xs = np.linspace(-2, 2, 2001); f = basis(4.0, M, xs).T @ v
    rs = np.linspace(0, 80, 801); fh = np.abs([np.trapezoid(f*np.exp(1j*r*xs), xs) for r in rs])
    a2, p2, pr2, _ = h.operators(mods['ZETA'], 4.0, M)
    out[f'EPS_L4_M{M}'] = dict(lam_min=float(ev[0]), lam_2=float(ev[1]), peak_freq=float(rs[np.argmax(fh)]),
                              zeta_lam_min=float(np.linalg.eigvalsh(a2+p2-pr2)[0]),
                              rayleigh_parts=dict(arch=float(v@arch@v), pole=float(v@pole@v), prime=float(v@prime@v)))
    print(M, json.dumps(out[f'EPS_L4_M{M}']), flush=True)
# width scan for EPS first failure
for L in (3.25, 3.5, 3.75):
    arch, pole, prime, _ = h.operators(mods['EPS'], L, int(10*L)+4)
    out[f'EPS_lam_min_L{L}'] = float(np.linalg.eigvalsh(arch+pole-prime)[0]); print(L, out[f'EPS_lam_min_L{L}'], flush=True)
json.dump(out, open('harness_checks.json', 'w'), indent=1)
