"""Round 13c: attack O9 -- is there a source of normalizability built from primes, pole, arch only?
Each candidate 'conservation law' is tested WITHOUT zeros; a hypothetical FE-symmetric off-line quartet
(beta=0.1, gamma0=28) is used only as a toy to see whether the law could exclude it.
(a) KG flux: J(U) = int_{u0}^U (Phi_g^2 - W) du; conservation <=> J bounded.
(b) Functional-equation time reversal: Phi_g(-u) == Phi_{g~}(u) from the prime side; does the toy also satisfy it?
(c) Unitary scattering S(s) = xi(2s-1)/xi(2s) (Lax-Phillips/Pavlov-Faddeev): |S(1/2+it)| = 1 from zeta; toy too?
(d) Arch-only initial energy: is W(g) = (arch+pole)(g) (primes dropped)?
"""
import numpy as np, math, json, sys
import mpmath as mp
sys.path.insert(0, '.')
from akg_kg import vonmangoldt_table, make_g, W_parts, Phi, gv, w
NS, LAM = vonmangoldt_table(int(math.exp(14.05))+1)
out = {}
du = 0.004
# (a) flux
fl = []
for L in (0.5, 1.0, 2.0):
    g = make_g(L, 'even'); W = W_parts(g, L, NS, LAM)['W']
    us = np.arange(2.0, 13.0+du/2, du); F, _ = Phi(g, L, us, NS, LAM)
    J = np.concatenate([[0], np.cumsum(0.5*(F[1:]**2+F[:-1]**2 - 2*W)*du)])
    # toy quartet added to the same field
    gt = gv(g, L, np.array([28-0.1j, 28+0.1j]))
    T = F + 2*np.real(gt[0]*np.exp(28j*us)*np.exp(0.1*us)) + 2*np.real(gt[1]*np.exp(28j*us)*np.exp(-0.1*us))
    JT = np.concatenate([[0], np.cumsum(0.5*(T[1:]**2+T[:-1]**2 - 2*W)*du)])
    blocks = [float(np.max(np.abs(J[(us >= a) & (us < a+2.75)]))/W) for a in (2.0, 4.75, 7.5, 10.25)]
    blocksT = [float(np.max(np.abs(JT[(us >= a) & (us < a+2.75)]))/W) for a in (2.0, 4.75, 7.5, 10.25)]
    fl.append(dict(L=L, W=W, maxabsJ_over_W_by_quarter=blocks, J_end_over_W=float(J[-1]/W),
                   toy_maxabsJ_over_W_by_quarter=blocksT, toy_J_end_over_W=float(JT[-1]/W)))
out['a_flux'] = fl
# (b) FE time reversal, from the prime side, for an asymmetric g
L = 1.0; g = make_g(L, 'asym'); gr = lambda x: g(-x); gr.poly = g.poly.copy(); gr.poly.coef = gr.poly.coef*np.array([(-1)**k for k in range(len(gr.poly.coef))]); gr.a = g.a
us = np.array([3.0, 6.5, 10.0])
gm = gv(g, L, np.array([0.5j, -0.5j]))
lg = np.log(NS); c = LAM/np.sqrt(NS)
def Phi_neg(u):  # Phi_g(-u) for u > L/2: explicit formula with G = g(. + u) supported in x<0
    i0, i1 = np.searchsorted(lg, [u-L/2, u+L/2])
    ps = np.sum(c[i0:i1]*g(-(lg[i0:i1]) + u))
    xs = np.linspace(u-L/2, u+L/2, 4001); ar = np.trapezoid(np.array([w(t) for t in xs])*g(-xs+u), xs)
    return (math.exp(u/2)*gm[0] + math.exp(-u/2)*gm[1]).real - ps - ar
F_rev, _ = Phi(gr, L, us, NS, LAM)
out['b_FE'] = dict(Phi_g_minus_u=[float(Phi_neg(u)) for u in us], Phi_grev_u=[float(x) for x in F_rev],
                   note='toy quartet {1/2+-0.1+-28i} is FE- and conjugation-symmetric by construction, so it satisfies the same reversal law')
# (c) Lax-Phillips scattering unitarity on Re s = 1/2, from zeta (no zeros); toy multiplies xi by P(s), P(s)=P(1-s), real coeffs
mp.mp.dps = 30
xi = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
rho = [mp.mpf('0.6')+28j, mp.mpf('0.6')-28j, mp.mpf('0.4')+28j, mp.mpf('0.4')-28j]
P = lambda s: mp.fprod([s - r for r in rho])
rows = []
for t in (0.3, 3.0, 7.07, 14.0, 20.0):
    s = mp.mpf('0.5') + 1j*t
    S = xi(2*s-1)/xi(2*s); St = S*P(2*s-1)/P(2*s)
    rows.append(dict(t=t, absS=float(abs(S)), absS_toy=float(abs(St))))
out['c_unitarity'] = rows
# (d) arch-only initial energy
r = json.load(open('akg_kg_results.json'))
out['d_arch_only'] = [dict(L=x['L'], g=x['g'], W=x['W'], arch_plus_pole=x['arch']+x['pole'], prime=x['prime'], rel_gap=x['prime']/x['W']) for x in r['runs']]
json.dump(out, open('o9_conservation_results.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
