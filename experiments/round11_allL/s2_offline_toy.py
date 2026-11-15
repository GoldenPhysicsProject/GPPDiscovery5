"""Round 11, script 2 (NUMERICAL toy, hypothetical data). How a single off-line tau-pair shows up in the tau-Gram,
as a function of the window L.  Measure = the first NZ on-line ordinates (comparison only) + ONE hypothetical
off-line quartet rho = 1/2 +- beta + i*(+-gamma0). Form: Q_L(c) = sum_rho PsiHat(z_rho) F_c(z_rho) F_c(-z_rho),
z_rho = (rho-1/2)/i, x_i on a grid in [0,L].  Reports lambda_min(Q_L)/||Q_L|| vs L, and lambda_min for the
on-line part alone (always >= 0: Gram/Bochner).  Also the 2x2 'PSS block' for the pair {rho, tau rho}."""
import numpy as np, mpmath as mp, json
NZ = 150
gam = np.array([float(mp.zetazero(k).imag) for k in range(1, NZ+1)])
PsiHat = lambda z: (np.sin(z/2)/(z/2))**4
def Qmat(xs, zs):
    D = xs[:, None]-xs[None, :]
    M = np.zeros_like(D, dtype=complex)
    for z in zs: M += PsiHat(z)*np.exp(1j*z*D)
    return M.real, np.abs(M.imag).max()
on = list(gam)+list(-gam)
res = []
for gamma0, beta in ((28.0, 0.1), (28.0, 0.02), (28.0, 0.25), (100.3, 0.1)):
    quart = [gamma0-1j*beta, gamma0+1j*beta, -gamma0-1j*beta, -gamma0+1j*beta]
    for L in (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0):
        h = min(0.04, 2*np.pi/gamma0/8); n = int(round(L/h))+1
        xs = np.linspace(0, L, n)
        Mon, _ = Qmat(xs, on); Moff, im = Qmat(xs, quart)
        ev = np.linalg.eigvalsh(Mon+Moff); ev0 = np.linalg.eigvalsh(Mon)
        nrm = np.abs(ev).max()
        res.append(dict(gamma0=gamma0, beta=beta, L=L, n=n, lmin_rel=ev[0]/nrm, lmin=ev[0], lmin_online_only=ev0[0],
                        offline_block_eigs=list(np.linalg.eigvalsh(Moff)[[0, -1]]), imag_residual=im))
        print({k: (float('%.4g' % v) if isinstance(v, float) else v) for k, v in res[-1].items() if k != 'offline_block_eigs'}, flush=True)
json.dump(res, open('s2_results.json', 'w'), indent=1, default=float)
