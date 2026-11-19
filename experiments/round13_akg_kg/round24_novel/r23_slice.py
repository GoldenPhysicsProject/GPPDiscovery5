"""Round 23: the gamma = 0 slice of GAP-22, D(L) := (pi L/2) rho_{e^L}(0).
Closed form for zeta (PROVED, inverse Laplace of xi'/xi(1/2+sigma)/sigma^2 term by term):
  D(L) = 4(e^{L/2} - 2 + e^{-L/2}) + c1 L + C2 - T(L) - S(L),
  c1 = (psi(1/4) - log pi)/2,  C2 = psi'(1/4)/4,  T(L) = sum_{k>=0} e^{-(2k+1/2)L}/(2k+1/2)^2,
  S(L) = sum_{n <= e^L} Lambda(n) n^{-1/2} (L - log n).
Laplace transform: int_0^inf e^{-sigma L} D(L) dL = (xi'/xi)(1/2+sigma)/sigma^2 (Re sigma > 1/2).
Under RH, D(L) = sum_rho (1 - cos(gamma L))/gamma^2 (sum over all zeros).  No zeros are inputs; zeros used only to
validate (illustration). Fakes: D_M(L) via the general explicit formula (r22_fejer.rho at gamma0 = 0)."""
import json, math, numpy as np, mpmath as mp
from r22_fejer import rho
from r22_models import lam_table
mp.mp.dps = 30
c1 = float((mp.digamma(0.25) - mp.log(mp.pi))/2); C2 = float(mp.psi(1, 0.25)/4)
def D_zeta(L, lam):
    X = math.exp(L); n = np.arange(2, int(X) + 1); lm = lam[2:int(X) + 1]; nz = lm != 0
    S = float(np.sum(lm[nz]*n[nz]**-0.5*(L - np.log(n[nz]))))
    k = np.arange(0, 60); T = float(np.sum(np.exp(-(2*k + 0.5)*L)/(2*k + 0.5)**2))
    return 4*(math.exp(L/2) - 2 + math.exp(-L/2)) + c1*L + C2 - T - S
if __name__ == '__main__':
    out = dict(c1=c1, C2=C2, two_c0_mean_under_RH=0.0462099862308379)
    lamZ = lam_table('ZETA', 10**7 + 1)
    # cross-check closed form vs r22 general explicit formula, and vs zeros (illustration)
    gam = np.array([float(mp.zetazero(k).imag) for k in range(1, 401)])
    chk = []
    for L in (0.8, 1.5, 2.3, 4.0, 6.0, 9.0):
        a = D_zeta(L, lamZ); b = float(np.pi*L/2*rho('ZETA', math.exp(L), np.array([0.0]), lamZ)[0])
        z = float(2*np.sum((1 - np.cos(gam*L))/gam**2))
        chk.append(dict(L=L, closed_form=a, r22_formula=b, zero_side_400=z)); print(chk[-1], flush=True)
    out['crosscheck'] = chk
    Ls = np.arange(math.log(2), math.log(1e7), 0.0005)
    D = np.array([D_zeta(L, lamZ) for L in Ls])
    i = int(np.argmin(D))
    out['zeta'] = dict(L_range=[float(Ls[0]), float(Ls[-1])], min=float(D[i]), argmin_L=float(Ls[i]), argmin_X=float(math.exp(Ls[i])),
                       mean=float(D.mean()), mean_L_ge_8=float(D[Ls >= 8].mean()),
                       min_by_decade={str(d): float(D[(Ls >= math.log(10**d)) & (Ls < math.log(10**(d + 1)))].min()) for d in range(0, 7)})
    print('zeta', out['zeta'], flush=True)
    np.save('r23_D_zeta.npy', np.vstack([Ls, D]))
    for name in ('FTH', 'EPS', 'DH'):
        lam = lam_table(name, 10**5 + 1)
        Lf = np.arange(math.log(2), math.log(1e5), 0.002)
        Df = np.array([np.pi*L/2*rho(name, math.exp(L), np.array([0.0]), lam)[0] for L in Lf])
        neg = Lf[Df < 0]
        out[name] = dict(first_negative_L=float(neg[0]) if len(neg) else None,
                         first_negative_X=float(math.exp(neg[0])) if len(neg) else None, min=float(Df.min()),
                         frac_negative=float(np.mean(Df < 0)))
        print(name, out[name], flush=True)
    json.dump(out, open('r23_slice.json', 'w'), indent=1)
