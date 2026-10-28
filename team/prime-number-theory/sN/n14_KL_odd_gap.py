# N.14 NUMERICAL cross-check (of S&M's pending Arb certification) of the odd L^2 gap delta_o^0 = 1 - kappa^odd_L
# of CC's K_L at L = 1/2 log 2, plus the plain band projector (Golden Eye's 0.909/0.424/0.053).
# K_L: <eta|K xi> = (1/(2 eps'(1+))) int int eta(x) xi(x+v) (Q eps)(e^{|v|}) dx dv  (CC 2006.13771 (104)),
# Q eps(rho) = sum_n lambda(n)^2/(1-lambda(n)^2) C_n(rho)  (CC (99)), eps'(1+) = sum t(n) (CC Lemma 5.4),
# xi_n = even prolates for F_{e_R} on [-1,1] (c = 2*pi), first 11 terms (CC Lemma F.1: remainder <= 2.4e-12).
# Norm of each number: all L^2(-L, L). No Arb remainder: NUMERICAL.
import numpy as np
c = 2*np.pi
K = 120                                           # Legendre degree cutoff (even k only)
ks = np.arange(0, K, 2)
def x2kk(k): return (2*k*k + 2*k - 1)/((2*k + 3)*(2*k - 1))
def x2k2(k): return (k + 1)*(k + 2)/((2*k + 3)*np.sqrt((2*k + 1)*(2*k + 5)))
M = np.diag(ks*(ks + 1) + c*c*x2kk(ks)) + np.diag(c*c*x2k2(ks[:-1]), 1) + np.diag(c*c*x2k2(ks[:-1]), -1)
chi, D = np.linalg.eigh(M)
xg, wg = np.polynomial.legendre.leggauss(200)
Pn = np.array([np.polynomial.legendre.Legendre.basis(k)(xg)*np.sqrt((2*k + 1)/2) for k in ks])
NT = 11
psi = [D[:, n] @ Pn for n in range(NT)]           # values on Gauss nodes, int psi^2 = 1
def xian(n, x, deriv=False):                      # analytic continuation via xi = (1/nu) F xi
    x = np.atleast_1d(x)[:, None]
    ker = (-2*np.pi*xg[None, :]*np.sin(2*np.pi*x*xg[None, :])) if deriv else np.cos(2*np.pi*x*xg[None, :])
    return (ker*(wg*psi[n])[None, :]).sum(1)/nu[n]
nu = []
for n in range(NT):
    p0 = D[:, n] @ np.array([np.polynomial.legendre.Legendre.basis(k)(0.0)*np.sqrt((2*k + 1)/2) for k in ks])
    nu.append((wg*psi[n]).sum()/p0)
nu = np.array(nu); lam2 = nu**2
t = lam2/(1 - lam2)*np.array([xian(n, 1.0)[0]**2 for n in range(NT)])
print("lambda(n)^2:", lam2[:5])
print("t(n) (CC: 11.9719, 8.77574, 2.20528, 0.0433983, 0.000125459; normalization cancels in kappa):", t[:5])
print("ratios t(n)/t(0):", t[:5]/t[0], " CC ratios:", np.array([11.9719, 8.77574, 2.20528, 0.0433983, 0.000125459])/11.9719)
epsp = t.sum()
xq, wq = np.polynomial.legendre.leggauss(80)
def Qeps(rho):
    a0 = 1/rho; xs = a0 + (1 - a0)*(xq + 1)/2; ws = wq*(1 - a0)/2
    tot = 0.0
    for n in range(NT):
        d1 = xian(n, xs, True); d2 = xian(n, rho*xs, True)
        Cn = np.sqrt(rho)*np.sum(ws*xs*d1*rho*xs*d2) + rho**-1.5*xian(n, 1/rho, True)[0]*xian(n, 1.0)[0] \
             - rho**1.5*xian(n, 1.0)[0]*xian(n, rho, True)[0]
        tot += lam2[n]/(1 - lam2[n])*Cn
    return tot
for v in [1e-6, 0.1, 0.3, 0.5, np.log(2)]:
    print(f"kappa({v:.4f}) = Qeps(e^v)/(2eps') = {Qeps(np.exp(v))/(2*epsp): .6f}")
# Nystrom on (-L, L), L = 1/2 log 2
L = 0.5*np.log(2); Nn = 400
xn, wn = np.polynomial.legendre.leggauss(Nn); xn = L*xn; wn = L*wn
vgrid = np.linspace(0, 2*L, 1201); kv = np.array([Qeps(np.exp(v)) for v in vgrid])/(2*epsp)
kap = lambda v: np.interp(np.abs(v), vgrid, kv)
def spec(kernel):
    A = kernel(xn[:, None] - xn[None, :]); S = np.sqrt(wn)
    H = S[:, None]*A*S[None, :]
    R = np.eye(Nn)[::-1]                                     # x -> -x (nodes symmetric)
    Pe, Po = (np.eye(Nn) + R)/2, (np.eye(Nn) - R)/2
    ev = np.sort(np.linalg.eigvalsh(Pe @ H @ Pe))[::-1][:3]; od = np.sort(np.linalg.eigvalsh(Po @ H @ Po))[::-1][:3]
    return ev, od
T = 6.289836
ev, od = spec(lambda d: np.where(np.abs(d) < 1e-14, T/np.pi, np.sin(T*d)/(np.pi*np.where(d == 0, 1, d))))
print("band projector |t|<=6.289836 on (-L,L): even", ev, " odd", od)
ev, od = spec(kap)
print("K_L (CC kernel): top even (CC: 1.05177)", ev, " top odd", od)
print("delta_o^0 = 1 - kappa^odd =", 1 - od[0])
# C.14.A margin on |t| < 6.29 with sigma_S >= 0 only (indicative): 2 eps'(1+) delta/(1/4+t^2) - c2, min at t = 6.29
c2 = np.log(2)/np.sqrt(2); E2 = 2*22.9965
print("C.14.A indicative: 2eps'(1+) delta/(1/4+T^2) at T=6.29 (with sigma_S>=0 only) =", E2*(1 - od[0])/(0.25 + T*T), " vs c2 =", c2)
