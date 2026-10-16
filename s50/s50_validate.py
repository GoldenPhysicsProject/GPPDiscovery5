# S50 (Grok): (i) direct mpmath validation of the zero-form L_n at one planted point; (ii) actual DH L_n, n=1..4;
# (iii) F_theta (theta=1/4, all zeros off line) L_n near low and high zeros.
import json, math, numpy as np, mpmath as mp
G = np.loadtxt('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt', usecols=1); G2 = G**2
res = {}
def Ln_from_taylor(fk, n):  # fk = Taylor coefficients f^(k)/k!, L_n = sum (-1)^{n+k} f_k f_{2n-k}
    return sum((-1)**(n + k)*fk[k]*fk[2*n - k] for k in range(2*n + 1))
mp.mp.dps = 40
def xi(s): return s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
F = lambda v: mp.re(xi(mp.mpf(1)/2 + 1j*mp.sqrt(-v)))
gp = 0.5*(G[22] + G[23]); d = 0.3; al = d*d - gp*gp; be = 2*d*gp; x0 = al + 0.3*be
a = mp.mpc(d, gp)**2
tay = mp.taylor(F, mp.mpf(x0), 8, chop=False)
q = [1, -2*mp.re(1/a), 1/abs(a)**2]  # (1-u/a)(1-u/abar) expanded at u=0 -> re-expand at x0
# Taylor of quadratic at x0: Q(x0+h) = Q(x0) + Q'(x0) h + q2 h^2
Q0 = q[0] + q[1]*x0 + q[2]*x0**2; Q1 = q[1] + 2*q[2]*x0; Q2 = q[2]
fk = [sum(tay[i]*[Q0, Q1, Q2][k - i] for i in range(max(0, k - 2), k + 1)) for k in range(9)]
direct = [Ln_from_taylor(fk, n)/fk[0]**2 for n in range(1, 5)]
# zero form at x0 (unscaled)
w = 1/(x0 + G2)**2; E = np.zeros(6); E[0] = 1
for wj in w: E[1:] = E[1:] + wj*E[:-1]
Dv = (al - x0)**2; p0 = (Dv + be*be)**2; p1 = 2*(Dv - be*be)
zf = [(p0*E[n] + p1*E[n - 1] + (E[n - 2] if n >= 2 else 0))/p0 for n in range(1, 5)]
res['validation'] = dict(gamma=float(gp), delta=d, x0=x0, direct=[mp.nstr(v, 12) for v in direct], zero_form=[f'{v:.11e}' for v in zf],
                         rel_diff=[float(abs(direct[i] - zf[i])/abs(zf[i])) for i in range(4)])
print(res['validation'], flush=True)
# (ii) DH actual (S49 definition), off-line zero 0.808517+85.699348i -> delta=0.308517
mp.mp.dps = 30
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
chl = [0, 1, 1j, -1j, -1]; chb = [0, 1, -1j, 1j, -1]
def L(s, ch): return 5**(-s)*sum(ch[r]*mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5))
def Lam(s): return (5/mp.pi)**((s + 1)/2)*mp.gamma((s + 1)/2)*((1 - 1j*k)/2*L(s, chl) + (1 + 1j*k)/2*L(s, chb))
fu = lambda v: mp.re(Lam(mp.mpf(0.5) + 1j*mp.sqrt(-v)))
d0, g0 = mp.mpf('0.308517'), mp.mpf('85.699348'); aD = float(d0**2 - g0**2); bD = float(2*d0*g0)
rows = []
for s in [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9]:
    tk = mp.taylor(fu, mp.mpf(aD + s*bD), 8, chop=False)
    rows.append(dict(s=s, L=[mp.nstr(Ln_from_taylor(tk, n)/tk[0]**2, 6) for n in range(1, 5)]))
    print('DH', rows[-1], flush=True)
res['DH_actual'] = dict(alpha=aD, beta=bD, rows=rows, note='x = alpha + s*beta; L_n/f^2')
# (iii) F_theta theta=1/4: every u-zero (theta+i g_j)^2; L_n/f^2 = [t^n] prod_j P_j(t)/P_j(0)
th = 0.25
def Ftheta_L(xs, Kmax=6):
    A = th*th - G2; B = 2*th*G
    E = np.zeros((len(xs), Kmax + 1)); E[:, 0] = 1
    for Aj, Bj in zip(A, B):
        Dv = (Aj - xs)**2; p0 = (Dv + Bj*Bj)**2
        c1 = 2*(Dv - Bj*Bj)/p0; c2 = 1/p0
        E[:, 2:] = E[:, 2:] + c1[:, None]*E[:, 1:-1] + c2[:, None]*E[:, :-2]
        E[:, 1] = E[:, 1] + c1*E[:, 0]
    return E
ft = {}
for j in [0, 99, 9999, 99000]:
    Aj, Bj = th*th - G2[j], 2*th*G[j]
    xs = Aj + Bj*np.linspace(-1, 1, 41)
    E = Ftheta_L(xs)
    ft[str(j + 1)] = dict(gamma=float(G[j]), eps=float(th/(2*math.pi/math.log(G[j]/(2*math.pi)))) if G[j] > 2*math.pi*1.5 else None,
                          min_L=[float(E[:, n].min()) for n in range(1, 7)], L_at_alpha=[float(E[20, n]) for n in range(1, 7)])
    print('Ftheta zero', j + 1, ft[str(j + 1)], flush=True)
res['F_theta_quarter'] = ft
json.dump(res, open('s50_validate.json', 'w'), indent=1, default=str)
