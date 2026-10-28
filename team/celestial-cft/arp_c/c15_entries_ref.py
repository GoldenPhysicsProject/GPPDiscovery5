"""C.15(c) reference table (mpmath, 30 digits, NUMERICAL): B_jk closed form, A_jk = Q'_T(f_j,f_k) at ell=0.1, T=60.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np, mpmath as mp
from c15_even_schur import Bmat
mp.mp.dps = 30
a = mp.log(2); c2 = a/mp.sqrt(2)
def Om(t): return mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*t)) - mp.log(mp.pi)
def om(k, L): return (k + mp.mpf(1)/2)*mp.pi/L
def ehat(k, t, L):
    w = om(k, L)
    if abs(w - t) < mp.mpf(10)**-20: return mp.sqrt(L)
    return (-1)**k*2*w*mp.cos(t*L)/(mp.sqrt(L)*(w*w - t*t))
def Bjk(j, k, L):
    wj, wk = om(j, L), om(k, L)
    def C(nu, ph):
        if nu == 0: return (2*L - a)*mp.cos(ph)
        return (mp.sin(nu*(L - a) + ph) - mp.sin(-nu*L + ph))/nu
    return ((wj*wk + mp.mpf(1)/4)*C(wj - wk, -wk*a) + (mp.mpf(1)/4 - wj*wk)*C(wj + wk, wk*a))/(2*L)
def Ajk(j, k, L, T):
    qT = Om(T) - 2*c2
    f = lambda t: (mp.mpf(1)/4 + t*t)*(Om(t) - 2*c2*mp.cos(a*t) - qT)*ehat(j, t, L)*ehat(k, t, L)
    pts = sorted(set([0] + [float(x) for x in mp.linspace(0, T, 31)] + [float(om(i, L)) for i in (j, k) if om(i, L) < T] + [T]))
    I = mp.quad(f, pts)/mp.pi
    nj, nk = mp.sqrt(om(j, L)**2 + mp.mpf(1)/4), mp.sqrt(om(k, L)**2 + mp.mpf(1)/4)
    return (qT if j == k else 0) + I/(nj*nk)
ell = mp.mpf('0.1'); L = (a + ell)/2; T = 60
print("ell = 0.1, L =", mp.nstr(L, 20), ", T = 60, q_T = Omega(60) - 2c2 =", mp.nstr(Om(T) - 2*c2, 15))
print("omega_k (k=0..3):", [mp.nstr(om(k, L), 15) for k in range(4)])
print("B_jk (unnormalised, closed form), j,k = 0..3:")
for j in range(4): print("  ", [mp.nstr(Bjk(j, k, L), 15) for k in range(4)])
Bf = Bmat(4, float(L)); print("  float quadrature check, max diff:", float(max(abs(Bf[j][k] - float(Bjk(j, k, L))) for j in range(4) for k in range(4))))
print("A_jk = (Q'_T)(f_j, f_k), H-normalised f_k = e_k/sqrt(omega_k^2 + 1/4), j,k = 0..3:")
for j in range(4): print("  ", [mp.nstr(Ajk(j, k, L, T), 12) for k in range(4)])
