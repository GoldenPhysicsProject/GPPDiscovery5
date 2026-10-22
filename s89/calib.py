"""S89: calibrate zero-side Legendre vectors against the S87 geometric builder (float check).
u_i(z) = 2 sqrt(L) sqrt(2i+1) Re/Im[ e^{i zL/2} i^i j_i(zL/2) ]  (even: cos -> Re, odd: sin -> Im)."""
import sys, numpy as np, mpmath as mp
from scipy.special import spherical_jn
sys.path.insert(0, "../s87"); import hp_forms as H
g = np.loadtxt('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt', usecols=1)
def U(gg, L, N, par):
    a = gg*L/2; out = np.empty((len(gg), N))
    for i in range(N):
        z = np.exp(1j*a)*(1j**i)*spherical_jn(i, a)
        out[:, i] = 2*np.sqrt(L)*np.sqrt(2*i+1)*(z.real if par == 'e' else z.imag)
    return out
for L in (0.5, 1.0):
    N = 8; H.setup(40); d = H.build('ZETA', L, N)
    for par in 'eo':
        A = np.array(mp.matrix(d['A'+par]).tolist(), float); P = np.array(mp.matrix(d['P'+par]).tolist(), float)
        Q = A - P + sum(np.outer(np.array(p, float), np.array(p, float)) for p in d['p'+par])
        Uz = U(g, L, N, par); G = 2*Uz.T@Uz
        print(L, par, 'maxabs Q-G', np.abs(Q-G).max(), 'Q00', Q[0,0], 'G00', G[0,0], 'Q11', Q[1,1], G[1,1])
print('--- odd sign test, half-sum zeros')
for L in (0.5, 1.0):
    N = 8; H.setup(40); d = H.build('ZETA', L, N)
    for par in 'eo':
        A = np.array(mp.matrix(d['A'+par]).tolist(), float); P = np.array(mp.matrix(d['P'+par]).tolist(), float)
        pp = sum(np.outer(np.array(p, float), np.array(p, float)) for p in d['p'+par])
        Uz = U(g, L, N, par); G = Uz.T@Uz
        print(L, par, '|A-P+pp - G|', np.abs(A-P+pp-G).max(), '|A-P-pp - G|', np.abs(A-P-pp-G).max())
