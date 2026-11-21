"""Round 32 shared pieces: the pole/archimedean ('non-prime') part P_M(L) of the R31 floor D_M(L) = P_M(L) - S_M(L),
for ZETA (GPPVerify closed form, main c81c031), F_theta (any theta), EPS, DH (explicit formula, r22_fejer at g0 = 0),
and the Beurling main term (no FE, no archimedean place). No zeros."""
import math, numpy as np, mpmath
import r22_models as RM
from r22_fejer import K as _K, arch_I

c1 = float((mpmath.digamma(0.25) - mpmath.log(mpmath.pi))/2)
C2 = float(mpmath.nsum(lambda k: 1/(2*k + 0.5)**2, [0, mpmath.inf]))
_b = 2*np.arange(0, 400) + 0.5
def P_zeta(L):
    return 4*(math.exp(L/2) - 2 + math.exp(-L/2)) + c1*L + C2 - float(np.sum(np.exp(-_b*L)/_b**2))

def model(name, th=None):
    if name == 'FTH':
        return dict(logQ=-math.log(math.pi), arch=[(0.5, 0.25 - th/2, 0.5), (0.5, 0.25 + th/2, 0.5)], poles=[1 + th, 1 - th])
    return RM.MODELS[name]

def P_model(name, L, th=None):
    M = model(name, th); g0 = np.array([0.0])
    out = 2*(1/np.pi)*M['logQ']
    for (lmb, a, b) in M['arch']:
        out += 2*lmb*float(arch_I(a, b, L, g0)[0])
    for z in M['poles']:
        r = -1j*(z - 0.5)
        out += 2*float((_K(r - g0, L) + _K(r + g0, L)).real[0])
    return math.pi*L/2*out

def P_beur(L, x0=20.0):
    u0 = math.log(x0)
    if L <= u0: return 0.0
    return -2*(L - u0)*math.exp(u0/2) + 4*(math.exp(L/2) - math.exp(u0/2))
