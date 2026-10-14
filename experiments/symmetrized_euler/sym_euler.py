"""Zeros of symmetrized truncated Euler products, zeta vs fakes. NUMERICAL (float, numpy/scipy + mpmath polish).
Lambda_N(s) = Gamma_R(s) prod_{p<=N}(1-p^-s)^-1, Gamma_R(s)=pi^{-s/2}Gamma(s/2).
 (Z)   f = Lambda_N(s) + Lambda_N(1-s)                       [analytic symmetrization, single Euler product]
 (DH)  f = Lam_D(s) + Lam_D(1-s), Lam_D = (5/pi)^{s/2}Gamma((s+1)/2)[(1-ik)/2 P_chi + (1+ik)/2 P_chibar]   [sum of 2 Euler products]
 (DHc) same with a single character: Lam_chi(s) + eps(chi) Lam_chibar(1-s)     [single Euler product, conductor 5]
 (FT)  f = A(s)+A(1-s), A = Lambda_N(s-th)Lambda_N(s+th)     [shifted, non-tempered]
 (G)   Gonek harmonic version P + chi * conj(P): zeros off the line impossible where |chi|!=1, FOR ANY P (lemma).
Counts zeros with Re s > 1/2 + 1e-6 and 0 < t <= T by grid minima of |f| + Newton polish.
Run: python sym_euler.py N T
"""
import sys, json, os
import numpy as np
from scipy.special import loggamma
from mpmath import mp, mpc, findroot, gamma as mg, pi as mpi, sqrt as msq
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from common import write_result
N = int(sys.argv[1]) if len(sys.argv) > 1 else 5; T = float(sys.argv[2]) if len(sys.argv) > 2 else 60
primes = [p for p in range(2, N + 1) if all(p % q for q in range(2, int(p**.5) + 1))]
k = (np.sqrt(10 - 2*np.sqrt(5)) - 2) / (np.sqrt(5) - 1)
chi = {1: 1, 2: 1j, 4: -1, 3: -1j, 0: 0}
def P(s, ch=None):
    out = np.ones_like(s)
    for p in primes:
        c = 1 if ch is None else ch[p % 5]
        out = out / (1 - c * p ** (-s))
    return out
GR = lambda s: np.exp(-s/2*np.log(np.pi) + loggamma(s/2))
G5 = lambda s: np.exp(s/2*np.log(5/np.pi) + loggamma((s+1)/2))
chib = {a: np.conj(b) for a, b in chi.items()}
tau = sum(chi[n]*np.exp(2j*np.pi*n/5) for n in (1, 2, 3, 4)); eps = tau/(1j*np.sqrt(5))
th = 0.25
fams = {
 "Z":  lambda s: GR(s)*P(s) + GR(1-s)*P(1-s),
 "FT": lambda s: GR(s-th)*P(s-th)*GR(s+th)*P(s+th) + GR(1-s-th)*P(1-s-th)*GR(1-s+th)*P(1-s+th),
 "DH": lambda s: (lambda D: D(s) + D(1-s))(lambda u: G5(u)*((1-1j*k)/2*P(u, chi) + (1+1j*k)/2*P(u, chib))),
 "DHc": lambda s: G5(s)*P(s, chi) + eps*G5(1-s)*P(1-s, chib),
}
sig = np.arange(0.5005, 3.0, 0.01); tt = np.arange(0.3, T, 0.02)
S = sig[None, :] + 1j*tt[:, None]
res = {}
for name, f in fams.items():
    V = np.abs(f(S))
    m = (V[1:-1, 1:-1] < V[:-2, 1:-1]) & (V[1:-1, 1:-1] < V[2:, 1:-1]) & (V[1:-1, 1:-1] < V[1:-1, :-2]) & (V[1:-1, 1:-1] < V[1:-1, 2:])
    idx = np.argwhere(m); found = []
    for i, j in idx:
        s0 = complex(S[i+1, j+1])
        try:
            g = lambda z: complex(f(np.array([complex(z)]))[0])
            # secant polish in complex plane
            z0, z1 = s0, s0 + 1e-3
            for _ in range(60):
                g0, g1 = g(z0), g(z1)
                if g1 == g0: break
                z0, z1 = z1, z1 - g1*(z1-z0)/(g1-g0)
                if abs(z1-z0) < 1e-13: break
            if abs(g(z1)) < 1e-9*max(1, abs(g(z1+1e-4))/1e-4*1e-4) and z1.real > 0.5 + 1e-6 and 0 < z1.imag <= T:
                if all(abs(z1-w) > 1e-6 for w in found): found.append(z1)
        except Exception: pass
    found.sort(key=lambda z: z.imag); res[name] = [[round(z.real, 10), round(z.imag, 10)] for z in found]
    print(f"N={N} {name}: {len(found)} zeros with Re s>1/2 in 0<t<={T}:", res[name][:8])
write_result(os.path.join(os.path.dirname(__file__), f"result_N{N}_T{int(T)}.json"),
  {"experiment": "sym_euler", "label": "NUMERICAL", "N": N, "T": T, "primes": primes, "offline_zeros_right": res},
  kind="float", method="numpy grid minima + secant polish")
