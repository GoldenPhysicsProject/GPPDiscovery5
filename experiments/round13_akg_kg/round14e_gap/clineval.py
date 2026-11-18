"""Round 14E: complex F(1/2+it), t in [0,250] step 0.01, for ZETA, EPS (x^2+5y^2), DH, F_theta (theta=0.2), from Dirichlet
coefficients via Hurwitz sums (integer data). No zeros. Writes cline_<name>.npy (t, F)."""
import sys, numpy as np, mpmath as mp
from multiprocessing import Pool
mp.mp.dps = 15
T, DT = 250.0, 0.01
def L(s, q, chi): return q**(-s)*mp.fsum(chi[r]*mp.zeta(s, mp.mpf(r)/q) for r in range(1, q) if chi[r] != 0)
c4 = {r: [0, 1, 0, -1][r % 4] for r in range(4)}
c5 = {r: [0, 1, -1, -1, 1][r % 5] for r in range(5)}
c20 = {r: c4[r % 4]*c5[r % 5] for r in range(20)}
KAP = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
chi = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
dhc = {r: ((1-1j*KAP)/2*chi[r] + (1+1j*KAP)/2*mp.conj(chi[r])) for r in range(5)}
def F(name, s):
    if name == 'ZETA': return mp.zeta(s)
    if name == 'EPS': return mp.zeta(s)*L(s, 20, c20) + L(s, 4, c4)*L(s, 5, c5)
    if name == 'DH': return L(s, 5, dhc)
    if name == 'FTH': return mp.zeta(s-0.2)*mp.zeta(s+0.2)
def job(args):
    name, ts = args
    return [complex(F(name, mp.mpc(0.5, t))) for t in ts]
if __name__ == '__main__':
    name = sys.argv[1]
    t = np.arange(0, T+DT/2, DT); chunks = np.array_split(t, 32)
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 2) as p:
        v = np.concatenate([np.array(x) for x in p.map(job, [(name, c) for c in chunks])])
    np.save(f'cline_{name}.npy', np.vstack([t.astype(complex), v]))
    print(name, 'done', v[:3], flush=True)
