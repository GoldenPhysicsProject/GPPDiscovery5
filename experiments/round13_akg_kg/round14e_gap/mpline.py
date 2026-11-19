"""Round 14E2: F(1/2+it) at 90 digits on t = 0, 0.1, ..., 260 for ZETA and F_theta = zeta(s-0.2)zeta(s+0.2) (mpmath zeta:
Dirichlet-series data via Euler-Maclaurin; no zeros).  -> mp80/line_<NAME>.pkl"""
import sys, os, pickle, mpmath as mp
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); DPS = 90; H = mp.mpf('0.1'); NT = 2601
def F(args):
    name, k = args; mp.mp.dps = DPS; s = mp.mpc(mp.mpf(1)/2, k*mp.mpf(1)/10)
    v = mp.zeta(s) if name == 'ZETA' else mp.zeta(s - mp.mpf('0.2'))*mp.zeta(s + mp.mpf('0.2'))
    return (mp.nstr(v.real, DPS+3), mp.nstr(v.imag, DPS+3))
if __name__ == '__main__':
    name = sys.argv[1]
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as p: v = p.map(F, [(name, k) for k in range(NT)], chunksize=20)
    pickle.dump(dict(name=name, dps=DPS, h='0.1', NT=NT, F=v), open(os.path.join(HERE, 'mp80', f'line_{name}.pkl'), 'wb'))
    print(name, 'done', v[0], v[-1], flush=True)
