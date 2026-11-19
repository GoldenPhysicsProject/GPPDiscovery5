"""Round 14E (iii) driver: optimal in-window probe at the Lehmer pair (#6709, #6710, t ~ 7005.06) vs a median-gap pair and the
first two zeros. Zeros = CHECK inputs (they choose the probe and give the zero side); W(g*) is recomputed from integer data."""
import json, numpy as np, mpmath as mp
from multiprocessing import Pool
from lehmer import *
def zz(k):
    mp.mp.dps = 20; return float(mp.zetazero(k).imag)
if __name__ == '__main__':
    with Pool(8) as p:
        zs = p.map(zz, range(6560, 6861)); lo = p.map(zz, range(1, 120))
    k0 = 6560; d = np.diff(zs); i = int(np.argmin(d)); mid = np.arange(100, 200); j = int(mid[np.argmin(np.abs(d[mid]-np.median(d)))])
    out = []
    for L in [3.0, 5.0, 7.0, 8.0]:
        for lab, Z_, ii in [(f'Lehmer #{k0+i},{k0+i+1}', zs, i), (f'median-gap #{k0+j},{k0+j+1}', zs, j), ('first two zeros #1,2', lo, 0)]:
            r = optimal(L, Z_, ii, half=15.0); r['pair'] = lab; out.append(r)
            print(lab, 'L=%g gap=%.4f W_data=%.4g W_zero=%.4g q=%.3f Rq<=%.4g ideal=%.3g' % (L, r['gap'], r['W_data'], r['W_zero_side'], r['q'], r['Rq_upper_bound'], r['ideal_bound']), flush=True)
    json.dump(dict(zeros_k0=k0, rows=out), open('lehmer_opt.json', 'w'), indent=1)
