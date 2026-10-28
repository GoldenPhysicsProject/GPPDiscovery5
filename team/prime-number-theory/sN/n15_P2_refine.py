# N.15 NUMERICAL: refinement of the Ritz upper bounds for the Weil form at 2L = log 3 (P2 = full form there)
# and at 2L = log 2 + 0.3, to see whether the small margins converge or keep dropping. Imports n15_P2_margins.
import sys, numpy as np; sys.argv = ['x']
src = open(__file__.replace('n15_P2_refine.py', 'n15_P2_margins.py')).read().split('c2, c3, c4 =')[0]
exec(src)
c2 = np.log(2)/np.sqrt(2)
for lab, twoL in [("log2+0.3", np.log(2)+0.3), ("log3", np.log(3))]:
    for N in [40, 80, 160, 320]:
        r = margins(twoL/2, N, [(2, c2)])
        print(f"2L={lab:9s} N={N:4d}  even <= {r[0]: .6f}   odd <= {r[1]: .6f}", flush=True)
