"""S84: on-line DH zeros in (T0,T1) via sign changes of the real Hardy function (see dh_zeros.py). NUMERICAL."""
import sys, json, mpmath as mp
from dh_zeros import Z
T0, T1 = float(sys.argv[1]), float(sys.argv[2]); h = 0.05
mp.mp.dps = 25
ts = [T0 + i*h for i in range(int((T1-T0)/h)+1)]; vals = [mp.re(Z(mp.mpf(t))) for t in ts]
mp.mp.dps = 40; on = []
for i in range(len(ts)-1):
    if vals[i]*vals[i+1] < 0: on.append(mp.findroot(lambda t: mp.re(Z(t)), (mp.mpf(ts[i]), mp.mpf(ts[i+1])), solver='anderson'))
json.dump([mp.nstr(x, 35) for x in on], open("dh_zeros_%d_%d.json" % (T0, T1), "w"))
print(T0, T1, len(on))
