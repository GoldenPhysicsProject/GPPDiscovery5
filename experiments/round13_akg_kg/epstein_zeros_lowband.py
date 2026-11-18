"""Supplement to epstein_zeros.py (test data only): the first pass searched off-line zeros only for sigma >= 0.55.
Gap 13.10 -> 17.60 in the on-line list suggests a near-line off-line pair. Search sigma in (0.5, 0.55], t in [0.5, 60]
via findroot from a grid, and confirm counts with the argument principle on the rectangle [0.5+1e-3, 1.6] x [0.5, 60]."""
import mpmath as mp, json
mp.mp.dps = 20
src = open('epstein_zeros.py').read().split('# coefficient check')[0]
ns = {}; exec(src, ns); Z = ns['Z']
found = []
for t0 in [x*0.5 for x in range(1, 120)]:
    for s0 in (0.51, 0.53, 0.56):
        try:
            z = mp.findroot(Z, mp.mpc(s0, t0))
        except Exception:
            continue
        if 0.5005 < z.real < 1.6 and 0.3 < z.imag < 60.5 and abs(Z(z)) < 1e-12:
            if all(abs(z - f) > 1e-6 for f in found): found.append(z)
found.sort(key=lambda z: z.imag)
def argcount(a, b, c, d, n=4000):
    path = [mp.mpc(a, c), mp.mpc(b, c), mp.mpc(b, d), mp.mpc(a, d), mp.mpc(a, c)]
    tot = 0
    for p, q in zip(path[:-1], path[1:]):
        prev = Z(p)
        for k in range(1, n+1):
            cur = Z(p + (q-p)*k/n); tot += mp.im(mp.log(cur/prev)); prev = cur
    return float(tot/(2*mp.pi))
cnt = argcount(0.501, 1.6, 0.5, 60.5)
out = dict(offline_right_half=[[float(z.real), float(z.imag)] for z in found], argument_principle_count=cnt)
print(json.dumps(out)); json.dump(out, open('epstein_zeros_lowband.json', 'w'), indent=1)
