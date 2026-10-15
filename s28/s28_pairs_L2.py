"""Grok S28: L = 2 grouped torus lift (NUMERICAL). sup over T^n of a sum <= sum of sups over groups, so
A <= sum_groups A_group with each group's lift built from that group's (cross-)invisible characters.
Groups: consecutive pairs of primes (2,3),(5,7),(11,13),...; compare with per-prime (axis-only, same kmax).  RH OPEN."""
import json, math, sys
sys.path.insert(0, '.')
from s28_torus_lp import solve, primes_upto
L = 2.0; D = int(sys.argv[1]) if len(sys.argv) > 1 else 4
ps = [p for p in primes_upto(math.exp(2 * L) + 1) if math.log(p) < 2 * L]
groups = [ps[i:i + 2] for i in range(0, len(ps), 2)]
tot_pair = 0.0; tot_axis = 0.0; out = []
for g in groups:
    rp = solve(L, D=D, plist=g, npts=8000, rounds=12, verbose=False)
    ra = solve(L, D=0, plist=g, npts=8000, rounds=12, verbose=False)
    tot_pair += rp['A_max_est']; tot_axis += ra['A_max_est']
    out.append(dict(group=g, A_pair=rp['A_max_est'], A_axis=ra['A_max_est'], raw=rp['per_prime_raw']))
    print(g, f"raw {rp['per_prime_raw']:.4f} axis-only {ra['A_max_est']:.4f} pair-torus {rp['A_max_est']:.4f}", flush=True)
print(f'L=2 D={D}: per-prime(axis, finite kmax) {tot_axis:.4f}  pair-torus {tot_pair:.4f}  (S27 per-prime 13.61; 2e^L = {2*math.exp(2):.3f})')
json.dump(dict(L=L, D=D, groups=out, A_axis=tot_axis, A_pair=tot_pair), open(f'pairs_L2_D{D}.json', 'w'), indent=1)
