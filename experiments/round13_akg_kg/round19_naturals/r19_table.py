"""Round 19, Task 2 coefficient checks on n <= 3000 (CHECKED by computation; the FE/degree columns are cited theory).
ZETA a_n = 1; EPS a_n = r_Q(n)/2, Q = x^2+5y^2 (Epstein, zeta_Q = sum' Q^{-s}/2); DH periodic mod 5;
F_theta = zeta(s-theta)zeta(s+theta), theta = 1/4: a_n = sum_{d|n} d^theta (n/d)^{-theta}.
Tests: all-ones on N, multiplicative (coprime), completely multiplicative, Ramanujan |a_p| <= degree."""
import math, json, numpy as np
N = 3000
def eps():
    r = np.zeros(N+1)
    for x in range(-60, 61):
        for y in range(-30, 31):
            v = x*x+5*y*y
            if 0 < v <= N: r[v] += 1
    return r/2
xi = (math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)
dhc = {0: 0.0, 1: 1.0, 2: xi, 3: -xi, 4: -1.0}
th = 0.25
A = dict(ZETA=(np.ones(N+1), 1), EPS=(eps(), 2), DH=(np.array([dhc[n % 5] for n in range(N+1)]), 1),
         FTHETA=(np.array([0.0]+[sum(d**th*(n//d)**(-th) for d in range(1, n+1) if n % d == 0) for n in range(1, N+1)]), 2))
P = [p for p in range(2, N+1) if all(p % q for q in range(2, int(p**0.5)+1))]
out = {}
for k, (a, deg) in A.items():
    mult = max(abs(a[m*n]-a[m]*a[n]) for m in range(2, 55) for n in range(2, 55) if math.gcd(m, n) == 1 and m*n <= N)
    cmult = max(abs(a[m*n]-a[m]*a[n]) for m in range(2, 55) for n in range(2, 55) if m*n <= N)
    ram = max(abs(a[p]) for p in P) / deg
    ex = None
    if mult > 1e-9:
        ex = next([m, n, float(a[m*n]), float(a[m]*a[n])] for m in range(2, 55) for n in range(2, 55) if math.gcd(m, n) == 1 and m*n <= N and abs(a[m*n]-a[m]*a[n]) > 1e-9)
    row = dict(all_ones_on_N=bool(np.all(np.abs(a[1:]-1) < 1e-12)), max_multiplicativity_defect=float(mult),
               max_complete_multiplicativity_defect=float(cmult), max_ap_over_degree=float(ram), example_defect=ex)
    out[k] = row; print(k, row)
json.dump(out, open('r19_table.json', 'w'), indent=1)
