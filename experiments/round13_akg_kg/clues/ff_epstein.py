"""Clue probe B: an exactly solvable, fully FINITE 'Epstein' in a function field.
A = F_q[E - O], E elliptic over F_q, Cl(A) = Pic^0(E) = E(F_q), h = #E(F_q), a = q + 1 - h.
Partial (ideal-class) zeta Z_c(T) = sum_{I in class c} T^{deg I}. Riemann-Roch (g = 1) gives
  Z_0(T) = (1 - qT + qT^2)/(1 - qT),   Z_c(T) = T/(1 - qT) (c != 0),
  sum_c Z_c = (1 - aT + qT^2)/(1 - qT) = zeta_A(T)   [Euler product = average over all classes].
Z_0 has the functional equation (numerator palindromic under T -> 1/(qT)) but no Euler product, and its zeros
are real for q >= 4: off the RH circle |T| = q^{-1/2}. The full zeta has them on the circle (Hasse).
Pencil N_lam(T) = 1 - qT + qT^2 + lam (h-1) T: on the circle iff |q - lam (h-1)| <= 2 sqrt q; at the window edge the
two circle zeros collide at T = +-q^{-1/2} and leave (cf. 14F collision lambda* in (0.1400, 0.1425)).
Brute-force check of the class counts for deg 1 and deg 2 over small primes."""
import itertools, json, math, cmath

def run(p, A_, B_):
    # F_{p^2} = F_p[r]/(r^2 - nr), nr a non-residue; elements (x0, x1)
    nr = next(x for x in range(2, p) if pow(x, (p-1)//2, p) == p-1)
    add = lambda u, v: ((u[0]+v[0]) % p, (u[1]+v[1]) % p)
    sub = lambda u, v: ((u[0]-v[0]) % p, (u[1]-v[1]) % p)
    mul = lambda u, v: ((u[0]*v[0] + nr*u[1]*v[1]) % p, (u[0]*v[1] + u[1]*v[0]) % p)
    def inv(u):
        n = (u[0]*u[0] - nr*u[1]*u[1]) % p; ni = pow(n, p-2, p); return ((u[0]*ni) % p, (-u[1]*ni) % p)
    A, B = (A_ % p, 0), (B_ % p, 0)
    els = [(x, y) for x in range(p) for y in range(p)]
    pts = [(x, y) for x in els for y in els if mul(y, y) == add(add(mul(mul(x, x), x), mul(A, x)), B)]
    O = None
    def eadd(P, Q):
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0] and add(P[1], Q[1]) == (0, 0): return None
        if P == Q: lam = mul(add(mul((3, 0), mul(P[0], P[0])), A), inv(mul((2, 0), P[1])))
        else: lam = mul(sub(Q[1], P[1]), inv(sub(Q[0], P[0])))
        x3 = sub(sub(mul(lam, lam), P[0]), Q[0]); y3 = sub(mul(lam, sub(P[0], x3)), P[1]); return (x3, y3)
    rat = [P for P in pts if P[0][1] == 0 and P[1][1] == 0]
    h = len(rat) + 1; a = p + 1 - h
    frob = lambda P: ((P[0][0], (-P[0][1]) % p), (P[1][0], (-P[1][1]) % p))
    classes = [None] + rat
    cnt1 = {str(c): 0 for c in classes}; cnt2 = {str(c): 0 for c in classes}
    for P in rat: cnt1[str(P)] += 1
    for P, Q in itertools.combinations_with_replacement(rat, 2): cnt2[str(eadd(P, Q))] += 1
    seen = set()
    for R in pts:
        if R in rat or R in seen: continue
        Rs = frob(R); seen.add(R); seen.add(Rs); cnt2[str(eadd(R, Rs))] += 1
    ok1 = all(cnt1[str(c)] == (0 if c is None else 1) for c in classes)
    ok2 = all(v == p for v in cnt2.values())
    z0 = sorted([abs(r) for r in [(p + cmath.sqrt(p*p - 4*p))/(2*p), (p - cmath.sqrt(p*p - 4*p))/(2*p)]])
    lo, hi = (p - 2*math.sqrt(p))/(h - 1), (p + 2*math.sqrt(p))/(h - 1)
    return dict(p=p, curve=f"y^2=x^3+{A_}x+{B_}", h=h, a=a, hasse_ok=a*a <= 4*p,
                deg1_counts_ok=ok1, deg2_counts_ok=ok2, Z0_zero_moduli=z0, RH_circle=p**-0.5,
                pencil_window_lam=[lo, hi], lam1_in_window=lo <= 1 <= hi)
def nonsing(p, A, B): return (4*A**3 + 27*B**2) % p != 0
out = [run(p, A, B) for p, A, B in [(5, 1, 1), (7, 1, 3), (11, 1, 1), (13, 1, 6), (17, 2, 3), (19, 1, 4)] if nonsing(p, A, B)]
for o in out: print(o)
json.dump(out, open('ff_epstein.json', 'w'), indent=1)
