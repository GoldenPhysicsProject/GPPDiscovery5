"""Round 31, task 3: does Daniel's doubling D(z) = 2z carry arithmetic up the ladder, or only the stage?
 L1  CHECK: the d=1 -> d=2 lift is the substitution r -> r^2 on R+; Mellin: M[f(r^2)](Delta) = (1/2) M[f](Delta/2).
     Applied to the theta kernel: the rung-1 transform of Theta(r^2)-1 is (1/2)*Mellin_0 at s = Delta/2, i.e. Delta = 2s,
     critical line Re s = 1/2 -> celestial principal series Re Delta = 1. Dirichlet side: n^{-s} = (sqrt n)^{-Delta}:
     the same prime torus with circle lengths (log p)/2 -- a rescaled flow on an isomorphic torus.
 L2  CHECK: Gr(2,4)(F_q) by brute force (RREF enumeration) = 1+q+2q^2+q^3+q^4 for q = 2,3,5, and the Euler factor at p
     prod_k (1 - p^{k-s})^{-m_k}, (m_k) = (1,1,2,1,1): per-prime 'Satake radii' p^{k}, normalised by the centre p^2:
     p^{-2}, p^{-1}, 1, 1, p, p^2  -> stretched exactly like F_theta with theta in {1, 2}, plus multiplicity 2 in the middle.
 L3  F_theta lifted: zero lines Re s = 1/2 +- theta go to Re Delta = 1 +- 2 theta (still off the lifted line)."""
import json, math, itertools, numpy as np, mpmath
mpmath.mp.dps = 30
out = {}
# L1
def Theta(t):
    # Jacobi inversion Theta(t) = t^{-1/2} Theta(1/t) for small t (keeps the nome small)
    if t < 1: return t**-0.5*mpmath.jtheta(3, 0, mpmath.exp(-mpmath.pi/t))
    return mpmath.jtheta(3, 0, mpmath.exp(-mpmath.pi*t))
def mellin(f, z): return mpmath.quad(lambda r: f(r)*r**(z - 1), [0, 1, mpmath.inf])
rows = []
for Delta in (mpmath.mpc(3, 0.7), mpmath.mpc(2.6, -4), mpmath.mpc(4.2, 9)):
    # rung 0 (d=1): g(r) = Theta(r^2) - 1 (Gaussian lattice sum), M[g](s) = pi^{-s/2} Gamma(s/2) zeta(s) = Lambda(s)
    # lift (d=2):    g(r^2) = Theta(r^4) - 1,  M[g(r^2)](Delta) = (1/2) Lambda(Delta/2)
    lhs = mellin(lambda r: Theta(r**4) - 1, Delta)
    rhs = mellin(lambda r: Theta(r**2) - 1, Delta/2)/2
    xi_side = mpmath.pi**(-Delta/4)*mpmath.gamma(Delta/4)*mpmath.zeta(Delta/2)
    rows.append(dict(Delta=str(Delta), rung1=complex(lhs), rung0_half=complex(rhs), xi_form=complex(xi_side/2),
                     rel_err=float(abs(lhs - rhs)/abs(rhs)), rel_err_xi=float(abs(lhs - xi_side/2)/abs(rhs))))
out['L1_mellin_lift'] = rows; print(rows)
# L2
def gr24_count(q):
    # 2-dim subspaces of F_q^4 = 2x4 matrices in reduced row echelon form of rank 2
    cnt = 0
    for piv in itertools.combinations(range(4), 2):
        free = [(i, j) for i in range(2) for j in range(4) if j > piv[i] and j not in piv]
        cnt += q**len(free)
    return cnt
L2 = [(q, gr24_count(q), 1 + q + 2*q*q + q**3 + q**4) for q in (2, 3, 4, 5, 7)]
out['L2_point_counts'] = L2; print(L2)
# Euler factor check at p=2: log of prod_k (1-p^{k-s})^{-m_k} vs sum_r N(p^r) p^{-rs}/r  (Hasse-Weil, q = p^r), at s = 7.3
p, s = 2, mpmath.mpf('7.3'); m = [1, 1, 2, 1, 1]
lhs = -sum(mk*mpmath.log(1 - mpmath.mpf(p)**(k - s)) for k, mk in enumerate(m))
rhs = mpmath.nsum(lambda r: gr24_count(int(p**int(r))) * mpmath.mpf(p)**(-r*s)/r if r < 40 else (1 + p**r + 2*p**(2*r) + p**(3*r) + p**(4*r))*mpmath.mpf(p)**(-r*s)/r, [1, mpmath.inf])
out['L2_euler_factor_p2'] = dict(lhs=float(lhs), rhs=float(rhs), err=float(abs(lhs - rhs)))
out['L2_normalised_radii'] = 'p^{-2}, p^{-1}, 1, 1, p, p^2'
# zero strips of zeta(s-k): k < Re s < k+1; FE line of Lambda_Gr is Re s = 5/2
out['L2_zero_strips'] = [dict(k=k, strip=(k, k + 1), mult=m[k], touches_FE_line_5_2=(k < 2.5 < k + 1)) for k in range(5)]
# L4: the integer lattices of the division algebras (R, C, H, O) and their Dirichlet series (classical):
#   Z      : #{x in Z : x^2 = n}            -> 2 zeta(2w)            (theta weight 1/2)
#   Z[i]   : r_2(n) = 4 sum_{d|n} chi_{-4}(d) -> 4 zeta(w) L(w,chi_-4)   (weight 1)  -- two species, SAME line
#   Hurwitz: 24 sum_{d|n, d odd} d          -> 24 (1-2^{1-w}) zeta(w) zeta(w-1)   (weight 2)  -- stretch p^{1/2}
#   E8 (Coxeter-Dickson octaves, norm 2n): 240 sigma_3(n) -> 240 zeta(w) zeta(w-3)  (weight 4)  -- stretch p^{3/2}
# zeta(w) zeta(w-k) = zeta(u + k/2) zeta(u - k/2), u = w - k/2: exactly F_theta with theta = k/2.
def count_Zi(n):
    M = int(math.isqrt(n)) + 1
    return sum(1 for a in range(-M, M + 1) for b in range(-M, M + 1) if a*a + b*b == n)
def count_hurwitz(n):
    M = int(math.isqrt(n)) + 1; c = 0
    for x in itertools.product(range(-M, M + 1), repeat=4):
        if sum(v*v for v in x) == n: c += 1
    for x in itertools.product(range(-M, M + 1), repeat=4):     # half-integer points (2y+1)/2
        if sum((2*v + 1)**2 for v in x) == 4*n: c += 1
    return c
def count_E8(n):
    # E8 = {x in Z^8 or (Z+1/2)^8 : sum x even}, norm sum x^2 = 2n; count by recursion on coordinates
    from functools import lru_cache
    tgt = 8*n  # work with 2x (integers: even; half-integers: odd), sum (2x)^2 = 4*2n
    def cnt(parity):
        vals = [v for v in range(-int(math.isqrt(tgt)) - 1, int(math.isqrt(tgt)) + 2) if v % 2 == parity]
        @lru_cache(None)
        def f(i, rem, ssum):
            if i == 8: return 1 if rem == 0 and ssum % 4 == 0 else 0
            return sum(f(i + 1, rem - v*v, (ssum + v) % 4) for v in vals if v*v <= rem)
        return f(0, tgt, 0)
    return cnt(0) + cnt(1)
def chi4(d): return [0, 1, 0, -1][d % 4]
def divs(n): return [d for d in range(1, n + 1) if n % d == 0]
L4 = []
for n in range(1, 9):
    L4.append(dict(n=n, Zi=(count_Zi(n), 4*sum(chi4(d) for d in divs(n))),
                   Hurwitz=(count_hurwitz(n), 24*sum(d for d in divs(n) if d % 2)),
                   E8=(count_E8(n), 240*sum(d**3 for d in divs(n)))))
out['L4_division_lattices'] = L4; print(L4)
out['L4_stretch'] = dict(R=dict(weight=0.5, L='zeta', theta=0), C=dict(weight=1, L='zeta*L(chi_-4)', theta=0),
                         H=dict(weight=2, L='zeta(w)zeta(w-1)', theta=0.5), O=dict(weight=4, L='zeta(w)zeta(w-3)', theta=1.5))
# L3
TH = 0.25
out['L3_Ftheta_lift'] = dict(rung0_lines=(0.5 - TH, 0.5 + TH), rung1_lines=(1 - 2*TH, 1 + 2*TH), rung1_FE_line=1.0)
print(out['L2_euler_factor_p2'], out['L2_zero_strips'])
json.dump(out, open('r31_lift.json', 'w'), indent=1, default=str)
