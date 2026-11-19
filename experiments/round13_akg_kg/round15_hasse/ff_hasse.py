"""Round 15 part (d): the finite-field Hasse Gram, side by side.
On C x C (C/F_q genus g): Gamma_j = graph of F^j, H = C x pt, V = pt x C; primitive part Gamma'_j = Gamma_j - H - q^j V.
Intersection theory (PROVED, standard):  Gamma'_j . Gamma'_k = -q^j a_{k-j} (j<=k),  a_m = sum_i alpha_i^m = q^m + 1 - N_m,  a_0 = 2g.
Hasse Gram  T[j,k] := -Gamma'_j.Gamma'_k / q^{(j+k)/2} = sum_i omega_i^{k-j},  omega_i = alpha_i/sqrt q    (Toeplitz in k-j).
Hodge index (primitive part negative definite) => T >= 0 for every K.  For g = 1 the 2x2 block is [[2, a/sqrt q],[a/sqrt q, 2]] and
  (m, -n sqrt q) T (m, -n sqrt q)^T = 2(m^2 - a m n + q n^2) = 2 deg(m - nF)    (checked symbolically and numerically below).
Rows: (1) elliptic curves of probe B (zeta: Euler product)   (2) the class-0 partial zeta Z_0 (probe B 'Epstein', no Euler product)
(3) genus-2 curves y^2 = f(x) by brute-force counts over F_p, F_p^2   (4) a synthetic genus-2 'fake' with FE but omega off the circle,
whose 2x2 Hasse block PASSES while larger Toeplitz blocks fail: one Frobenius-shift block is not enough once g >= 2."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import _paths  # noqa: E402  (run from anywhere)
import json, math, cmath, itertools, numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'clues'))

def toeplitz(om, K):
    c = [sum(w**m for w in om) for m in range(K+1)]
    T = np.array([[c[abs(k-j)] if k >= j else np.conj(c[abs(k-j)]) for k in range(K+1)] for j in range(K+1)])
    T = 0.5*(T + T.conj().T); ev = np.linalg.eigvalsh(T)
    return [float(x) for x in ev], int(np.sum(ev < -1e-9*max(1, abs(ev).max())))

def count_points(p, f):
    """#C(F_p), #C(F_p^2) for y^2 = f(x), deg f = 5 (one point at infinity)."""
    nr = next(x for x in range(2, p) if pow(x, (p-1)//2, p) == p-1)
    mul = lambda u, v: ((u[0]*v[0] + nr*u[1]*v[1]) % p, (u[0]*v[1] + u[1]*v[0]) % p)
    def ev(x):
        r = (0, 0)
        for c in f: r = mul(r, x); r = ((r[0]+c) % p, r[1])
        return r
    sq = {}
    for y in itertools.product(range(p), repeat=2): v = mul(y, y); sq[v] = sq.get(v, 0) + 1
    N1 = 1 + sum(1 + (1 if pow(ev((x, 0))[0], (p-1)//2, p) == 1 else (0 if ev((x, 0))[0] == 0 else -1)) for x in range(p))
    N2 = 1 + sum(sq.get(ev(x), 0) for x in itertools.product(range(p), repeat=2))
    return N1, N2
def nonsing(p, f):
    # gcd(f, f') = 1 over F_p (polys as coefficient lists, highest first)
    def trim(a):
        while a and a[0] % p == 0: a = a[1:]
        return a
    def mod(a, b):
        a = trim([x % p for x in a]); b = trim([x % p for x in b]); inv = pow(b[0], p-2, p)
        while len(a) >= len(b):
            c = a[0]*inv % p; a = trim([(x - c*y) % p for x, y in zip(a, b + [0]*(len(a)-len(b)))])
        return a
    a, b = f[:], [c*(len(f)-1-i) % p for i, c in enumerate(f[:-1])]
    while trim(b): a, b = b, mod(a, b)
    return len(trim(a)) == 1

out = dict(elliptic=[], Z0=[], genus2=[], fake_genus2=[])
for r in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'clues', 'ff_epstein.json'))):
    q, a = r['p'], r['a']; al = (a + cmath.sqrt(a*a - 4*q))/2; om = [al/math.sqrt(q), al.conjugate()/math.sqrt(q)]
    T2 = np.array([[2, a/math.sqrt(q)], [a/math.sqrt(q), 2]])
    degok = all(abs(np.array([m, -n*math.sqrt(q)])@T2@np.array([m, -n*math.sqrt(q)]) - 2*(m*m - a*m*n + q*n*n)) < 1e-9
                for m in range(-4, 5) for n in range(-4, 5))
    # deg(1 - F^k) = N_k = q^k + 1 - (alpha^k + conj): positive for all k (finite check)
    Nk = [round((q**k + 1 - 2*(al**k).real)) for k in range(1, 7)]
    ev2, g2 = toeplitz(om, 1); ev8, g8 = toeplitz(om, 8)
    out['elliptic'].append(dict(q=q, curve=r['curve'], a=a, hasse_2x2_eig=ev2, ghosts_K1=g2, toeplitz_K8_rank=int(np.sum(np.array(ev8) > 1e-9)),
                                ghosts_K8=g8, deg_form_identity=degok, N_k=Nk))
    # Z_0: numerator 1 - qT + qT^2 = (1 - b1 T)(1 - b2 T), b1 + b2 = q, b1 b2 = q
    b1 = (q + cmath.sqrt(q*q - 4*q))/2; omz = [b1/math.sqrt(q), (q - b1)/math.sqrt(q)]
    ev2z, g2z = toeplitz(omz, 1); ev8z, g8z = toeplitz(omz, 8)
    out['Z0'].append(dict(q=q, hasse_offdiag=math.sqrt(q), hasse_2x2_eig=ev2z, ghosts_K1=g2z, ghosts_K8=g8z, omega=[abs(x) for x in omz]))
    print('E', q, a, ev2, g8, degok, '| Z0', ev2z, g2z, g8z)
for p, f in [(5, [1, 0, 0, 0, 1, 1]), (7, [1, 0, 0, 0, 1, 3]), (11, [1, 0, 2, 0, 1, 1]), (13, [1, 1, 0, 0, 0, 5])]:
    if not nonsing(p, f): print('singular', p, f); continue
    N1, N2 = count_points(p, f); a1 = p + 1 - N1                     # s1 = sum alpha = a1
    s2 = p*p + 1 - N2                                               # sum alpha^2
    e2 = (a1*a1 - s2)/2                                             # second elementary symmetric
    P = [1, -a1, e2, -p*a1, p*p]                                    # reciprocal char poly x^4 - s1 x^3 + e2 x^2 - p s1 x + p^2
    al = np.roots(P); om = list(al/math.sqrt(p))
    assert abs(p + 1 - N1 - sum(al).real) < 1e-6 and abs(p*p + 1 - N2 - sum(al**2).real) < 1e-6
    ev1, g1 = toeplitz(om, 1); ev10, g10 = toeplitz(om, 10)
    out['genus2'].append(dict(p=p, f=f, N1=N1, N2=N2, abs_omega=[float(abs(x)) for x in om], ghosts_K1=g1, ghosts_K10=g10,
                              lam_min_K10=ev10[0], rank_K10=int(np.sum(np.array(ev10) > 1e-8))))
    print('G2', p, N1, N2, [round(abs(x), 6) for x in om], g1, g10)
for phi, r in [(2.2, 2.0), (2.6, 1.5), (1.0, 1.2)]:
    om = [cmath.exp(1j*phi), cmath.exp(-1j*phi), r, 1/r]
    first = None; lst = []
    for K in range(1, 13):
        ev, g = toeplitz(om, K); lst.append(dict(K=K, lam_min=ev[0], ghosts=g))
        if g and first is None: first = K
    ck = [abs(sum(w**k for w in om)) for k in range(1, 40)]
    first_k = next((k+1 for k, c in enumerate(ck) if c > 4 + 1e-12), None)   # one-shift block [Delta, Gamma_{F^k}] fails: |c_k| > 2g
    print('   one-shift 2x2 blocks (Delta, F^k): first failing k =', first_k)
    out['fake_genus2'].append(dict(one_shift_first_fail_k=first_k, note='first ghost at K+1 = 2g = 4: the Toeplitz block must resolve all 2g modes', phi=phi, r=r, c1=float((2*math.cos(phi) + r + 1/r)), passes_2x2=lst[0]['ghosts'] == 0, first_ghost_K=first, rows=lst))
    print('fake g2', phi, r, 'c1=%.3f' % (2*math.cos(phi)+r+1/r), '2x2 ghosts', lst[0]['ghosts'], 'first ghost K', first)
json.dump(out, open('ff_hasse.json', 'w'), indent=1, default=float)
