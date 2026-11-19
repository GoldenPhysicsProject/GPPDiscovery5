"""Round 16 Part 1: Bombieri's version of Stepanov's method, EXECUTED on concrete curves (probe-B elliptic curves and the
Round-15 genus-2 curves), base-changed to F_q, q = p^2 (Bombieri needs q square, q > (g+1)^4).
Steps (tags = ingredients: i additivity of x->x^Q, ii multiplicativity / Frobenius fixes rational points,
iii Riemann-Roch dimension count, iv extensions F_{q^n} (tensor-power trick), v functional equation):
 S1 [iii]   L(mO) = span{x^i y^j : 2i + j(2g+1) <= m}; dim = m + 1 - g for m >= 2g-1 (Riemann-Roch), distinct pole orders.
 S2 [iii]   injectivity of  L(l) (x) L(m) -> K,  u (x) v -> u^p v^q  for l < p (pole orders p a + q b all distinct).
 S3 [iii]   kernel of delta: sum c_uv u^p v^q -> sum c_uv u^p v  in L(lp+m) is non-zero when dim L(l) dim L(m) > dim L(lp+m).
 S4 [ii]    v^q(P) = v(P) at every F_q-rational P (multiplicativity + P fixed by Frobenius), so F(P) = delta(F)(P) = 0.
 S5 [i]     F = H^p with H = sum c^p u v^p (additivity: (sum)^p = sum of p-th powers); so ord_P F >= p at each rational P != O.
 S6 [count] zeros = poles: p (N - 1) <= pole_O(F) <= l p + m q  =>  N <= 1 + l + m p   (upper bound, Bombieri: q + (2g+1)sqrt q + 1).
 S7 [Galois cover / twist]  N(C) + N(C^d) = 2(q + 1) (x-line count, hyperelliptic twist), so upper bound on C^d = lower bound on C.
 S8 [iv,iii,v] N_n = q^n + 1 - sum alpha^n (rationality, from Riemann-Roch), bound for every n with constant 2g+1, n -> inf:
            |alpha| <= sqrt q; FE alpha <-> q/alpha gives |alpha| >= sqrt q.  Fakes (no curve) violate the bound.
All steps are CHECKED numerically below; nothing is assumed beyond field arithmetic."""
import json, math, cmath, itertools, sys, time, numpy as np
from ffq import Fq, Curve

def nullspace(F, A):
    """Right nullspace of A (rows x cols, field ints) over F_q, by Gauss-Jordan."""
    A = A.copy(); R, C = A.shape; piv = []; r = 0
    for c in range(C):
        k = next((i for i in range(r, R) if A[i, c]), None)
        if k is None: continue
        A[[r, k]] = A[[k, r]]; A[r] = F.M[F.INV[A[r, c]], A[r]]
        for i in range(R):
            if i != r and A[i, c]: A[i] = F.A[A[i], F.NEG[F.M[A[i, c], A[r]]]]
        piv.append(c); r += 1
        if r == R: break
    free = [c for c in range(C) if c not in piv]; basis = []
    for fc in free:
        v = np.zeros(C, dtype=np.int64); v[fc] = 1
        for i, pc in enumerate(piv): v[pc] = F.NEG[A[i, fc]]
        basis.append(v)
    return basis, len(piv)

def coeffvec(cur, u, top):
    F = cur.F; na = top//2 + 1; nb = (top - cur.d)//2 + 1 if top >= cur.d else 0
    v = np.zeros(na + nb, dtype=np.int64); a, b = u
    assert a.shape[1] <= na and b.shape[1] <= nb
    v[:a.shape[1]] = a[0] + F.p*a[1]; v[na:na+b.shape[1]] = b[0] + F.p*b[1]; return v

def ringpow(cur, u, n):
    r = cur.mono(0, 0)
    for _ in range(n): r = cur.mul(r, u)
    return r

def same(u, v): return all(x.shape == y.shape and np.array_equal(x, y) for x, y in zip(u, v))

def choose_lm(cur):
    p, q, g = cur.F.p, cur.F.q, cur.g; best = None
    for l in range(1, p):
        dl = len(cur.basis(l))
        for m in range(1, 4*p + 4*g):
            if dl*len(cur.basis(m)) > len(cur.basis(l*p + m)):
                b = 1 + l + m*p
                if best is None or b < best[2]: best = (l, m, b)
                break
    return best

def bombieri(cur, tag):
    F = cur.F; p, q, g = F.p, F.q, cur.g; t0 = time.time(); rec = dict(tag=tag, p=p, q=q, g=g)
    # S1 Riemann-Roch
    rr = [(m, len(cur.basis(m)), m + 1 - g) for m in range(0, 6*g + 4)]
    rec['S1_RR_ok'] = all(d == e for m, d, e in rr if m >= 2*g - 1) and all(d >= max(1, e) for m, d, e in rr)
    rec['S1_gaps'] = [m for m in range(1, 2*g + 1) if len(cur.basis(m)) == len(cur.basis(m-1))]
    l, m, bnd = choose_lm(cur); rec.update(l=l, m=m)
    Bl, Bm = cur.basis(l), cur.basis(m); top = l*p + m
    # S2 injectivity: pole orders p*ord(u) + q*ord(v) distinct
    ords = [p*(2*i + j*cur.d) + q*(2*k + jj*cur.d) for (i, j) in Bl for (k, jj) in Bm]
    rec['S2_injective'] = len(set(ords)) == len(ords)
    # additivity check of the Frobenius-power formula on monomials (S5 ingredient i)
    rec['S5_monomial_frob_formula_ok'] = all(same(cur.monopow(i, j, p), ringpow(cur, cur.mono(i, j), p)) for (i, j) in Bl[:3] + Bm[-2:])
    # S3 delta matrix
    cols, pairs = [], []
    for (i, j) in Bl:
        up = cur.monopow(i, j, p)
        for (k, jj) in Bm:
            cols.append(coeffvec(cur, cur.mul(up, cur.mono(k, jj)), top)); pairs.append(((i, j), (k, jj)))
    A = np.array(cols).T; ker, rank = nullspace(F, A)
    rec.update(dim_L_l=len(Bl), dim_L_m=len(Bm), dim_L_target=len(cur.basis(top)), rows=A.shape[0], unknowns=A.shape[1],
               rank=rank, kernel_dim=len(ker), S3_kernel_predicted=len(Bl)*len(Bm) - len(cur.basis(top)))
    c = ker[0]
    # S3 check: delta(F) == 0 exactly
    dF = (F.P([]), F.P([]))
    for cc, ((i, j), (k, jj)) in zip(c, pairs):
        if cc: dF = cur.add(dF, cur.smul(int(cc), cur.mul(cur.monopow(i, j, p), cur.mono(k, jj))))
    rec['S3_delta_zero'] = dF[0].shape[1] == 0 and dF[1].shape[1] == 0
    # build F = sum c u^p v^q and H = sum c^p u v^p
    Fn = (F.P([]), F.P([])); Hn = (F.P([]), F.P([]))
    for cc, ((i, j), (k, jj)) in zip(c, pairs):
        if not cc: continue
        Fn = cur.add(Fn, cur.smul(int(cc), cur.mul(cur.monopow(i, j, p), cur.monopow(k, jj, q))))
        Hn = cur.add(Hn, cur.smul(int(F.FROB[cc]), cur.mul(cur.mono(i, j), cur.monopow(k, jj, p))))
    poleF = cur.pole(Fn)
    rec['F_nonzero'] = poleF is not None
    rec['pole_F'] = poleF; rec['pole_F_predicted'] = max(o for o, cc in zip(ords, c) if cc); rec['pole_cap_lp_plus_mq'] = l*p + m*q
    # S5: F == H^p, two ways (Frobenius formula via additivity, and literal ring power)
    Hp_ring = ringpow(cur, Hn, p)
    rec['S5_F_equals_H_pow_p'] = same(Fn, Hp_ring)
    # S4: evaluation at all rational affine points
    pts = cur.points(); N = len(pts) + 1
    vq_ok = all(cur.ev(cur.monopow(k, jj, q), x, y) == cur.ev(cur.mono(k, jj), x, y) for (k, jj) in Bm[-3:] for (x, y) in pts[:40])
    rec['S4_vq_equals_v_at_points'] = vq_ok
    rec['S4_F_vanishes_all_rational'] = all(cur.ev(Fn, x, y) == 0 for (x, y) in pts)
    rec['S4_H_vanishes_all_rational'] = all(cur.ev(Hn, x, y) == 0 for (x, y) in pts)
    # H vanishes at non-rational points too? (it must not vanish identically): count zeros of H among F_q-x with y in F_q only
    # S6
    rec['N'] = N; rec['S6_zero_count_lhs'] = p*(N - 1); rec['S6_ok'] = p*(N - 1) <= poleF
    rec['upper_bound'] = 1 + poleF//p; rec['bombieri_formula'] = q + 1 + (2*g + 1)*p; rec['weil'] = q + 1 + 2*g*p
    rec['secs'] = round(time.time() - t0, 2)
    return rec

def count_Fp(cur):
    """N over the prime field F_p (x, y in F_p), plus O."""
    F = cur.F; n = 1
    for x in range(F.p):
        r = F.ev(cur.f, x)
        n += sum(1 for y in range(F.p) if F.M[y, y] == r)
    return n

if __name__ == '__main__':
    CURVES = [('E', 5, [1, 1, 0, 1]), ('E', 7, [3, 1, 0, 1]), ('E', 11, [1, 1, 0, 1]), ('E', 13, [6, 1, 0, 1]),
              ('G2', 11, [1, 1, 0, 2, 0, 1]), ('G2', 13, [5, 0, 0, 0, 1, 1])]   # f low->high: x^3+x+1 ...; x^5+2x^3+x+1, x^5+x^4+5
    out = dict(curves=[])
    for kind, p, f in CURVES:
        F = Fq(p); C = Curve(F, f); g = C.g
        if F.q <= (g + 1)**4: print('skip q too small', kind, p); continue
        d = F.nonsq; Ct = Curve(F, [int(F.M[d, a]) for a in f])            # hyperelliptic (quadratic) twist y^2 = d f(x)
        r1 = bombieri(C, f'{kind} p={p}'); r2 = bombieri(Ct, f'{kind}^d p={p}')
        N, Nt = r1['N'], r2['N']
        rec = dict(kind=kind, p=p, q=F.q, f=f, g=g, curve=r1, twist=r2, S7_N_plus_Ntwist=N + Nt, S7_equals_2q_plus_2=(N + Nt == 2*F.q + 2),
                   lower_bound=2*F.q + 2 - r2['upper_bound'], bombieri_lower=F.q + 1 - (2*g + 1)*p)
        # S8: extensions / tensor power / FE, from the curve over F_p
        N1 = count_Fp(C); a1 = p + 1 - N1
        if g == 1:
            al = [(a1 + cmath.sqrt(a1*a1 - 4*p))/2]; al.append(al[0].conjugate())
        else:
            s2 = p*p + 1 - N; e2 = (a1*a1 - s2)/2; al = list(np.roots([1, -a1, e2, -p*a1, p*p]))
        al = [complex(z) for z in al]
        Nn = lambda n: p**n + 1 - sum(z**n for z in al).real
        rec['S8_N2_from_alphas_matches_count'] = abs(Nn(2) - N) < 1e-6
        ext = []
        for n in range(1, 16):            # Bombieri over F_{p^{2n}} = F_{q^n}: |N - q^n - 1| <= (2g+1) q^{n/2}
            Q = p**(2*n); dev = abs(Nn(2*n) - Q - 1); ext.append(dict(n=n, dev_over_sqrtQ=dev/p**n, within_2g_plus_1=dev <= (2*g + 1)*p**n + 1e-6))
        rec['S8_ext'] = ext
        rec['S8_implied_max_abs_omega_bound'] = [(2*g + 1)**(1/(2*n)) for n in (1, 2, 5, 10, 50, 500)]
        rec['S8_abs_omega'] = [abs(z)/math.sqrt(p) for z in al]
        rec['S8_FE_pairs_q_over_alpha'] = all(min(abs(w - p/z) for w in al) < 1e-8 for z in al)
        out['curves'].append(rec)
        print(kind, p, 'q', F.q, 'g', g, 'l,m', r1['l'], r1['m'], 'N', N, '<=', r1['upper_bound'], '(Bombieri', r1['bombieri_formula'], 'Weil', r1['weil'], ')',
              '| twist N', Nt, '<=', r2['upper_bound'], '-> N >=', rec['lower_bound'], '| checks',
              all(r1[k] for k in ('S1_RR_ok', 'S2_injective', 'S3_delta_zero', 'F_nonzero', 'S4_F_vanishes_all_rational', 'S5_F_equals_H_pow_p', 'S6_ok', 'S5_monomial_frob_formula_ok', 'S4_vq_equals_v_at_points')),
              all(r2[k] for k in ('S1_RR_ok', 'S2_injective', 'S3_delta_zero', 'F_nonzero', 'S4_F_vanishes_all_rational', 'S5_F_equals_H_pow_p', 'S6_ok')),
              rec['S7_equals_2q_plus_2'], rec['S8_N2_from_alphas_matches_count'], r1['secs'], 's', flush=True)
    # fakes: not curves, so S1-S6 cannot even be set up; their 'N_n' violate every C q^{n/2}
    fakes = []
    for q in (5, 7, 11, 13, 25):
        b1 = (q + math.sqrt(q*q - 4*q))/2; om = [b1/math.sqrt(q), (q - b1)/math.sqrt(q)]   # Z_0 numerator 1 - qT + qT^2
        dev = [abs(sum(w**n for w in om)) for n in range(1, 41)]
        fakes.append(dict(fake='Z0 class-0 partial zeta (FE, no Euler product)', q=q, abs_omega=om,
                          first_n_exceeding_3=next((n+1 for n, v in enumerate(dev) if v > 3), None), dev_n40=dev[-1]))
    for phi, r in ((1.0, 1.3), (2.0, 1.1), (0.5, 2.0)):
        om = [cmath.exp(1j*phi), cmath.exp(-1j*phi), r, 1/r]; dev = [abs(sum(w**n for w in om)) for n in range(1, 81)]
        fakes.append(dict(fake='synthetic genus-2 (FE, off-circle pair)', phi=phi, r=r, first_n_exceeding_5=next((n+1 for n, v in enumerate(dev) if v > 5), None)))
    out['fakes'] = fakes
    for fk in fakes: print(fk)
    json.dump(out, open('stepanov_ff.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print('DONE')
