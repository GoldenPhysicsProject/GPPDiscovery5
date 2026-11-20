"""Round 29 (B): toy run of the sqrt(q) multiplicity trade-off, Stepanov's method for y^2 = f(x) over F_p.
Classical source (PROVED): S. A. Stepanov (1969); Iwaniec-Kowalski, Analytic Number Theory, Sec. 11.6, Prop. 11.18 and
Thm 11.13; Bombieri (1973). Here we EXECUTE the construction (Part 1, CHECKED) and evaluate the exact dimension-count bound
for large q (Part 2, NUMERICAL), and compare with true point counts (Part 3, CHECKED).

Part 1, for small primes p and every multiplicity l, by exact linear algebra over F_p:
  r = f^l * sum_{j<J} (r_j + s_j g) X^{jp},  g = f^c, c = (p-1)/2,  deg r_j, s_j <= c - m   (I-K (11.23))
  impose sigma^(k) = 0 for k < l (I-K (11.25)): then E^k r(x) = 0 at every x in S_a = {f(x)=0 or g(x)=a}.
  Find the smallest J with a nonzero solution, build r, CHECK r != 0, CHECK ord_x(r) >= l at every x in S_a, and record
     ASK  := sum_{x in S_a} ord_x(r)     (forced divisibility, >= l |S_a|)
     MAIN := deg r                        (size)
  The 'product formula' for the nonzero polynomial r (sum over all places of ord = 0, the place at infinity carries -deg r)
  gives ASK <= MAIN. This is the function-field analogue of |N| >= 1 and of Alpoge's ASK - Main <= error.
Part 2: bound B_a(l) = floor(deg-bound / l) with the smallest J allowed by the I-K equation count A > B; error bound
  E(l) = 2 B(l) - q on |N_aff - q|; optimum l* ~ sqrt(q).
Part 3: true affine point counts N_aff = sum_x (1 + legendre(f(x))) vs E(l*) and the Hasse-Weil 2g sqrt(q).
Run:  python3 toy_fq.py   (about a minute; writes toy_fq.json and toy_fq.png)."""
import json, math, time
import numpy as np

def trim(a):
    a = np.asarray(a, dtype=np.int64); n = len(a)
    while n > 1 and a[n-1] == 0: n -= 1
    return a[:n].copy()

def pmul(a, b, p):
    # exact convolution mod p (split to avoid int64 overflow)
    a = np.asarray(a, dtype=np.int64) % p; b = np.asarray(b, dtype=np.int64) % p
    if len(a) * p * p < 2**62:
        return trim(np.convolve(a, b) % p)
    out = np.zeros(len(a) + len(b) - 1, dtype=object)
    for i, ai in enumerate(a.tolist()):
        if ai: out[i:i+len(b)] += ai * b.astype(object)
    return trim(np.array([int(v) % p for v in out], dtype=np.int64))

def ppow(a, e, p):
    r = np.array([1], dtype=np.int64); b = trim(a)
    while e:
        if e & 1: r = pmul(r, b, p)
        b = pmul(b, b, p); e >>= 1
    return r

def pdivmod(a, b, p):
    a = trim(a) % p; b = trim(b); db = len(b) - 1; inv = pow(int(b[-1]), p - 2, p)
    if len(a) - 1 < db: return np.array([0], dtype=np.int64), a
    q = np.zeros(len(a) - db, dtype=np.int64); a = a.copy()
    for i in range(len(a) - 1, db - 1, -1):
        c = (int(a[i]) * inv) % p
        if c:
            q[i - db] = c; a[i - db:i + 1] = (a[i - db:i + 1] - c * b) % p
    return trim(q), trim(a[:max(db, 1)])

def hasse(a, k, p):
    # k-th Hasse derivative: coefficient of X^{n-k} is binom(n, k) a_n  (mod p)
    a = trim(a); n = len(a)
    if k >= n: return np.array([0], dtype=np.int64)
    out = [(math.comb(i, k) % p) * int(a[i]) % p for i in range(k, n)]
    return trim(np.array(out, dtype=np.int64))

def peval(a, x, p):
    r = 0
    for c in reversed(a.tolist()): r = (r * x + c) % p
    return r

def nullvec(M, p):
    """One nonzero vector in the right nullspace of M over F_p, or None."""
    A = M.copy() % p; R, C = A.shape; piv = []; r = 0
    for c in range(C):
        rows = np.nonzero(A[r:, c])[0]
        if len(rows) == 0: continue
        k = r + rows[0]
        if k != r: A[[r, k]] = A[[k, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        nz = np.nonzero(A[:, c])[0]; nz = nz[nz != r]
        if len(nz): A[nz] = (A[nz] - np.outer(A[nz, c], A[r])) % p
        piv.append(c); r += 1
        if r == R: break
    free = [c for c in range(C) if c not in set(piv)]
    if not free: return None, len(piv)
    v = np.zeros(C, dtype=np.int64); fc = free[0]; v[fc] = 1
    for i, pc in enumerate(piv): v[pc] = (-A[i, fc]) % p
    return v, len(piv)

def order_at(r, x0, p):
    lin = np.array([(-x0) % p, 1], dtype=np.int64); k = 0; a = trim(r)
    while True:
        q, rem = pdivmod(a, lin, p)
        if np.any(rem % p): return k
        a = q; k += 1

def stepanov_build(p, f, a, l, Jmax=60):
    m = len(f) - 1; c = (p - 1) // 2; D0 = c - m
    fl = ppow(f, l, p); flc = ppow(f, l + c, p)
    # P[k][i] = E^k(f^l X^i) / f^(l-k),  Q[k][i] = E^k(f^(l+c) X^i) / f^(l+c-k)
    P, Q = [], []
    for k in range(l):
        dl, dlc = ppow(f, l - k, p), ppow(f, l + c - k, p); Pk, Qk = [], []
        for i in range(D0 + 1):
            Xi = np.zeros(i + 1, dtype=np.int64); Xi[i] = 1
            q1, r1 = pdivmod(hasse(pmul(fl, Xi, p), k, p), dl, p)
            q2, r2 = pdivmod(hasse(pmul(flc, Xi, p), k, p), dlc, p)
            assert not np.any(r1) and not np.any(r2)        # Lemma 11.20 (exact divisibility)
            Pk.append(q1); Qk.append(q2)
        P.append(Pk); Q.append(Qk)
    for J in range(1, Jmax + 1):
        # unknown index: (j, i, 0) -> r_j coeff i ; (j, i, 1) -> s_j coeff i
        cols = 2 * J * (D0 + 1); rowsk = [D0 + k * (m - 1) + J + 1 for k in range(l)]
        off = np.concatenate([[0], np.cumsum(rowsk)]); M = np.zeros((off[-1], cols), dtype=np.int64)
        for k in range(l):
            for j in range(J):
                for i in range(D0 + 1):
                    for t, poly in ((0, P[k][i]), (1, (a * Q[k][i]) % p)):
                        col = (j * (D0 + 1) + i) * 2 + t; n = len(poly)
                        assert j + n <= rowsk[k]
                        M[off[k] + j: off[k] + j + n, col] = (M[off[k] + j: off[k] + j + n, col] + poly) % p
        v, rank = nullvec(M, p)
        if v is not None: break
    else:
        return None
    g = ppow(f, c, p); r = np.array([0], dtype=np.int64)
    for j in range(J):
        rj = trim(np.array([v[(j * (D0 + 1) + i) * 2] for i in range(D0 + 1)]))
        sj = trim(np.array([v[(j * (D0 + 1) + i) * 2 + 1] for i in range(D0 + 1)]))
        term = (np.concatenate([rj, np.zeros(max(0, len(g) + len(sj) - 1 - len(rj)), dtype=np.int64)]) + 0)
        gs = pmul(sj, g, p); L = max(len(term), len(gs)); tt = np.zeros(L, dtype=np.int64)
        tt[:len(rj)] += rj; tt[:len(gs)] += gs; tt = trim(tt % p)
        shift = np.zeros(j * p + 1, dtype=np.int64); shift[-1] = 1
        acc = pmul(tt, shift, p); L = max(len(r), len(acc)); s = np.zeros(L, dtype=np.int64)
        s[:len(r)] += r; s[:len(acc)] += acc; r = trim(s % p)
    r = pmul(fl, r, p)
    return dict(J=J, rank=int(rank), unknowns=int(cols), equations=int(off[-1]), r=r)

def S_set(p, f, a):
    c = (p - 1) // 2; out = []
    for x in range(p):
        fx = peval(f, x, p)
        if fx == 0 or pow(fx, c, p) == a % p: out.append(x)
    return out

def count_aff(p, f):
    xs = np.arange(p, dtype=object); fx = np.zeros(p, dtype=object)
    for cf in reversed([int(t) for t in f]): fx = (fx * xs + cf) % p
    fx = np.array(fx, dtype=np.int64); c = (p - 1) // 2
    # vectorised modpow
    res = np.ones(p, dtype=np.int64); b = fx % p; e = c
    while e:
        if e & 1: res = (res * b) % p
        b = (b * b) % p; e >>= 1
    leg = np.where(fx == 0, 0, np.where(res == 1, 1, -1))
    return int(np.sum(1 + leg)), int(np.sum(fx == 0))

def ik_bound(q, m, l):
    """Exact I-K dimension count: smallest J with A > B; return deg bound (11.24) and B_a = floor(deg/l)."""
    c = (q - 1) // 2; D0 = c - m
    eqs0 = sum(D0 + k * (m - 1) for k in range(l))
    J = 1
    while 2 * (D0 + 1) * J <= eqs0 + l * J: J += 1
    deg = l * m + D0 + c * m + (J - 1) * q
    return J, deg, deg // l

def main():
    t0 = time.time(); out = dict(part1=[], part2=[], part3=[])
    curves = {3: [1, 2, 0, 1], 5: [3, 1, 0, 2, 0, 1]}            # f = x^3+2x+1 (genus 1), x^5+2x^3+x+3 (genus 2), low->high
    # ---------------- Part 1: explicit auxiliary polynomials (CHECKED) ----------------
    for p in (53, 101):
        for m, fc in curves.items():
            f = trim(np.array(fc) % p); Naff, N0 = count_aff(p, f)
            for a in (1, -1):
                S = S_set(p, f, a)
                for l in range(1, p // 8 + 1):
                    res = stepanov_build(p, f, a % p, l)
                    if res is None: continue
                    r = res['r']; deg = len(r) - 1
                    ords = [order_at(r, x, p) for x in S]
                    ok = bool(np.any(r)) and min(ords) >= l
                    ask = int(sum(ords))
                    rec = dict(p=p, m=m, a=a, l=l, J=res['J'], unknowns=res['unknowns'], equations=res['equations'],
                               rank=res['rank'], deg_r=deg, S=len(S), min_ord=min(ords), ASK=ask, MAIN=deg,
                               ASK_le_MAIN=ask <= deg, vanish_ok=ok, bound_S=deg // l, N_aff=Naff, N0=N0)
                    out['part1'].append(rec)
                    print('P1', {k: rec[k] for k in ('p', 'm', 'a', 'l', 'J', 'deg_r', 'S', 'min_ord', 'ASK', 'bound_S', 'vanish_ok')}, flush=True)
    # ---------------- Part 2: I-K exact bound vs l for large q (NUMERICAL) ----------------
    for q in (1009, 10007, 100003, 1000003):
        for m in (3, 5):
            for l in range(1, int(4 * math.sqrt(q)) + 2):
                if l > q // 8: break
                J, deg, B = ik_bound(q, m, l)
                out['part2'].append(dict(q=q, m=m, l=l, J=J, deg=deg, B=B, E=2 * B - q))
    # ---------------- Part 3: true counts vs best bound (CHECKED) ----------------
    for q in (1009, 10007, 100003, 1000003):
        for m, fc in curves.items():
            f = trim(np.array(fc) % q); Naff, N0 = count_aff(q, f)
            rows = [r for r in out['part2'] if r['q'] == q and r['m'] == m]
            best = min(rows, key=lambda r: r['E'])
            g = (m - 1) // 2
            out['part3'].append(dict(q=q, m=m, genus=g, N_aff=Naff, dev=Naff - q, abs_dev=abs(Naff - q),
                                     hasse=2 * g * math.sqrt(q), E_best=best['E'], l_best=best['l'],
                                     l_best_over_sqrtq=best['l'] / math.sqrt(q), E_over_sqrtq=best['E'] / math.sqrt(q),
                                     within=abs(Naff - q) <= best['E']))
            print('P3', out['part3'][-1], flush=True)
    out['elapsed_s'] = time.time() - t0
    json.dump(out, open('toy_fq.json', 'w'), indent=1, default=int)
    plot(out)

def plot(out):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    for q in (1009, 10007, 100003, 1000003):
        rows = [r for r in out['part2'] if r['q'] == q and r['m'] == 3]
        L = np.array([r['l'] for r in rows]); E = np.array([r['E'] for r in rows], dtype=float)
        ax[0].plot(L / math.sqrt(q), E / math.sqrt(q), label=f'q={q}')
    ax[0].set_ylim(0, 60); ax[0].set_xlabel('multiplicity l / sqrt(q)'); ax[0].set_ylabel('error bound E(l) / sqrt(q)')
    ax[0].set_title('Stepanov bound |N-q| <= E(l), genus 1 (I-K count)\ncollapse: E ~ q/l + l, optimum l ~ sqrt(q)'); ax[0].legend()
    for p, mk in ((53, 'o'), (101, 's')):
        rows = [r for r in out['part1'] if r['p'] == p and r['m'] == 3]
        Ls = sorted(set(r['l'] for r in rows)); E = []
        for l in Ls:
            B1 = [r['bound_S'] for r in rows if r['l'] == l and r['a'] == 1][0]
            Bm = [r['bound_S'] for r in rows if r['l'] == l and r['a'] == -1][0]
            N0 = rows[0]['N0']; E.append(max(2 * B1 - N0 - p, 2 * Bm - N0 - p))
        ax[1].plot(Ls, E, mk + '-', label=f'p={p}: built r (CHECKED)')
        ax[1].axhline(abs(rows[0]['N_aff'] - p), ls=':', color='k' if p == 101 else 'gray', label=f'p={p}: true |N-p|')
        ax[1].axhline(2 * math.sqrt(p), ls='--', color='C3' if p == 101 else 'C1', label=f'p={p}: Hasse 2 sqrt(p)')
    ax[1].set_xlabel('multiplicity l'); ax[1].set_ylabel('bound on |N_aff - p|')
    ax[1].set_title('Explicit auxiliary polynomials, y^2=x^3+2x+1'); ax[1].legend(fontsize=7)
    rows = [r for r in out['part1'] if r['m'] == 3 and r['a'] == 1 and r['p'] == 101]
    ax[2].plot([r['l'] for r in rows], [r['MAIN'] for r in rows], 'k-', label='MAIN = deg r (size)')
    ax[2].plot([r['l'] for r in rows], [r['ASK'] for r in rows], 'C0o-', label='ASK = sum_{x in S} ord_x r (divisibility)')
    ax[2].plot([r['l'] for r in rows], [r['l'] * r['S'] for r in rows], 'C2^--', label='forced: l |S|')
    ax[2].set_xlabel('multiplicity l'); ax[2].set_title('p=101, a=1: forced divisibility <= size (r != 0)'); ax[2].legend(fontsize=8)
    plt.tight_layout(); plt.savefig('toy_fq.png', dpi=130)

if __name__ == '__main__':
    main()
