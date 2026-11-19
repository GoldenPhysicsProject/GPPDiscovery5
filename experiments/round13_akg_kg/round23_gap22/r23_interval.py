"""Round 23, route (b): computer-assisted proof (ball arithmetic, python-flint arb) that
    D(L) > 0  for all L in [log 2, log X1],   i.e.  rho_X(0) > 0 for every real X in [2, X1].
D(L) = A(L) - S(L),  A(L) = 4(e^{L/2}-2+e^{-L/2}) + c1 L + C2 - T(L),  S(L) = alpha_j L - beta_j on [log q_j, log q_{j+1}],
alpha_j = sum_{n<=q_j} Lam(n) n^{-1/2}, beta_j = sum_{n<=q_j} Lam(n) n^{-1/2} log n  (q_j = consecutive prime powers).
Lemma (convexity, PROVED by hand): A''(L) = e^{L/2}+e^{-L/2} - e^{-L/2}/(1-e^{-2L}) > 0 for L >= log 2, since
(e^L+1)(1-e^{-2L}) is increasing and equals 9/4 > 1 at L = log 2.  S is affine on each gap, so D is convex on each gap.
A convex function on [a,b] is >= max of its two endpoint tangents; the minimum of that max is the certified lower bound.
At L = log q every exponential is algebraic in q: e^{L/2} = sqrt q, T(log q) = sum_k q^{-(2k+1/2)}/(2k+1/2)^2 (+ rigorous tail).
Inputs: only integers (the sieve) and the constants gamma, pi, log 2, Catalan G (c1 = (psi(1/4)-log pi)/2 = (-gamma-pi/2-3log2-log pi)/2,
C2 = psi'(1/4)/4 = (pi^2+8G)/4). No zeta zeros are used."""
import sys, math, json, time, numpy as np
from flint import arb, ctx
ctx.prec = 160
X1 = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
c1 = (-arb.const_euler() - arb.pi()/2 - 3*arb(2).log() - arb.pi().log())/2
C2 = (arb.pi()**2 + 8*arb.const_catalan())/4
K = 40
def AA(q=None, L=None):
    """returns (A, A') at L = log q (q integer) or at a general arb L."""
    if q is not None:
        r = arb(q).sqrt(); L = arb(q).log(); x = 1/arb(q)     # x = e^{-L}
    else:
        r = (L/2).exp(); x = (-L).exp()
    T = arb(0); T1 = arb(0); e = 1/r                           # e = e^{-L/2} = e^{-(2k+1/2)L} at k=0
    x2 = x*x
    for k in range(K):
        b = arb(2*k) + arb(0.5); T += e/(b*b); T1 += e/b; e *= x2
    b = arb(2*K) + arb(0.5); tail = e/(b*b*(1 - x2)); tail1 = e/(b*(1 - x2))
    T += arb(0, tail.upper()); T1 += arb(0, tail1.upper())   # adds [-tail,tail] (tail>0, enclosure is wider, still valid)
    A = 4*(r - 2 + 1/r) + c1*L + C2 - T
    A1 = 2*(r - 1/r) + c1 + T1
    return A, A1, L
def prime_powers(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    ps = np.nonzero(s)[0]
    out = [(int(p), int(p)) for p in ps]
    for p in ps[ps <= int(N**0.5)]:
        q = int(p)*int(p)
        while q <= N: out.append((q, int(p))); q *= int(p)
    out.sort(); return out
if __name__ == '__main__':
    t0 = time.time(); pp = prime_powers(X1); print('prime powers', len(pp), time.time() - t0, flush=True)
    alpha = arb(0); beta = arb(0)
    worst = None; worst_gap = None; nsub = 0; checked = 0
    pts = pp + [(X1, None)] if pp[-1][0] != X1 else pp
    # gap list: [log q_j, log q_{j+1}] with S = alpha_j L - beta_j, alpha_j includes q_j
    q0, p0 = pts[0]
    lp = arb(p0).log(); r = arb(q0).sqrt(); alpha += lp/r; beta += lp*arb(q0).log()/r
    Aa, A1a, La = AA(q=q0)
    Da = Aa - (alpha*La - beta); D1a = A1a - alpha
    rec = []; nextck = 10
    for q, p in pts[1:]:
        Ab, A1b, Lb = AA(q=q)
        Db = Ab - (alpha*Lb - beta); D1b = A1b - alpha
        # certified lower bound of the convex D on [La, Lb]
        if D1a.lower() >= 0: lb = Da
        elif D1b.upper() <= 0: lb = Db
        else:
            # tangents: Da + D1a (L-La), Db + D1b (L-Lb) ; min over L of max(...) >= value at crossing (enclosed)
            # lower bound valid for any slopes enclosure: use min(Da,Db) - (positive part) ... compute crossing in arb
            den = D1a - D1b
            if den.upper() < 0:   # convexity: D'(b-) > D'(a+), so den < 0
                Ls = (Db - Da + D1a*La - D1b*Lb)/den
                lb = Da + D1a*(Ls - La)
            else:
                lb = arb(-1)
        if lb.lower() <= 0 or float((Lb - La).upper()) > 0.01:
            # subdivide the gap (wide gaps at small q, or if the two-tangent bound is not positive): bisect with general-L evaluation, 64 pieces
            nsub += 1; lbs = []
            M = 64
            for i in range(M):
                a_ = La + (Lb - La)*i/M; b_ = La + (Lb - La)*(i + 1)/M
                A_a, A1_a, _ = AA(L=a_); A_b, A1_b, _ = AA(L=b_)
                da = A_a - (alpha*a_ - beta); db = A_b - (alpha*b_ - beta)
                sa = A1_a - alpha; sb = A1_b - alpha
                if sa.lower() >= 0: l = da
                elif sb.upper() <= 0: l = db
                else:
                    Ls = (db - da + sa*a_ - sb*b_)/(sa - sb); l = da + sa*(Ls - a_)
                lbs.append(l.lower())
            lbv = min(lbs)
        else:
            lbv = lb.lower()
        lbv = float(lbv)
        if worst is None or lbv < worst: worst = lbv; worst_gap = (None, q)
        if lbv <= 0:
            print('FAIL at gap ending', q, lbv, flush=True); break
        checked += 1
        if p is not None:
            lp = arb(p).log(); r = arb(q).sqrt(); alpha += lp/r; beta += lp*Lb/r
            Da = Ab - (alpha*Lb - beta); D1a = A1b - alpha   # value continuous; right slope drops by Lam/sqrt q
        La = Lb
        if q >= nextck or q == X1:
            nextck *= 10
            rec.append(dict(X=q, worst_lower_bound_so_far=worst, time=time.time() - t0)); print(rec[-1], flush=True)
    res = dict(X1=X1, gaps_checked=checked, subdivided_gaps=nsub, certified_min_lower_bound=worst,
               worst_gap_right_end=worst_gap[1], all_positive=bool(worst > 0), prec_bits=ctx.prec, checkpoints=rec,
               time_s=time.time() - t0)
    print(res, flush=True)
    json.dump(res, open(f'r23_interval_X{X1:.0e}.json'.replace('+', ''), 'w'), indent=1)
