"""Independent check of ||C_p|| <= s on L^2(-L,L) (the input s of the Mori route), by a weighted Schur test
with a step weight found by power iteration (any w > 0 gives a valid bound; the iteration only finds a good w).
  ||C_p|| <= ess sup_x R_w(x),  R_w(x) = sum_n c_n [1_I(x+l_n) w(x+l_n) + 1_I(x-l_n) w(x-l_n)] / w(x),  c_n = Lambda(n)/sqrt(n).
(Schur test: 2|f(x)g(y)| <= |f(x)|^2 w(y)/w(x) + |g(y)|^2 w(x)/w(y).)
Rigorous evaluation: for x in cell j, x +- l lies in an interval covering cells [lo, hi] with lo, hi from Arb enclosures of
l/h, widened by one cell on each side; w(x+-l) <= max of w on those cells; 1_I <= 1; final sum in Arb with exact dyadic weights.
usage: cp_norm_check.py L Ncells [iters]      Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, math, time
import numpy as np
from flint import arb, ctx
ctx.prec = 128
Ls = sys.argv[1]; Nc = int(sys.argv[2]); iters = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
p_, _, q_ = Ls.partition('/'); L = arb(int(p_)) / (int(q_) if q_ else 1)
def pp(n):
    for p in range(2, n + 1):
        if n % p == 0:
            k = n
            while k % p == 0: k //= p
            return p if k == 1 else None
e2L = (2 * L).exp(); terms = []
n = 2
while bool(arb(n) < e2L):
    p = pp(n)
    if p: terms.append((n, arb(p).log() / arb(n).sqrt(), arb(n).log()))
    n += 1
assert bool(arb(n) > e2L)
h = 2 * L / Nc
wins = []   # (c_upper_arb, offset_lo, offset_hi) per signed shift
for (n, c, ell) in terms:
    for sg in (1, -1):
        r = sg * ell / h
        lo = int(math.floor(float(r.lower()))) - 1; hi = int(math.floor(float(r.upper()))) + 2
        wins.append((n, c, lo, hi))
def winmax(w, lo, hi):
    # out[j] = max(w[j+lo .. j+hi]) over indices inside [0,Nc), 0 if none
    pad = max(abs(lo), abs(hi)) + 2
    wp = np.concatenate([np.zeros(pad), w, np.zeros(pad)])
    out = np.full(Nc, 0.0)
    for o in range(lo, hi + 1):
        out = np.maximum(out, wp[pad + o: pad + o + Nc])
    return out
cf = [float(c.upper()) for (_, c, _, _) in wins]
w = np.ones(Nc)
t0 = time.time(); hist = []
for it in range(iters):
    Cw = sum(cf[k] * winmax(w, lo, hi) for k, (_, _, lo, hi) in enumerate(wins))
    R = Cw / w
    if it % 250 == 0: hist.append((it, float(R.max()))); print(it, R.max(), flush=True)
    wn = Cw + 2.0 * w; wn = 0.5 * (wn + wn[::-1]); w = wn / wn.max()
w = np.maximum(w, 1e-300)
# rigorous recheck in Arb: weights as exact dyadic numbers
Cw = sum(cf[k] * winmax(w, lo, hi) for k, (_, _, lo, hi) in enumerate(wins))
R = Cw / w; order = np.argsort(-R)
best = arb(0)
cands = order[:2000]
for j in cands:
    s = arb(0)
    for (n, c, lo, hi) in wins:
        idx = [j + o for o in range(lo, hi + 1) if 0 <= j + o < Nc]
        if idx: s += c * arb(float(max(w[i] for i in idx)))
    v = s / arb(float(w[j])); best = arb.max(best, arb(v.upper()))
# the float ranking can mis-order by rounding only; certify the rest by a float bound with margin
rest_max = float(R[order[2000:]].max()) * (1 + 1e-12) if Nc > 2000 else 0.0
bound = max(float(best.upper()), rest_max * (1 + 1e-10))
res = dict(L=Ls, Ncells=Nc, iters=iters, primes=[t[0] for t in terms], schur_bound=bound, arb_top=best.str(12),
           float_rest=rest_max, wmin_over_wmax=float(w.min() / w.max()), seconds=round(time.time() - t0, 1), hist=hist)
print(json.dumps(res)); json.dump(res, open('cp_norm_%s_N%d.json' % (Ls.replace('/', '_'), Nc), 'w'), indent=1)
