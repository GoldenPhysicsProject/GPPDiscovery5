"""Grok S77: Arb (python-flint) check of Astra's zero-damping energy identity
and of the terminal-state / future-tail lemma, using exact per-cell antiderivatives.

Objects (t = log x, Astra 7b064a7):
  J0(t)  = int_0^t (e^{-u/2} psi(e^u) - e^{u/2}) du
         = 2 S1 - 2 psi(x)/sqrt(x) - 2 sqrt(x) + 2,           S1 = sum_{n<=x} Lambda(n)/sqrt(n)
  j(t)   = J_{1/2}(t) = 2 int_0^t q - e^{-t/2} J0(t),        q(u) = e^{-u} psi(e^u) - 1
  c      = -2(1+gamma)
  r(t)   = e^{t/2}(j-c) = 2 sqrt(x) (S2 - log x + 2 + gamma) - 2 S1 - 2,  S2 = sum Lambda(n)/n
  w(t)   = J0 - c e^{-t/2},   x_state = r - w.
Checks:  W(T) = int_0^T r^2 = int_0^T w^2 + x_state(T)^2   (Astra identity)
         x_state(T) -> -2 Jinf, W(T)/T -> Jinf^2 (RH prediction), Jinf = 2 - 2 zeta'/zeta(1/2).
All cell integrals are closed-form in Arb balls (prec 192), so the finite-X values are CERTIFIED
enclosures for the stated X (prime sums exact up to ball rounding). The asymptotic comparisons are
NUMERICAL evidence only. RH OPEN.
"""
import json, sys, math
from flint import arb, ctx
import mpmath as mp
ctx.prec = 192
X = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
# prime powers up to X
sieve = bytearray([1])*(X+1); sieve[0]=sieve[1]=0
for i in range(2, int(X**0.5)+1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(range(i*i, X+1, i)))
pp = []
for p in range(2, X+1):
    if sieve[p]:
        lp = arb(p).log(); n = p
        while n <= X:
            pp.append((n, lp)); n *= p
pp.sort(key=lambda z: z[0])
gam = arb.const_euler(); c = -2*(1+gam)

def E(k, t):  # e^{k t}
    return (k*t).exp()
def A0(k, t):  # antiderivative of e^{kt}
    return E(k, t)/k
def A1(k, t):
    return E(k, t)*(t/k - 1/(k*k))
def A2(k, t):
    return E(k, t)*(t*t/k - 2*t/(k*k) + 2/(k**3))

def cell(P, S1, S2, a, b):
    """int_a^b r^2 and int_a^b w^2 with psi=P, S1, S2 constant on [a,b)."""
    al = 2*(S2 + 2 + gam); be = arb(-2); de = -2*S1 - 2
    D = lambda F, k: F(arb(k), b) - F(arb(k), a)
    r2 = (al*al*D(A0, 1) + be*be*D(A2, 1) + de*de*(b-a) + 2*al*be*D(A1, 1)
          + 2*al*de*D(A0, 0.5) + 2*be*de*D(A1, 0.5))
    Ac = 2*S1 + 2; Bc = -2*P - c
    w2 = (Ac*Ac*(b-a) + Bc*Bc*D(A0, -1) + 4*D(A0, 1) + 2*Ac*Bc*D(A0, -0.5)
          - 4*Ac*D(A0, 0.5) - 4*Bc*(b-a))
    return r2, w2

def rw(P, S1, S2, t):
    y = (t/2).exp()
    r = 2*y*(S2 - t + 2 + gam) - 2*S1 - 2
    w = (2*S1 - 2*P/y - 2*y + 2) - c/y
    return r, w

mp.mp.dps = 40
Jinf = 2 - 2*mp.zeta(0.5, derivative=1)/mp.zeta(0.5)
checkpoints = sorted(set([10**k for k in range(2, 12) if 10**k <= X] + [X]))
out = {"X": X, "Jinf_numerical": mp.nstr(Jinf, 20), "rows": []}
P = arb(0); S1 = arb(0); S2 = arb(0); W = arb(0); Ec = arb(0); a = arb(0)
idx = 0; ci = 0
events = [n for n, _ in pp]
def report(T_int):
    global out
    T = arb(T_int).log()
    # integrate from a to T on current cell
    r2, w2 = cell(P, S1, S2, a, T)
    Wt = W + r2; Et = Ec + w2
    r, w = rw(P, S1, S2, T); xs = r - w
    resid = Wt - Et - xs*xs
    Tf = float(T.mid())
    row = {"X": T_int, "T": Tf, "W": Wt.str(12, radius=True), "Ec": Et.str(12, radius=True),
           "xT": xs.str(12, radius=True), "xT2": (xs*xs).str(12, radius=True),
           "identity_residual": resid.str(5, radius=True),
           "W_over_T": float(Wt.mid())/Tf,
           "xT_plus_2Jinf": float(xs.mid()) + 2*float(Jinf),
           "W_minus_Jinf2_T": float(Wt.mid()) - float(Jinf)**2*Tf}
    out["rows"].append(row); print(row, flush=True)
for i, (n, lp) in enumerate(pp):
    while ci < len(checkpoints) and checkpoints[ci] < n:
        report(checkpoints[ci]); ci += 1
    b = arb(n).log()
    r2, w2 = cell(P, S1, S2, a, b)
    W += r2; Ec += w2
    P += lp; S1 += lp/arb(n).sqrt(); S2 += lp/arb(n); a = b
while ci < len(checkpoints):
    report(checkpoints[ci]); ci += 1
json.dump(out, open(f"/workspace/GPPDiscovery5/s77astra/zero_damping_X{X}.json", "w"), indent=1)
