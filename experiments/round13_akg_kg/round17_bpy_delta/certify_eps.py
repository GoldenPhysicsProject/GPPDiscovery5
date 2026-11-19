"""Certified sign: Phi_EPS(0) < 0, i.e. the Epstein (x^2+5y^2) BPY-type law is NOT a positive measure.
Phi_E(0) = sum'_{(m,n)} (a^2 - 2a) e^{-a},  a = pi (m^2+5n^2)/sqrt5.  Terms with a > 2 (Q >= 2) are positive, so
Phi_E(0) <= S(Q0) + T(Q0) with S the exact finite sum over Q <= Q0 (arb ball arithmetic) and the tail bounded using
r_Q <= 2(2 sqrt(Q/5)+1) and a^2 e^{-a} decreasing for a > 2.  Output certify_eps.json"""
import json, math, os
from flint import arb, ctx
ctx.prec = 300
Q0 = 400
r = [0]*(Q0 + 1); M = int(math.isqrt(Q0)) + 1
for m in range(-M, M + 1):
    for n in range(-M, M + 1):
        Q = m*m + 5*n*n
        if 0 < Q <= Q0: r[Q] += 1
c = arb.pi()/arb(5).sqrt()
S = sum((r[Q]*((c*Q)**2 - 2*c*Q)*(-(c*Q)).exp() for Q in range(1, Q0 + 1) if r[Q]), arb(0))
# tail: sum_{Q>Q0} 2(2 sqrt(Q/5)+1) (cQ)^2 e^{-cQ} <= sum_{Q>Q0} 2(2Q+1)(cQ)^2 e^{-cQ} <= (geometric) bound
q = Q0 + 1; first = 2*(2*q + 1)*(c*q)**2*(-(c*q)).exp()
ratio = arb(1.01)**3*(-c).exp()      # consecutive-term ratio bound for Q >= 400: ((Q+1)/Q)^3 e^{-c} < 1.01^3 e^{-c}
T = first/(1 - ratio)
upper = S + T
out = dict(partial_sum=S.str(30, radius=True), tail_bound=T.str(5, radius=True), upper_bound=upper.str(30, radius=True),
           certified_negative=bool(upper < 0), label='CHECKED (ball arithmetic, 300 bits)')
print(out)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'certify_eps.json'), 'w'), indent=1)
