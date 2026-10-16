# S49 (Grok): Arb-certified Taylor coefficients of xi at 1/2 and the first real-diagonal
# Fisher margin  -m'_xi(0) = (a1^2 - 2 a0 a2)/a0^2,  f(u) = xi(1/2+sqrt u) = sum a_k u^k.
# Also the Turan form (gamma_1^2 - gamma_0 gamma_2)/gamma_0^2, gamma_k = k! a_k, and
# the Csordas-Norfolk-Varga form b_1^2 - (1/3) b_0 b_2 with b_m = xi^{(2m)}(1/2).
import json, flint
from flint import arb, acb, acb_series, ctx
ctx.prec = 400
N = 12  # series length in z
z = acb_series([acb(0.5), acb(1)], prec=N)          # s = 1/2 + z
half_s = acb_series([acb(0.25), acb(0.5)], prec=N)  # s/2
pi = acb(arb.pi())
xi = acb(0.5) * z * (z - 1) * (-(half_s * pi.log())).exp() * half_s.gamma() * acb_series.zeta(z)
c = xi.coeffs()
# odd coefficients must vanish (evenness)
odd = [c[k] for k in range(1, 9, 2)]
a = [c[2 * k] for k in range(5)]                 # a_k = xi^{(2k)}(1/2)/(2k)!
for k in range(5):
    assert a[k].imag.contains(0)
a = [x.real for x in a]
import math
b = [a[k] * math.factorial(2 * k) for k in range(5)]   # b_k = xi^{(2k)}(1/2)
g = [a[k] * math.factorial(k) for k in range(5)]       # Jensen/Turan gamma_k
M0 = (a[1] ** 2 - 2 * a[0] * a[2]) / a[0] ** 2
M0b = (3 * b[1] ** 2 - b[0] * b[2]) / (12 * b[0] ** 2)
T0 = (g[1] ** 2 - g[0] * g[2]) / g[0] ** 2
CNV = b[1] ** 2 - b[0] * b[2] / 3
# second Turan (k=2) and the next Hadamard sums for completeness
T1 = g[2] ** 2 - g[1] * g[3]
# power sums over u-roots: -log f(u) expansion: log f = log a0 + sum ... ; S2 = sum u_rho^{-2}
# m(u) = f'/f = -sum_k S_{k+1} u^k with S_j = sum_rho u_rho^{-j}; so -m'(0) = S_2 (here S_2 = sum 1/u_rho^2)
S1 = -a[1] / a[0]                         # sum 1/u_rho  (= -sum 1/gamma^2 under RH)
S2 = M0
print('odd coeffs contain 0:', all(x.real.contains(0) and x.imag.contains(0) for x in odd))
for k in range(5):
    print(f'xi^({2*k})(1/2) =', b[k].str(35, radius=True))
print("-m'_xi(0) =", M0.str(40, radius=True))
print("check (3b1^2-b0b2)/(12b0^2) =", M0b.str(40, radius=True))
print('Turan (g1^2-g0g2)/g0^2 =', T0.str(40, radius=True), ' equal to -m\'(0):', (T0 - M0).contains(0))
print('CNV b1^2 - b0 b2/3 =', CNV.str(30, radius=True), ' >0:', CNV > 0)
print('second Turan g2^2-g1g3 =', T1.str(25, radius=True), ' >0:', T1 > 0)
print('S1 = sum 1/u_rho = -a1/a0 =', S1.str(40, radius=True))
res = dict(prec_bits=ctx.prec, xi_derivs={f'xi^({2*k})(1/2)': b[k].str(40, radius=True) for k in range(5)},
           minus_mprime0=M0.str(45, radius=True), minus_mprime0_positive=bool(M0 > 0),
           turan_equal=bool((T0 - M0).contains(0)), CNV_b1sq_minus_b0b2_over3=CNV.str(30, radius=True),
           second_turan=T1.str(30, radius=True), S1_sum_inv_u=S1.str(40, radius=True), label='CERTIFIED (Arb ball arithmetic)')
json.dump(res, open('s49_turan_cert.json', 'w'), indent=1)
