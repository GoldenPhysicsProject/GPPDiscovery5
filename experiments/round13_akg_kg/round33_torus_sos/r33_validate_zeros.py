"""Round 33 VALIDATION ONLY (zeros are test data, never proof inputs). Under RH + simple zeros the windowed mean squares
have slopes  kappa_mu  = sum_rho 1/(|rho-1/2|^2 |zeta'(rho)|^2)  (Moebius side, M(e^u) - 1/zeta(1/2))
             kappa_lam = sum_rho 1/|rho-1/2|^2 = (log xi)''(1/2) = R0   (psi side, (Lambda-1) vector)
First 300 zeros (both signs of gamma), plus a crude tail estimate for kappa_lam."""
import json, mpmath
mpmath.mp.dps = 20
km, kl = 0.0, 0.0; rows = []
for k in range(1, 301):
    r = mpmath.zetazero(k); g = float(r.imag); d = abs(complex(mpmath.zeta(r, derivative=1)))
    km += 2/(g*g*d*d); kl += 2/(g*g)
    if k in (10, 50, 100, 200, 300): rows.append(dict(K=k, gamma=g, kappa_mu=km, kappa_lam=kl)); print(rows[-1], flush=True)
json.dump(dict(note='validation only; zeros are test data', rows=rows,
               R0_integer_side=0.04620998623083794, kappa_lam_tail_hint='sum_{gamma>T} 2/gamma^2 ~ (1/pi) log(T/2pi)/T'),
          open('r33_validate_zeros.json', 'w'), indent=1)
