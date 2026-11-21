"""Round 32 code VALIDATION ONLY (zeros are test data, never inputs to any claim): under RH + simple zeros the
Moebius floor is m1~(L) = sum_rho e^{(rho-1/2)L}/(zeta'(rho)(rho-1/2)^2) + (trivial zeros, O(e^{-5L/2})).
Compare with the integer-side m1~ from r32_series.npz at a few L, first 300 zeros."""
import json, numpy as np, mpmath
mpmath.mp.dps = 20
z = np.load('r32_series.npz'); L = z['ZETA_Lg']; m = z['ZETA_m1t']
rh = []
for k in range(1, 301):
    r = mpmath.zetazero(k); g = float(r.imag); d = complex(mpmath.zeta(r, derivative=1)); rh.append((g, d))
out = []
for l in (6.0, 9.0, 12.0, 15.0, 18.0, 20.5):
    i = int(np.argmin(abs(L - l)))
    zs = sum(2*(np.exp(1j*g*L[i])/(d*(1j*g)**2)).real for g, d in rh)
    out.append(dict(L=float(L[i]), integer_side=float(m[i]), zero_side_300=float(zs), diff=float(m[i] - zs)))
    print(out[-1], flush=True)
json.dump(dict(note='validation only; zeros are test data', rows=out, tail_bound_hint='sum_{gamma>T} 1/(|zeta\'|gamma^2) ~ 1e-3'),
          open('r32_validate_zeros.json', 'w'), indent=1)
