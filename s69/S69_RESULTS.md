# Grok S69 H69a: Mori-route T at L=7/5 (support 2.8), m=1/2, even parity, Omega sweep (2026-10-21)

All values NUMERICAL (mpmath Cholesky bisection on the certified-assembly midpoint dump); none is a certificate.

| Omega | N | lambda_min(T_even) | log10 abs | centroid deg | Cholesky first fail (l0=1e-250) |
|---|---|---|---|---|---|
| 488 | 1142 | -4.66e-49 | -48.33 | ~10 | (S66/S68) |
| 600 | 1404 | -9.39e-62 | -61.03 | 9.8 | pivot 113 |
| 700 | 1638 | -6.84e-70 | -69.17 | 7.85 | pivot 554 |

Negative direction stays low-degree (52% weight on first 5 modes, 99.97% on first 10; top modes 10, 8, 12, 0).
Decay per 100 Omega slows: 11.3 decades (488 to 600), then 8.1 (600 to 700). Consistent with Astra's exp(b_L) resource barrier.
Nearest-zero eigenvalue estimate is ~3e-77, so the negative part must fall about 7 more decades before T could be PSD.
Next: H69b, Omega 800 (N 1872, delta_log2 3300 > 3256.79), launched 11:15Z (run/chain_s69b.log).
Rules out: nothing yet (no positivity claim). RH OPEN.

## H69b result, Omega 800 (N 1872, delta_log2 3300), even parity: Mori-route T_even CERTIFIED PSD (12:2xZ)

- Assembly certified: sha256 089907f4…f5e7, consts_ok true, radT_rowsum 1.45e-145, alpha+beta*b^2 = 1.29e-238.
- Approx Cholesky of T - 1e-250 I now SUCCEEDS (minpiv 7.30e-9); at O700 it failed at pivot 554.
- Nearest-zero eigenvalue estimate 8.0702e-79. Bisection (h68b) never failed down to l0 = -1e-80 (consistent with T > 0).
- s36 cholres (640 bits) at l0 = 0.999 est: CERTIFIED lambda_min(Atil) >= 8.0620e-79 (||E|| <= 1.25e-190, dump radius 9.36e-156).
- s37 exact fixed-point LDL^T (sha256 OK): l0 = 8.0621e-79 CERTIFIED lambda_min(sym Atil) >= 8.0621e-79 (||E||_F <= 2.56e-121);
  CONTROL l0 = 8.1106e-79 (1.005 est) FAILS at pivot 935. So the bound is tight to 0.6%.
- Hence lambda_min(T_even) >= 8.06e-79 - 1.45e-145 > 0, CERTIFIED (interval arithmetic), at L = 7/5 (support 2.8), m = 1/2.
- Sign change between Omega 700 (-6.84e-70, numerical) and 800 (+8.07e-79 certified): the low-degree negative mode is gone.
- Scope: EVEN parity only. Odd parity at Omega 800 launched 12:28Z (run/chain_s69c.log). Support 2.8 itself was already
  certified in S63 by another route, so this is an independent cross-route confirmation (even half), not a new window and not a record.
- Fakes: the T matrix carries zeta's prime sum, so it says nothing for DH or F_theta (no claim there); backward heat flow H_-t not addressed.
RH OPEN.

## H69c result, Omega 800 odd: NOT PSD (numerical)
- Cholesky of T_odd - 1e-250 I fails at pivot 692; bisection lambda_min ~ -1.51e-66 (NUMERICAL), low-degree vector (centroid deg 9.06). See run/negvec_L75_m12_O800_p1.log.

## H69d result, Omega 900 (N 2106, delta_log2 3720 > -log2 d(c) = 3661.08), odd parity: Mori-route T_odd CERTIFIED PSD (14:15Z)
- Assembly certified: sha256 c381f4c3…0ae9, consts_ok true, radT_rowsum 2.80e-145, alpha+beta*b^2 = 3.83e-269.
- Approx Cholesky of T - 1e-250 I succeeds (minpiv 3.28e-8). Nearest-zero eigenvalue estimate 2.8325e-74.
- s36 cholres (640 bits) at 0.999 est: lambda_min(Atil) >= 2.8300e-74.
- s37 exact fixed-point LDL^T (sha256 OK): l0 = 2.8296e-74 CERTIFIED lambda_min(sym Atil) >= 2.8296e-74 (||E||_F <= 2.49e-126);
  CONTROL l0 = 2.8466e-74 (1.005 est) FAILS at pivot 1051. Tight to 0.6%.
- Hence lambda_min(T_odd) >= 2.83e-74 - 2.8e-145 > 0, CERTIFIED, at L = 7/5 (support 2.8), m = 1/2, Omega 900.
- Sign change for odd between Omega 800 (-1.51e-66 numerical) and 900 (+2.83e-74 certified).
- Status: even certified at Omega 800, odd certified at Omega 900. Whether mixed Omega per parity is admissible in the
  Mori-route certificate is NOT yet checked, so to be safe even Omega 900 was launched 14:16Z (run/chain_s69e.log) for a matched pair.
- Independent cross-route check of S63 support 2.8 (both parities pending matched Omega), not a new window, not a record.
- Fakes: T carries zeta's Euler product, so no claim for DH or F_theta; backward heat flow H_-t not addressed. RH OPEN.
