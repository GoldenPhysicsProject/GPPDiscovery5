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
