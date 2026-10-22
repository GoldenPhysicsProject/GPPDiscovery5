# Does Q_L^{P_X} ≥ 0 survive past the cutoff? (open question §4.1 of finite_windows_dmv.md)

Goddess of Prime Number Theory, 2026-10-22. **RH is OPEN.** Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
**Convention:** L is the half-width; atoms q^k with k·log q < 2L enter the window form.

**Setup.**
- P_X = the rational primes ≤ X = 23, plus synthetic g-primes above X, given as exact rationals.
- Assembly: `precision/cosritz_px.py` (s62 Arb assembly; env GPRIMES replaces the prime atoms).
- Validation: with GPRIMES = rational primes ≤ 29, the matrix agrees with the ζ assembly to 100 digits.
- Tests: even sector, cos basis.
  - CERTIFIED PD: `ldl_certify.py` at 576 bits.
  - CERTIFIED non-PSD: ball LDLᵀ with a negative pivot after positive ones.
  - Otherwise NUMERICAL: 576-bit float LDLᵀ against a 640-bit assembly.

## A. Just past the first synthetic atom (q₁ = 24.7, so it enters at L₀ = (log 24.7)/2 = 1.60337)
| L | 2L − log q₁ | P_X (primes ≤ 23 + 24.7), K = 160 | ζ at same L |
|---|---|---|---|
| 1.6035 | 0.00026 | **PD, CERTIFIED**, λ = 2.24e-121 | λ = 2.30e-121 |
| 1.6043 | 0.0018 | **not PSD, CERTIFIED** (pivot −2.1e-15) | PD, λ = 1.46e-121 |
| 1.6100 | 0.0132 | **not PSD, CERTIFIED** (pivot −1.2e-12) | PD, λ = 4.3e-123 |
| 1.6200 | 0.033 | not PSD, CERTIFIED (two synthetic atoms 24.7, 28.3; only 24.7 inside) | PD, λ = 9.5e-126 |

Also at L = 1.7105, K = 200: synthetic {24.7, 28.3} in place of 29 gives a pivot of −1.8e-3 (CERTIFIED non-PSD). Deleting 29 gives −8.3e-14 (CERTIFIED).

**Answer for this P_X (NUMERICAL, with certified witnesses):** positivity survives only an O(10^{-3}) penetration of the first non-rational
atom into the window. It then fails, for L₀ < L* < L₀ + 0.0008. The synthetic atom near the window edge costs only through g(±log q₁), which
vanishes at the edge. That is why there is a short grace interval.
**CONJECTURE:** for every P_X ≠ ℙ above X, Q_L^{P_X} fails for L slightly above (log q₁)/2, where q₁ is the first g-prime where P_X differs from ℙ.

## B. Pinning: how far can the prime 29 move at L = 1.7105 before positivity is lost? (K = 200, ζ λ_min = 6.4381899e-153)
| shift δ of 29 | +1e-20 | −1e-20 | +1e-30 | −1e-30 | +1e-40 | −1e-40 | ±1e-50, ±1e-60, ±1e-100 | ±1e-3 |
|---|---|---|---|---|---|---|---|---|
| result | not PD | not PD | not PD | not PD | PD (CERT), λ = 6.4389e-153 | PD (CERT), λ = 6.4375e-153 | PD (CERT), λ unchanged to 8 digits | not PSD (CERT) |

The response is linear: dλ_min/dδ ≈ +7e-117. The extrapolated tolerance on the −δ side is |δ| ≈ λ/|dλ/dδ| ≈ **1e-36**
(NUMERICAL). So **the window form at L = 1.7105 pins the position of the prime 29 to about 36 decimal places.** The pinning is much weaker
than the 1e-153 margin would suggest, because the soft eigenvector's autocorrelation g_v is tiny near u = log 29, i.e. close to the window edge 2L = 3.421.
Interpretation (for "where do the integers enter"): at finite L, positivity is not stable under moving a prime. ζ's prime configuration
sits on the edge of the PSD cone, at the 10^{-36} level for the outermost prime.
This is consistent with the finite-window theorem, which says *agreement* below the cutoff is invisible but says nothing past it.

## Status and gaps
- Section-level only (K = 160/200). There is no tail bound, so "PD" means PD of the K-section. The non-PSD results are genuine witnesses (Ritz vectors), hence CERTIFIED for the full form.
- One P_X family, one X. The ±1e-30 and ±1e-20 failures are NUMERICAL (576-bit float pivots, ball LDLᵀ unresolved).
- Runs: `precision/runs/px/*.jsonl`. Scripts: `run_px.sh`, `run_px2*.sh`, `run_px3.sh`.
