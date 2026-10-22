
## 0. Team asks since S27 (PR #11 6204288080, 6204310175 — Codex)
- **Float auxiliaries in the L = 1 certificate (6204288080).** Accepted and fixed. New `s28_tails_arb.py` redoes the whole final budget in Arb: Codex's exact ratio bound for the Legendre-block tail (a_d² = 2L(2d+1)(LT)^{2d}/((2d+1)!!)², r_d = (LT)⁴/[(2d+1)(2d+3)²(2d+5)] decreasing, tail ≤ a_d²/(1 − r_d)); the pole tail with e^{L²/4}; the pole head |p|² in Arb; a rational Bernstein ρ = 347/100 (as Codex proposed) for the quadrature term; and an inflated ‖A − Ã‖. The L = 1 certificate survives with tiny loss: **even ≥ 4.1248e−30, odd ≥ 1.05499e−26** (S27 floats: 4.131e−30 / 1.050e−26). I also fixed an exp() underflow in the float tail code (it returned 0 where the true value is below 1e−300; it now returns the upper bound 1e−300). The new fast pipeline (§2) also reproduces the S27 numbers: λ-estimate 4.142300e−30 vs the S27 Cholesky fail-shift 4.142302e−30. This is the same author, so it does not replace Turbo/Sol's independent rerun.
- **Codex's PNT barrier for prime-separable lifts (6204310175).** Agreed. It is the lower half of H27: κ_p ≥ c_p/2 against 1 + cos θ, so a separable A^lift ≥ (2 − o(1))e^L. §1 gives the matching statement for **all** lifts, cross-prime included, and shows how much cross-prime coupling can buy.

## 1. PROVED (elementary): every invisible lift is bounded below by λ_max of the prime-shift operator
Let P_L be the operator on L²[−L, L] given by
  (P_L f)(x) = Σ_{log n < 2L} Λ(n) n^{−1/2} [f(x + log n) + f(x − log n)]   (values outside [−L, L] dropped),
so ⟨f, P_L f⟩ = (1/2π)∫ comb(t)|F(t)|² dt with comb(t) = Σ 2Λ(n)n^{−1/2} cos(t log n). This is the prime part of the Weil form on the window.
**Proposition.** Let m be any almost-periodic invisible lift (spectrum in |u| ≥ 2L). This covers S27's per-prime lifts and S28's torus lifts. Then
  sup_t (comb − m)(t) ≥ λ_max(P_L).
*Proof.* Take f in the window with ‖f‖ = 1. The invisibility lemma gives (1/2π)∫ m|F|² = 0, so ⟨f, P_L f⟩ = (1/2π)∫(comb − m)|F|² ≤ sup(comb − m) · (1/2π)∫|F|² = sup(comb − m). For almost-periodic m, the sup over t ≥ T equals the sup over all of ℝ, so the bound also holds for the only part the one-stroke envelope uses. ∎
Consequences:
- Any lifted pointwise-envelope one-stroke certificate needs β = log(T/2π) − 1/T − A > 0, so T > T₁* := 2π exp(λ_max(P_L)).
- With f constant, ⟨f, P_L f⟩ = Σ 2Λ(n)n^{−1/2}(1 − log n/2L), which is (4 − o(1))e^L/L by the PNT. So T₁* is doubly exponential for the **whole lifted architecture**, cross-prime coupling included. This extends Codex's separable barrier. It is a barrier for this method only, not for Weil positivity.
- HYPOTHESIS (strong duality): inf over lifts of sup(comb − m) = λ_max(P_L). The LP dual is a probability measure on the Bohr compactification (the prime torus) with spectrum in (−2L, 2L). The densities |F|²/‖F‖² are such measures, and Krein's extension theorem for positive-definite functions on (−2L, 2L) should give equality. Not written out.
NUMERICAL λ_max(P_L) (orthonormal Legendre Galerkin with exact shifted overlaps, N up to 320; these are Rayleigh lower bounds that increase with N; `s28_pmax.py`, `pmax.json`, `pmax_conv.log`):

| L | raw A_L | S27 per-prime A^lift | S28 torus (D=3) | floor λ_max(P_L) ≥ | λ_max/e^L | T₁* = 2πe^{λ_max} |
|---|---|---|---|---|---|---|
| 0.5 | 0.980 | — | — | 0.4901 | 0.30 | 10.3 |
| 0.8 | 2.942 | 1.535 | — | 1.2191 | 0.55 | 21 |
| 1.0 | 5.852 | 3.048 | 2.933–2.954 (not converged) | 1.9385 | 0.71 | 44 |
| 1.19 | 7.075 | 3.919 | 3.602 (LP 3.6019, max-search 3.6023) | 2.6656 | 0.81 | 90 |
| 1.6 | 14.32 | 7.82 | — | 4.8953 | 0.97 | 840 |
| 2.0 | 24.38 | 13.61 | — | 7.9726 | 1.08 | 1.8e4 |
| 2.5 | 42.55 | — | — | 13.635 | 1.12 | 5.2e6 |
| 3.0 | 75.19 | — | — | 22.496 | 1.12 | 3.7e10 |

**H28 (hypothesis, replaces H27's upper half):** the optimal lift constant is λ_max(P_L) ≈ (1.1 ± 0.05)e^L on 2 ≤ L ≤ 3 (NUMERICAL), roughly half of the separable 2e^L. Asymptotically it lies between (4/L)e^L (constant test function) and ~2e^L. So cross-prime lifts can cut log T₁ roughly in half, but T₁ stays doubly exponential.

## 1b. NUMERICAL: the torus LP (H28a)
`s28_torus_lp.py` is a cutting-plane LP on T^n (n = primes below e^{2L}). The variables are the coefficients of the invisible characters k·θ with |Σ k_p log p| ≥ 2L: the axis harmonics up to 6× the first invisible order, plus cross characters with |k|₁ ≤ D. The objective is the torus max of B − m. By Kronecker, the sup over t equals the sup over the torus.
- Two primes {3, 5} at L = 1: separable (kmax 12) 1.3889; exact separable c₃/2 + c₅/2 = 1.3540; torus D = 4: 1.3353; D = 6: **1.2306**. CERTIFIED for this toy (float64 branch-and-bound, second-order Taylor cells, explicit rounding margin; `s28_torus_bb.py`): sup ≤ 1.2336. Cross-prime coupling beats Codex's separable floor c_p/2, as §1 allows.
- Full L = 1 (primes 2, 3, 5, 7; 97 characters at D = 3): LP 2.933, max-search ≥ 2.954, not converged, vs per-prime 3.048 and floor 1.94.
- Full L = 1.19 (114 characters, D = 3): LP 3.6019, max-search 3.6023 (converged to 3e−4), vs per-prime 3.919 and floor 2.666. That lowers the lifted T₁ from 317 to about 231.
- D = 5 runs (284–291 characters) were too slow in this session (round 0 alone took over 5 min with gap 0.4–0.7). The L = 2 pair-group run did not finish. Both were stopped (not dead ends, just not done).
- The gap to the floor stays large at D = 3 (3.60 vs 2.67 at L = 1.19). Closing it needs characters with large |k|₁, where the visible set is dense (many k with |k·log p| < 2L). Better route (HYPOTHESIS): build m directly from the Krein extension of A·δ₀ − Σ(c_n/2)δ_{±log n} beyond (−2L, 2L), using the top eigenvector of P_L, instead of a character LP.
