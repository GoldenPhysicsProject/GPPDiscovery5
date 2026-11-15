# Round 11c: C-KG checked against O7, and C-Split checked against the explicit formula

Heart (Grok / GPPDiscovery5). Context: Turbo adopted AKG (PR #11, 6890778106). O7 and C-Split are from Round 11b (d1a4510).
Script: `o7_lowband.py`, output in `o7_lowband_results.json`.
- Zeros are a comparison target only. They come from Daniel's Platt table `zeros_14.dat` (GPP Drive "Zeros"), parsed with the Bober block format. The first 3000 zeros (γ ≤ 3533.33) agree with mpmath.zetazero to ≤ 2e-31 at n = 1, 2, 10, 100.
- Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN. AKG, H_ψ, RH: **OPEN**.

## 1. The exact split identity (CHECKED: classical explicit formula, same normalisation as K = pole − prime + arch)
Fix t₀ = 6.289835988836902779… (the zero of w(t) = Re ψ(1/4 + it/2) − log π). For any admissible even test function h of positive type, write arch = arch₋ + arch₊ for the parts of the arch integral over |t| < t₀ and |t| > t₀, and Z(h) = Σ_ρ ĥ(γ_ρ) for the zero side. Then
  **B(h) := arch₋ − prime + pole = Z(h) − arch₊(h).**
- B is C-Split's bracket. Since w > 0 on |t| > t₀ and ĥ ≥ 0, arch₊(h) > 0.
- **Zero-free band (CHECKED, classical verified computation, not RH):** ζ has no zeros with |Im ρ| < γ₁ = 14.1347…, and t₀ < γ₁.
- So Z(h) carries only the part of ĥ's mass that lies at |t| ≥ γ₁, while arch₊ carries everything above t₀.

## 2. Gaussian family: NUMERICAL, 40 digits, Platt zeros as comparison only
Take ĥ(t) = e^{−t²/(2s²)}.

| s | pole | prime | arch₋ | arch₊ | W = pole − prime + arch | zero side Z | B = arch₋ − prime + pole |
|---|---|---|---|---|---|---|---|
| 1 | 2.26630 | 1.06995 | −1.19635 | 3.00e-12 | 1.3e-31 | 8.3e-44 | **−3.00e-12** |
| 1.5 | 2.11426 | 0.65008 | −1.46417 | 8.14e-7 | 1.04538e-19 | 1.04538e-19 | **−8.14e-7** |
| 2 | 2.06349 | 0.40878 | −1.65481 | 1.0712e-4 | 2.8513e-11 | 2.8513e-11 | **−1.0712e-4** |
| 2.5 | 2.04040 | 0.24911 | −1.79266 | 1.3662e-3 | 2.28875e-7 | 2.28875e-7 | **−1.3660e-3** |
| 3 | 2.02797 | 0.14183 | −1.89265 | 6.5372e-3 | 3.02407e-5 | 3.02407e-5 | **−6.5069e-3** |
| 4 | 2.01569 | 0.03363 | −2.01980 | 4.1636e-2 | 3.88829e-3 | 3.88829e-3 | **−3.7747e-2** |

- The terms are O(1), yet they cancel down to the zero side: |W − Z| ≤ 3e-35 for s ≥ 1.5, and 1.3e-31 at s = 1 (prime and quadrature floor). This is the explicit formula confirmed to 30+ digits on low-frequency test functions.
- **B < 0 in every row, and B/arch₊ → −1** as the mass moves below γ₁ (s = 2: B/arch₊ = −0.99999973).

## 3. Consequences
- **O7 sharpened (CHECKED from §1 and the zero-free band).** At low frequency the prime side does not supply positivity with a margin. It cancels: −prime + pole + arch₋ = −arch₊ + Z. On test functions whose spectral mass lies below γ₁, the whole Weil form is the exponentially small zero side, and every O(1) term is consumed by the cancellation. O7's wording "low-frequency positivity from the prime side" should read "**low-frequency nullity**: primes + pole cancel arch₋ exactly, leaving −arch₊ + Z".
- **C-Split as stated: KILLED (CHECKED via §1; NUMERICAL via §2).**
  - **Margin requirement.** The bracket B = Z − arch₊ is negative on every test function whose mass in the band (t₀, γ₁) exceeds its zero-side weight. On the strictly low band (ĥ supported in |t| < t₀) it is exactly null where defined. "λ_min(B) bounded away from 0 in L" therefore fails.
  - **Convergence on the low band.** On that band the prime sum Σ Λ(n) n^{−1/2} h(log n) does not converge absolutely, because band-limited h cannot decay like e^{−|x|/2} (the same mechanism as J5).
  - **PW_L ∩ ker ℓ_c prediction (CONJECTURE, for the C-Split test in progress).** λ_min(B; N) < 0 as soon as the grid can resolve the band (t₀, γ₁), roughly 2π/L < γ₁ − t₀ = 7.84, i.e. from L ≈ 0.8 upward. On ker ℓ_c the pole part (−2Ψ²S² ≤ 0) only pushes B further down.
- **New obstruction O8 (CHECKED from §1):** any AKG norm satisfies ‖Φ_g‖² = Z(g), so it must vanish to the order of Z on test functions concentrated in the zero-free band |t| < 14.13. A construction from primes and arch has to reproduce the zero-free band implicitly. Any split that keeps positivity term by term must be placed at or above γ₁, not at t₀. Placing it at γ₁ uses only the classical verified band (no RH input), but it is still zero information, so it is diagnostic per Turbo's guardrail.
- **C-KG against O7: consistent, no conflict (CHECKED).** T_KG c = (ψ̂(γ)F_c(γ))_γ has no component below γ₁, which is exactly the low-band nullity of §3. T_KG is built from the γ's, so it remains diagnostic only (Turbo R11.5 guardrail), not a construction.
- **Cesàro candidate:** none exists in this lane. Any Cesàro/Fejér-weighted zero sum Σ (1 − |γ|/T)₊ ψ̂(γ)²|F_c(γ)|² is also zero-built, hence diagnostic only, and O8 applies to it unchanged.

## 4. What a viable all-L candidate must now satisfy
- (i) O2: it breaks translation invariance.
- (ii) O7/O8: it is null to exponential order on the zero-free band, so the low-band cancellation primes + pole + arch₋ = −arch₊ + Z must be built in exactly, not bounded.
- (iii) Its positive part lives at |t| ≥ γ₁ and equals the zero side there.
- Concretely: a prime-side Gram must reproduce −arch₊ on the band (t₀, γ₁) **exactly**. That is a sharp, testable requirement for AKG attack §8.1 (Connes semilocal space).

Credit: Daniel Toupin / Golden Physics Project (zero tables, argument). This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
