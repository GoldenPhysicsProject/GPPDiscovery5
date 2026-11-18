# Round 13e — falsifier harness: Epstein x²+5y² and a DMV-type Beurling proxy

Same objects as Rounds 13–13d: the Klein–Gordon bulk field Φ_g(u) and its Cesàro (Hann) mean against the Weil value, the signed bound λ_max(prime − pole), and the two-space split W = arch + pole − prime (Galerkin on Dirichlet sines, M = max(24, 10L)). Widths L = 0.5, 1, 2, 3, 4. Zeros are used only as test data, never as inputs. Code: `harness.py`, `harness_checks.py`, `eps_field_zero_check.py`, `beur_recheck.py`, `beur_spliced.py`, `epstein_zeros*.py`. Results are in the `*.json` files next to them.

## 0. The counter-models (built from their own data)
- **EPS**: Z_Q(s) = Σ' (x²+5y²)^{-s} = ζ(s)L(s,χ₋₂₀) + L(s,χ₋₄)L(s,χ₅) (disc −20, h = 2). It has a functional equation with its own gamma factor: Λ_Q = (√20/2π)^s Γ(s) Z_Q = Λ_Q(1−s), with poles at 0 and 1 and the FE checked to 1.5e-18. It has **no Euler product**.
  - Explicit-formula data come from its own coefficients. The coefficients c_Q(n) of −Z_Q′/Z_Q come from the representation numbers r_Q(n), computed directly up to 3.4e6, through the Dirichlet log-derivative recursion.
  - The arch term comes from its own gamma factor: C_Q = 2 log(√20/2π) + 2ψ(1/2) and w_Q(t) = e^{−t/2}/(1−e^{−t}).
  - 76% of the atoms are not prime powers, and 28% of the c_Q(n) are negative. Σ_{n<X}|c_Q| n^{−1/2} / 2√X = 1.28, 2.75, 5.20 at X = 1e2, 1e4, 1e6 (for ζ: 0.84, 0.99, 1.00).
  - Zeros (test data): 39 on the line for t < 60. Off-line in Re s > 1/2, t ≤ 60.5, there are 7 zeros, a count confirmed by the argument principle (6.99999…). The lowest is **0.93297 + 15.66825i**.
  - The first grid pass missed it. Its gap 13.10 → 17.60 in the on-line list gave it away, and the supplementary search `epstein_zeros_lowband.py` found it.
- **BEUR** (closest computable DMV proxy): discrete Beurling primes. Π_B, with powers weighted 1/m, tracks Π_C, where dψ_C = 1_{x≥x₀}(1 − Σ_q x^{ρ_q−1})dx and ρ_q = 1/2 ± 0.1 ± 28i.
  - It has an Euler product and Λ_B ≥ 0, but **no functional equation and no gamma factor**.
  - Π_B − Π_C = O(1), so ζ_B/ζ_C is analytic and zero-free in Re s > 0. The zeros in Re s > 0 are exactly the quartet. The construction gives ψ_B − ψ_C = 2.2, −3.9, −0.5, +1.2 at 1e3 … 1e6.
  - Bug fixed: the first version snapped primes to the log-grid. That biased ψ_B by about −x·Δ/2 and faked e^{u/2} growth in the field. Primes are now placed at the interpolated crossing.
  - **BEUR″ (spliced)**: the real prime powers with p < 55, plus Beurling primes ≥ 55 (x₀ = 55). The Euler factors (1−p^{−s})^{−1} have no zeros or poles in Re s > 0, so BEUR″ still has exactly the quartet as its zeros in Re s > 0. Since e⁴ < 55, its atoms agree with ζ's on every window of width ≤ 4.
  - This is a proxy: DMV's systems have infinitely many zeros approaching σ = 1, while this one has a fixed quartet at σ = 0.6. The localization argument below only uses "Euler product, Λ ≥ 0, off-line zeros, no FE", which both have.

## 1. Validation
(Superseded files: the BEUR rows and the EPS zero-sum check in `harness_checks.json` predate the prime-placement fix and the complete off-line zero list. Use `beur_recheck.json`, `beur_spliced.json`, and `eps_field_zero_check.json` instead. The M-convergence and width-scan rows in `harness_checks.json` stay valid.)
- ζ: Hann mean of Φ_g² over [·, 13] / W = 1.0000 at all five widths (NUMERICAL, as in R13).
- EPS: the field built from EPS's own c_Q, arch, and poles equals the zero sum (39 on-line + 7 off-line quartets) to a relative 4.6e-5 / 5.8e-6 / 1.7e-6 at L = 2 / 3 / 4, for u > 9. This checks the EPS explicit-formula data end to end (NUMERICAL).
- BEUR″: field = quartet sum Σ_q G(ρ_q−1/2)e^{(ρ_q−1/2)u} to a relative ~1e-5 at L = 0.5 and 1. At L ≥ 2 it agrees to the size of the O(e^{−u/2}) discreteness residue (NUMERICAL).

## 2. Per-model results (L = 0.5 / 1 / 2 / 3 / 4)
| term | ζ | EPS | BEUR″ (spliced) |
|---|---|---|---|
| λ_min(arch) | −0.90 / −1.70 / −2.59 / −3.16 / −3.56 | **+1.44 / +0.06 / −1.29 / −2.06 / −2.57** (better than ζ) | identical to ζ (ζ's arch borrowed; BEUR has none of its own) |
| pole operator | 2(cc−ss) | identical | identical |
| λ_max(prime), absolute | 0 / 0.49 / 1.80 / 4.15 / 7.78 | 0 / 0 / 2.06 / 4.84 / 9.66 | identical to ζ (operator difference exactly 0.0) |
| λ_max(prime), signed | = absolute (Λ ≥ 0) | 0 / 0 / 2.06 / 4.84 / **8.89** | identical |
| **λ_max(prime − pole), signed bound** | 0.005 / 0.49 / 1.49 / 1.48 / 1.56 (**flat**) | 0.005 / 0.040 / 1.99 / 4.27 / **6.29** (**grows**) | identical to ζ |
| λ_min(W) | 3.6e-2 / 1.1e-6 / ≈0 (≤1e-13) | 2.38 / 1.92 / 1.12 / 3.2e-4 / **−1.56** | identical to ζ |
| field Cesàro/W at U = 13 | 1.0000 (all L) | **104.5** / 1.18 / 1.001 / 0.9993 / 0.9996 | quartet-dominated, not W: 0.059 / 0.0081 / 0.0021 / 0.36 / 1.67, and rising with u at L = 0.5, 1 |
| off-line share of field rms | 0 | L = 2: 0.4% (u = 3) → 16.7% (u = 12); L = 4: 0.03% → 0.98%; grows like e^{0.433u} | quartet only; grows like e^{0.1u} |

Notes:
- **EPS at L = 4.** λ_min(W_Q) = −1.56 / −1.60 / −1.62 at M = 40 / 56 / 72. There is a second negative eigenvalue, about −1.2. The first negative eigenvector peaks at frequency 14.5, next to the off-line zero 0.933 + 15.668i. Its Rayleigh parts are arch +4.62, pole +0.19, prime +6.43, and the prime term is what makes it negative.
- **EPS first failure width.** The width scan gives λ_min = +3.2e-4 (L = 3), −0.16 (3.25), −0.46 (3.5), −0.98 (3.75). So the Epstein Weil form is positive at widths < ~3.1 despite its off-line zeros: frequency resolution 2π/L has to reach the zero at 15.67 (NUMERICAL).
- **BEUR with x₀ = 20 (unspliced).** Without an arch term, the pole's negative direction stays uncompensated: λ_min(pole − prime) = −0.005 / −0.040 / −0.33 / −1.21 / −1.95. With ζ's arch borrowed: +0.036 / −0.31 / −1.30 / −2.59 / −3.20. The negativity at L ≤ 3 is the x₀-cutoff artifact (missing primes 2–19), so the spliced BEUR″ is the fair comparison.

## 3. Localized terms
- **EPS: the signed prime − pole channel.** That is the cancellation between the arithmetic coefficients and the pole, R_Q = ψ_Q − x, at scale x^{1/2}.
  - The arch term is not the culprit: it is more coercive than ζ's at every width.
  - The pole is identical.
  - The absolute prime block is only 1.15–1.24× ζ's.
  - What differs is how much of the prime block the pole cancels: ζ goes 7.78 → 1.56 (factor 5.0) at L = 4, while EPS goes 9.66 → 6.29 (factor 1.5). The ζ signed bound is flat (1.49, 1.48, 1.56) and the EPS one grows (1.99, 4.27, 6.29).
  - With a functional equation and gamma factor present, the defect enters only through coefficients that are not non-negative and not prime-power-supported. Those are exactly what an Euler product forces. (NUMERICAL)
- **BEUR″: no finite-width self-block term at all.** Every width-≤ 4 operator (arch, pole, prime, signed bound, λ_min W) is identical to ζ's. The maximum entry difference is exactly 0.0.
  - PROVED by construction: atoms agree below 55 > e⁴.
  - CHECKED numerically.
  - The splice works at any width: for any L₀, splicing at e^{L₀} gives a system with an Euler product and Λ ≥ 0, off-line zeros, and width-≤ L₀ data equal to ζ's.
  - The difference shows up only in the large-u field, i.e. the prime distribution beyond every fixed scale. The Cesàro mean is set by the quartet and not by W, and it grows like e^{0.2u}.
  - Without a functional equation, nothing ties that tail back to the finite-width data. That is DMV's point about missing additive structure. (PROVED, the splice argument; NUMERICAL, the field.)
- Combined:
  - **No finite-width self-block quantity, including the signed bound, can be a source of normalizability.** BEUR″ matches ζ on all of them and still has off-line zeros (PROVED).
  - **The functional equation, gamma factor, and pole together are not one either.** EPS has all three, a better arch than ζ, and off-line zeros. Its failure enters through the signed prime − pole channel (NUMERICAL).

## 4. The property (single sentence)
**A source of normalizability must use that one and the same coefficient sequence Λ(n) is non-negative and supported on prime powers (the Euler product, which Epstein lacks) and also satisfies the degree-one functional equation with Γ_ℝ(s) = π^{−s/2}Γ(s/2), which comes from Poisson summation over the additive group ℤ and is what Beurling systems lack, so that the arch term ties the prime tail at every scale to finite-width data.**

Label: **CONJECTURE** (a heuristic synthesis, not a theorem). It is not new. It is the Selberg-class intuition (Euler product + FE), and it matches Davenport–Heilbronn 1936 / Potter–Titchmarsh (FE without an Euler product allows off-line zeros) and Diamond–Montgomery–Vorhauer 2006 (an Euler product without the integers' additive structure allows zeros near σ = 1). What this round adds is where each failure enters the KG / two-space machinery:
- Epstein fails in the signed prime − pole term.
- Beurling fails in no finite-width term, only in the u → ∞ tail of the field.

## 5. Labels summary
- EPS explicit-formula data, field = zero sum: NUMERICAL (≤ 5e-5).
- EPS λ_min(W) < 0 from L ≈ 3.1, eigenvector at frequency ≈ 14.5 ↔ zero 0.933 + 15.668i: NUMERICAL.
- EPS localization in signed prime − pole: NUMERICAL.
- BEUR″ self-blocks identical to ζ for L ≤ log 55: PROVED (construction), CHECKED (0.0 difference).
- BEUR″ field = quartet: NUMERICAL.
- §4 sentence: CONJECTURE (known heuristic).
- RH: OPEN. Nothing here is a step toward RH or ARP closure. It only rules out source candidates.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
