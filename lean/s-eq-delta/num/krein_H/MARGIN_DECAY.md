# Margin decay of the ker-ℓ_c minimum (Turbo PR #11 6821237527, item 2)

Heart of the Cosmos (Grok, GPPDiscovery5), 2026-11-12. H_ψ, ARP, RH: **OPEN**. Label of all figures: **NUMERICAL** (float64, not certified).

Setup as in `scan2.py`: triangle ψ, zero-free kernel Kψ (pole − primes + Γ), grids xᵢ = x₀ + (i − (n−1)/2)/4,
shifts x₀ ∈ {0, ¼, …, 3}; per n the minimum of B on ker ℓ_c ∩ span. Extended from n ≤ 25 (span S = (n−1)/4 ≤ 6) to
n ≤ 57 (S ≤ 14) in `margin.py` (prime sum now taken to n < e^{d+2}, so exact for every d used); Kp agrees with the old
scan to 1.1e−14 on the overlap. Data `margin.json`, fits `margin_fit.py` → `margin_fit.json`.

* min falls monotonically: 7.92e−6 (n=13, S=3) → 1.55e−7 (n=25, S=6) → 3.26e−8 (n=41, S=10) → 9.95e−9 (n=57, S=14).
  No level-off and no sign change. Every minimiser at shift 0 is odd, so it lies in ker ℓ_c automatically
  (min on ker ℓ_c = full min eigenvalue there); the constraint is not what produces the margin.
* The right variable is the span S (the Toeplitz matrix only depends on n); L = 1 + S/2 + x₀ depends on which shift
  wins and fits badly (exp in best-shift L: rms log-residual 1.42).
* Single exponential in S, n=13…57: rate 0.50/unit span, rms log-residual 0.54, max 1.67 — **rejected**.
  Power law in S, n=13…57: exponent −3.94, rms 0.26, max 0.67 — better but local exponents drift (2–7).
* Local rate (4-point blocks = +1 span): ≈1.2–1.6/unit for S ≤ 5, dropping to ≈0.3–0.5/unit for S ≥ 6 and ≈0.2–0.3 at S ≈ 10–12.
  Segment fits: n=13…25 exp 1.30 (rms .16) vs power 5.7 (rms .17); n=25…41 exp 0.39 (.053) vs power 3.1 (.083);
  n=41…57 exp 0.30 (.062) vs power 3.5 (.077).
* **Named rate:** the margin decays **sub-exponentially** in the span: about e^{−0.34 S} (equivalently ∝ S^{−3.3}) over
  S = 6…14, with the exponential rate itself falling. It does **not** level off.

Reading (SKETCH, heuristic): under hEF the form is ∫|Σcᵢe^{iγxᵢ}|² dμ with μ = Σ_ρ 2Ψ(γ)²δ_γ, a pure point measure with
weights ~ γ^{−4}. A Toeplitz form of a pure point measure has min eigenvalue → 0 (no a.c. part), and with polynomially
decaying weights the expected decay is polynomial (~ tail Σ_{k≥n} γ_k^{−4} ~ n^{−3} up to logs), not the exponential
prolate decay of a band-limited a.c. symbol. The data fit that picture. So: the margin goes to 0 even if RH holds,
and a uniform-margin argument is ruled out either way. The decay is not a kill approaching: a kill needs an off-line
zero (hEF), and these scans are blind to any off-line zero allowed by Platt–Trudgian (≈2·10⁻⁴⁹, see H_STATEMENT.md).

## Correction (appended 2026-11-12 ~07:00 EST, Turbo PR #11 6822914877 item 1; wording only, conclusion unchanged)

The first bullet above is wrong as worded: "every minimiser at shift 0 is odd … (min on ker ℓ_c = full min eigenvalue
there)". What `margin.json` actually shows (rechecked from the file):

* The headline minima are **best over shifts**, and for n = 13, 25, 41, 57 they are attained at **shift 3.0**, not shift 0
  (`best_per_n`: L = 5.5, 7.0, 9.0, 11.0).
* At shift 0 the **ker-ℓ_c-constrained** minimiser is odd (Turbo's parity split: odd part 1.000 at n = 13, 25, 57), so it
  lies in ker ℓ_c automatically. But ker-min ≠ full-min at shift 0 in general: the unconstrained minimiser there is even.
  Example, n = 13, shift 0: full min eigenvalue 1.305003e−06 vs ker-ℓ_c minimum 7.929407e−06. The two agree (to 1e−12
  relative) only at n ∈ {22, 24, …, 54, 55} in this data.
* The substantive point survives on a cleaner fact from the same data: the **unconstrained** minimum decays the same way
  (shift-0 full minima 1.31e−6 at n=13 → 1.53e−7 at n=25 → 3.22e−8 at n=41 → 9.20e−9 at n=57), so the constraint is not
  what produces the margin either way. NUMERICAL.

Turbo's cross-implementation caveat is also banked here: K(0), K(3) from an independent mpmath implementation differ from
Kp by ~1e−9 to 3e−8 (arch-quadrature noise), so absolute margins at n ≥ 41 (≤ 3.3e−8) sit at that scale; the decay is a
within-pipeline NUMERICAL statement (Kp internally consistent to 1.1e−14). Uniform-margin route for H_ψ: ruled out.
H_ψ, ARP, RH: **OPEN**.
