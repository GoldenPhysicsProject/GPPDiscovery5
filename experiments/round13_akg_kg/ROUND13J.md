# Round 13j — the golden unit: K = Q(√5), ring Z[φ], unit φ, regulator log φ, discriminant 5

Code: `golden.py`, `zero_diffraction.py`. Data: `golden.json`, `golden_T_atoms.json`, `zero_diffraction.json`. Zeros are test data only.

## 1. Weil form of ζ_K and the unit circle (Tate)
- **Completed zeta.** Λ_K(s) = 5^{s/2} Γ_ℝ(s)² ζ_K(s) = Λ_K(1−s), with ζ_K = ζ·L(s, χ₅).
  - −ζ_K′/ζ_K = Σ Λ_K(n) n^{−s}, where Λ_K(n) = Λ(n)(1 + χ₅(n)) ≥ 0. Split p contribute 2 log p, inert p contribute only at p² (2 log p), and 5 contributes log 5.
- **Weil form.** W_K(g) = [2(ψ(1/4) − log π) + log 5] A(0) + 2∫ 2w(t)[A(0) − A(t)] dt + pole(ζ) − 2 Σ Λ_K(n) n^{−1/2} A(log n), with w(t) = e^{−t/2}/(1 − e^{−2t}).
  - The two real places give two copies of the ζ arch term, and the discriminant gives log 5.
  - **Exactly W_K = W_ζ + W_{χ₅}**, where W_{χ₅} has arch constant ψ(1/4) − log π + log 5 and no pole. Galerkin check: |W_K − W_ζ − W_χ| ≤ 2e-13 (PROVED, CHECKED).
- **Validation (NUMERICAL).** The ζ_K field built from these data has Cesàro mean / W = 0.9994 at L = 1 and 0.9999 at L = 2 (u ≤ 13). This checks the arch constant, the factor 2w, and the pole term.
- **Tate compactification.**
  - C_K = 𝕀_K/K^× ≅ C_K^1 × ℝ_{>0}, and C_K^1 is compact.
  - Its archimedean norm-one part {(x, y) : |xy| = 1}/±φ^ℤ is a circle of length R = log φ.
  - **This compactifies the norm-one direction, not the norm (scale) direction ℝ_{>0}.** The scale direction is where the e^{±u/2} pole weights and the normalizability question live, and it stays non-compact (PROVED, structural).
  - Fourier analysis on the circle splits the Weil form on C_K into a sum over Hecke characters with angular infinity type |α/α′|^{πik/log φ}. Positivity on C_K ⇔ GRH for every such Hecke L-function. The k = 0 component is W_K.
- **Where log φ enters.** Only in the residue of ζ_K at s = 1 (2 log φ/√5 = 0.4304) and in the spacing π/log φ of the angular characters. **It does not appear anywhere in W_K.** The explicit formula uses only the order of the pole, not its residue (PROVED).
- **Two-space numbers for ζ_K, L = 0.5 / 1 / 2 / 3 / 4 / 6 / 8** (NUMERICAL):
  - λ_min(W_K) ≥ 0 at every width, at the float floor from L = 2 on.
  - λ_min(arch_K) = −0.19 / −1.80 / −3.58 / −4.71 / −5.52 / −6.60 / −7.26. That is about twice ζ's, as two real places should give.
  - Signed λ_max(prime − pole) = 0.005 / 0.04 / 1.23 / 2.52 / 4.07 / 4.51 / 4.80.
    - This is roughly ζ's value (≈ 1.5–1.6) plus the χ₅ part (≈ 3.2, flat), since χ₅ has no pole to cancel against.
    - The two-space exponent is unchanged: κ = 1/2.
- **Verdict on 1: WRONG as a source.** The compact circle of length log φ gives no positivity or normalizability that ζ lacked, because it acts on the wrong factor.
  - What the unit group does do is make Hecke's theta-integral proof of the functional equation for K converge. That is an input to the functional equation, and the functional equation is not positivity: Epstein has one too.

## 2. Does positivity for ζ_K give positivity for ζ?
- W_K ≥ 0 for all g ⇔ RH for ζ_K ⇔ RH(ζ) and GRH(L(χ₅)), because ζ_K's zeros are the union of the two zero sets (Weil's criterion applied to ζ_K) (CHECKED). So it does imply RH(ζ), but from a strictly stronger hypothesis.
- For a single g, W_K(g) ≥ 0 does not give W_ζ(g) ≥ 0. You would need W_χ(g) ≤ W_K(g), and the sign of W_χ is itself GRH-hard.
- The L(χ₅) part does not spoil anything numerically: λ_min(W_χ) = +0.71 at L = 0.5, then at the float floor. It just adds a second, equally hard positivity problem.
- **Honest verdict: nothing is gained.** It is the sum of two Weil forms, as stated.

## 3. Falsifier over Z[φ]
- **No class-group Epstein exists with this unit.** Q(√5) has narrow class number 1 (N(φ) = −1). The norm-form "Epstein" Σ_{α mod units} |N α|^{−s} of x² + xy − y² is ζ_K itself (CHECKED).
- **The Epstein of Z[φ]'s own trace form** is Tr((x + yφ)²) = 2x² + 2xy + 3y², which has discriminant −20. **It is exactly the non-principal class of 13e** (genus theory; it contains L(χ₅)):
  Z_T = ζ L(χ₋₂₀) − L(χ₋₄) L(χ₅).
  - Its smallest norm is 2, so a(1) = 0. Normalized as T′(s) = 2^s Z_T/1, it is a generalized Dirichlet series in (m/2)^{−s}.
  - T′ has a functional equation with Γ(s) and conductor 5: Λ′(s) = (√5/2π)^s Γ(s) T′(s) = Λ′(1−s).
  - −T′′/T′ coefficients come from its own representation numbers. Checks: at s = 6 and 8 they agree with mpmath to 1e-15, and Σ r_T n^{−3} matches the L-function value to 6e-14 (CHECKED).
  - Its zeros reach far right. findroot gives 1.994 + 23.020i (test data). The coefficients grow geometrically, e.g. c((3/2)^k) ~ (−2)^k.
- **Results at L = 0.5 / 1 / 2 / 3 / 4 / 8** (NUMERICAL):
  - λ_min(W_T) = +0.87 / **−0.07** / −5.0 / −39 / −187 / −1.1e5.
  - Arch λ_min is ordinary: 0.05 / −1.33 / −2.68 / −3.44 / −3.95 / −4.97.
  - |prime| λ_max = 0.66 / 1.62 / 9.1 / 53 / 233 / 1.8e5. Signed prime − pole is 0.66 / 1.11 / 6.5 / 42 / 189 / 1.1e5.
  - **Localized term: the signed prime − pole channel**, i.e. the non-Euler coefficients. This is the same as 13e, only much stronger. The arch term and the functional equation are fine.
- **Fake prime system with the same unit.** I spliced ζ_K's Λ_K atoms below 55 with Beurling atoms from 55 on, keeping K's arch, discriminant, and implicitly its unit, which the Weil form never sees.
  - It is identical to ζ_K for L ≤ 4 (λ_min at the float floor).
  - It fails right after the splice: −0.71 at L = 4.5, −0.89 at 5, −1.17 at 6, −3.25 at 8.
  - **Localized term: no finite-width term.** The difference is only in the tail beyond the splice, the same as 13e. The unit does not help, because W_K contains no unit data.

## 4. Dyson's quasicrystal idea and the Fibonacci set
- **The zero set as a "crystal."** By the Guinand–Weil formula, μ = Σ_γ δ_γ minus its smooth density (1/2π) log(|t|/2π) has Fourier transform supported on {0} ∪ {±log pᵏ}, with weights −Λ(n)n^{−1/2}.
  - Numerical check (4519 zeros, γ < 5000, Gaussian taper T = 1500): peaks at log 2, 3, 4, 5, 7 are −139, −180, −98, −204, −209 in units of the background rms.
  - Measured peak ratio log 3 / log 2 = 1.2957, against the predicted (Λ(3)/√3)/(Λ(2)/√2) = 1.2941 (NUMERICAL).
- **It is not a Fourier quasicrystal in the strict sense.**
  - The density grows like log t, so the support is not uniformly discrete. The spectrum {log n} accumulates.
  - With off-line zeros, μ is not even a measure on ℝ: RH is needed just to state the criterion.
  - Lev–Olevskii (Invent. Math. 200, 2015) prove that uniformly discrete support and spectrum with positive masses force periodicity. Olevskii–Ulanovskii (CR Math. 2020) and Kurasov–Sarnak (2020) classify unit-mass Fourier quasicrystals as real-rooted exponential-polynomial zero sets, which have bounded density.
  - **ζ's zeros fall outside these classified classes.** Dyson's program needs a generalized class, and it remains OPEN. Reference: Lagarias, "Mathematical quasicrystals and the problem of diffraction" (2000).
- **The Fibonacci cut-and-project set** {m + nφ : m + nφ′ ∈ window} has a pure-point, positive diffraction supported on a **rank-2** module (1/√5)Z[φ].
  - ζ's "diffraction" is supported on log ℕ. The log p are linearly independent over ℚ, so this module has **infinite rank**.
  - The correct cut-and-project picture for ζ is Bohr's lift to the infinite torus 𝕋^{primes}, not ℤ² with φ.
- **Is there a peak at log φ?** None: |D| at log φ, 2 log φ, 3 log φ, log √5 is 0.37, 0.72, 0.97, 0.64 rms, i.e. background (NUMERICAL). **φ has no role here.** Its "spectrum" is a rank-2 module, while ζ's is rank-∞ and multiplicative.

## 5. Numerology flags
- **log φ = 0.4812 ≈ 1/2.** This is a coincidence. κ = 1/2 in 13b is exact: it comes from 2 sinh(L′/2), i.e. the pole at s = 1. The slopes 0.48–0.55 measured in 13b scatter around 1/2, not around log φ. **FLAG: numerology.**
- **"φ^{n+1} = φ^n + φ^{n−1} couples addition and multiplication."** True, but this is just the minimal polynomial x² = x + 1 of a unit. Every algebraic integer satisfies such a linear recurrence, e.g. (1+√2)^{n+1} = 2(1+√2)^n + (1+√2)^{n−1}, and the ring ℤ already couples + and ×.
  - The coupling ζ needs (13e) is the ring ℤ itself: Euler product over its primes plus Poisson summation over ℤ. Z[φ] has the same structure as a ring, and the unit only enters the proof of the functional equation.
  - **FLAG: right as intuition about ring structure, wrong as a specific role for φ.**
- **Physics first.** The only physical role the unit plays is compactifying a direction (the norm-one circle), and that is not the direction where normalizability fails.

## Labels
- W_K = W_ζ + W_χ: PROVED, CHECKED.
- ζ_K numerics: NUMERICAL. Field validation: NUMERICAL.
- Tate and unit-circle statements: PROVED (structural).
- §2: CHECKED.
- Trace-form Epstein coefficients: CHECKED. Its Weil numerics: NUMERICAL.
- Splice: PROVED (construction) + NUMERICAL.
- Diffraction: NUMERICAL. Quasicrystal classification: LITERATURE.
- Dyson program: OPEN. RH: OPEN.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
