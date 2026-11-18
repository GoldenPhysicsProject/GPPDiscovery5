# Round 14A — adelic Pasterski–Shao bulk: Tate / Weil representation and the Weil form

Lane: `experiments/round13_akg_kg/principal_series/` (branch `team/principal-series`). Code: `r14a_adelic_ps.py`. Data: `r14a_adelic_ps.json`. Counter-model data reused from Round 13e (`../ROUND13E.md`, `../epstein_zeros*.json`, `../beur_spliced.json`). Zeros are used only as test data, never as proof inputs. Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN.

**Question (from Heart, Round 13i follow-up).** Build the bulk from the adeles of ℚ instead of the primon Fock space. Use the Tate representation of GL₁(𝔸) on S(𝔸) and L²(𝔸), the Connes quotient L²(C_ℚ), and the Weil representation of SL₂(𝔸), which carries Poisson summation. Mellin-transform the "plane waves" along |x| into weight-Δ states. Does the inherited inner product equal the Weil form, with zeros as Δ values, so that unitarity forces Re Δ = 1/2? If not, where exactly does it fail? Run every step against the Beurling and Epstein counter-models, plus the Davenport–Heilbronn function, which turns out to be the sharpest in-framework test.

**Answer in one line.** Unitarity gives Re Δ = 1/2 as the *whole continuous unitary axis*, kinematically from Haar measure. The zeros are not in the inherited norm; they show up only as common zeros of overlaps (cokernel/absorption, and only in a Sobolev-weighted space). Positivity of the Weil form is an **extra premise**: the prime operator must be dominated by arch + pole. That premise is Weil's criterion itself. Every step before it holds for at least one counter-model with off-line zeros, and the Davenport–Heilbronn function passes **every** adelic step over ℚ.

---

## 1. Setup and the Mellin ("conformal primary") states

- Let π(a)f(x) = |a|^{1/2} f(ax) for a ∈ 𝔸^× acting on L²(𝔸, dx). The action is unitary because d(ax) = |a|dx. (PROVED, standard.)
- **Local Mellin states.** At each place, dx_v = |x|_v d^×x_v, so L²(ℚ_v, dx) = L²(ℚ_v^×, |x| d^×x). Mellin–Plancherel then decomposes it over unitary characters χ|·|^{s−1/2} with **Re s = 1/2**. The local weight-Δ "plane wave" is |Δ, χ⟩_v = χ(x)|x|^{Δ−1}, and ⟨Δ,χ|Δ′,χ′⟩_v = δ_{χχ′}·δ(Im Δ − Im Δ′) (up to 2π, Haar normalisation), exactly diagonal on Re Δ = 1/2. (PROVED.)
- **The ideles are a null set in the adeles.** ∏_p vol(ℤ_p^×) = ∏(1 − 1/p) = 0. So L²(𝔸) does **not** decompose directly over idele-class characters, and a global Mellin state on L²(𝔸) does not exist. You have to pass to C_ℚ = 𝔸^×/ℚ^× (Connes) through E(f)(a) = |a|^{1/2} Σ_{q∈ℚ^×} f(aq). (PROVED for the null-set statement. This is the reason the quotient is forced, not chosen.)
- **Tate.** Z(f, χ|·|^s) = ∫_{𝔸^×} f(x)χ(x)|x|^s d^×x = ∫_{C_ℚ} E(f)(a) χ(a)|a|^{s−1/2} d^×a. For f₀ = e^{−πx²} ⊗ 1_ℤ̂ this is Λ(s) = π^{−s/2}Γ(s/2)ζ(s). (PROVED, Tate 1950.)
  - Theta-integral evaluation agrees with Λ(s) to rel. 1.3e-22 at s = 2.3+1.7i and to 4e-8 at 1.5+6i, the quadrature limit there. (NUMERICAL, `tate_vs_Lambda`.)
- **Relation to the primon gas (13i).** ℓ²(ℕ) ≅ ℓ²((𝔸_f^× ∩ ℤ̂)/ℤ̂^×): integral finite ideles mod units are ℕ by unique factorisation. The primon Fock space is therefore the finite-idele, integral-cone shadow of this construction. (PROVED.)

## 2. (1) What inner product is inherited

- **On L²(C_ℚ), the inherited norm is Plancherel, with no arithmetic.** Write C_ℚ = ℝ_{>0} × ℤ̂^×. Then ‖φ‖² = Σ_χ ∫_{Re s=1/2} |φ̂(χ, s)|² |ds|/2π, where χ runs over characters of ℤ̂^×. Mellin states are diagonal, δ-normalisable, and exactly on Re Δ = 1/2. (PROVED.)
- **Comparison with 13i.** The PV off-diagonal term and the "π·Res" weight both disappear. That is an improvement in diagonality. The cost is that the line now comes from Haar measure alone, not from the pole at 1, so it is *even more* kinematic than in 13i. Any locally compact abelian group with a modulus has this line. (PROVED.)
- **On test vectors E(f), f ∈ S₀(𝔸)** (f(0) = f̂(0) = 0), the norm is
  ‖E f‖² = Σ_χ ∫ |Z(f, χ|·|^{1/2+iτ})|² dτ/2π.
  Here ζ (or L(χ, ·)) enters only as a *multiplier* |Λ(1/2+iτ)|² on the unitary axis. (PROVED, Plancherel.)
  - CHECKED for f_∞ = g − 2g(2·), g = x²e^{−πx²}, f_p = 1_{ℤ_p}. Both sides equal 0.016180078160874 (rel. err 2e-31). The Mellin transform of E f at 1/2+3i matches ζ(s)(1−2^{1−s})π^{−(s+2)/2}Γ(s/2+1) to 12 digits. (NUMERICAL.)

## 3. (2) Does this inner product equal or dominate the Weil form? No.

- **Weil's quadratic form.** For h = φ * φ̃ on C_ℚ, the explicit formula reads
  W(φ) = ĥ(0) + ĥ(1) − Σ_n Λ(n) n^{−1/2}[h(n) + h(1/n)] + ∫ |φ̂(1/2+iτ)|² w(τ) dτ,
  with w(τ) = (1/2π)[Re ψ(1/4 + iτ/2) − log π].
  The spectral side is Σ_ρ ĥ(ρ). (PROVED, Weil 1952 in the Connes / Bombieri normalisation.)
- **W is not the Plancherel norm.** W is that norm reweighted by w(τ), plus a rank-2 pole term, minus a prime convolution operator P supported at log n.
- **The weight w is neither positive nor bounded.** w(0) = −0.855, w(5) = −0.037, w(10) = +0.074, w(100) = +0.440, w(10⁴) = +1.173, growing like (1/2π) log τ. So W/‖φ‖² is unbounded above on high-frequency packets, and the Plancherel norm cannot dominate W. (NUMERICAL, `arch_weight`; PROVED asymptotically from Stirling.)
- **Positivity of the inherited norm says nothing about the sign of W.** ‖E f‖² > 0 holds for every f in every counter-model below, and W is negative for some of them (13e: Epstein λ_min(W) = −1.56 at L = 4). (NUMERICAL, from 13e.)
- **Exact failure point.** W ≥ 0 ⇔ P ≤ A + Pole as quadratic forms on the relevant test space, where A is the arch multiplier w. The adelic construction supplies A (the archimedean local gamma factor), the pole (Poisson / f̂(0)), and P (the finite places). It supplies no inequality between them. **That inequality is Weil's criterion verbatim**, so it is the extra positivity premise. (PROVED as a reduction; the inequality itself is OPEN.)

## 4. (3) Do zeros appear as Δ values, and does unitarity force Re Δ = 1/2?

- **Unitarity places the whole continuous spectrum on Re Δ = 1/2.** The generator of the ℝ_{>0}-scaling on L²(C_ℚ) has purely continuous spectrum, all of the line. This is automatic and has nothing to do with the zeros. (PROVED.)
- **In plain L² the zeros are completely invisible: the closure of E(S₀) is all of L².** (PROVED, sketch.)
  - Take the unramified χ. Suppose h is orthogonal to every E(f). Then ∫ ζ(1/2+iτ) G(τ) h̄(τ) dτ = 0 for a dense class of Mellin transforms G of even C_c^∞(ℝ^×) functions with ∫f = 0.
  - Hence ζ(1/2+iτ)·h̄ = 0 a.e. Since ζ ≠ 0 a.e. on the line, h = 0.
  - The critical zeros are a null set. At the first one (test data) the multiplier is |M| = 5.8e-35, but that has no L² consequence.
- **Connes' remedy.** In a Sobolev-weighted L²_δ (weight (1 + log²|x|)^{δ/2}), point evaluation on the line becomes continuous. The cokernel of E then carries the **critical** zeros as absorption spectrum. Off-line zeros enter only as **resonances**, which no Hilbert-space unitarity argument can exclude. Connes reduces RH to the validity (positivity) of the global trace formula, which is equivalent to W ≥ 0. (CHECKED against the abstract of Connes, Selecta Math. 5 (1999), arXiv:math/9811068: "critical zeros … as an absorption spectrum, while eventual noncritical zeros appear as resonances … reduces the Riemann hypothesis to the validity of the trace formula". The detailed statements are not re-derived here.)
- **Verdict on (3).** Zeros are never Δ-eigenvalues of a unitary generator in this bulk. They are common zeros of the overlaps ⟨Δ|E f⟩ = Z(f, |·|^Δ). Unitarity constrains only the continuous principal series. Turning zeros into point spectrum of a self-adjoint generator would be a Hilbert–Pólya realisation. (GAP / OPEN.)

## 5. GL₂(𝔸) / Weil representation: does Poisson add positivity?

- **What the Weil representation r of SL₂(𝔸) on S(𝔸) does.** (PROVED, standard.)
  - It is unitary on L²(𝔸).
  - It contains the dilations: r(diag(a, a^{−1}))f(x) = χ_ψ(a)|a|^{1/2} f(ax). So the GL₁ Tate action of §1 sits inside it.
  - Its Weyl element acts as the Fourier transform.
  - θ(f) = Σ_{q∈ℚ} f(q) is SL₂(ℚ)-invariant, which is Poisson summation over ℚ ⊂ 𝔸. This is the additive structure that 13e/13i found missing from the Fock space.
- **What Poisson contributes is an identity, not an inequality.** θ(r(w)f) = θ(f) yields the functional equation Λ(s) = Λ(1−s) and the pole terms. It yields no new positive form. The Plancherel norms of r are again Haar norms, and θ is a distribution, not an L² vector. (PROVED.)
- **Theta-lift inner products** (Rallis inner-product formula, Siegel–Weil) are positive. They evaluate to L-values at edge points (s = 1, 2, …), not to sums over zeros. (CHECKED, literature recall, not re-derived.)
- **Epstein x²+5y² lives in exactly this framework.**
  - Take the Weil representation of SL₂(𝔸_ℚ) on S(𝔸²) for Q = x² + 5y², with f = Gaussian_Q ⊗ 1_{ℤ_p²}.
  - θ_Q(f) is SL₂(ℚ)-invariant by Poisson on ℤ², and its Mellin transform is (√20/2π)^s Γ(s) Z_Q(s). (PROVED, definition plus Poisson.)
  - Z_Q = ζ·L(χ₋₂₀) + L(χ₋₄)L(χ₅) matches the direct lattice sum at s = 2: 2.5024222413486568 for both. (NUMERICAL.)
  - The off-line zero 0.932969697485 + 15.6682495313i satisfies |Z_Q| = 1e-30. (NUMERICAL, test data.)
- **Conclusion for GL₂.** Adding the Weil representation adds Poisson but not positivity, and Epstein passes it. (PROVED for "passes".)

## 6. Counter-models step by step

Steps: S1 = unitary scaling and kinematic line Re Δ = 1/2. S2 = Mellin/Plancherel norm ‖E f‖² = ∫|Z|² ≥ 0. S3 = zeros invisible in L², and only critical zeros detected in L²_δ. S4 = local Tate gamma factor, local Fourier unitarity, and Poisson / functional equation. S5 = test vector isotypic for C_ℚ, i.e. an Euler product. S6 = P ≤ A + Pole (Weil).

| step | ζ | Beurling (spliced BEUR″, 13e) | Epstein x²+5y² | Davenport–Heilbronn (DH) |
|---|---|---|---|---|
| S1 | holds | holds (Mellin on ℝ_{>0} of Σ_{n∈ℕ_B} f(tn)) | holds | holds |
| S2 | holds | holds (multiplier ζ_B) | holds (multiplier Z_Q) | holds (multiplier DH) |
| S3 | holds | holds (ζ_B ≠ 0 a.e. on the line) | holds | holds |
| S4 | holds | **fails**: no additive group, no Poisson, no gamma factor | holds (Poisson on ℤ², Γ(s)) | **holds** (adelic over ℚ; FE ξ(s) = ξ(1−s) checked to 1e-30) |
| S5 | holds | holds (Euler product, Λ_B ≥ 0) | **fails** (no Euler product) | **fails**: c(2)c(3) = −0.0807 ≠ c(1) = 1 |
| S6 | OPEN | fails (quartet 1/2 ± 0.1 ± 28i) | fails (13e: λ_min W = −1.56 at L = 4) | fails (off-line zero 0.808517182457 + 85.6993484854i, |DH| = 2e-30) |

- **Labels.** S1–S3 rows are PROVED (one-line Mellin/Plancherel, valid for any Dirichlet series with a pole of the right kind or none). The S4/S5 entries are PROVED (structural) and CHECKED numerically where values are quoted. The S6 failures are NUMERICAL (test data).
- **The DH function is an exact GL₁(𝔸_ℚ) Tate integral.**
  - Construction: Z(f, |·|^s) = (arch)·DH(s), with f = x e^{−πx²} ⊗ f₅ ⊗ ⊗_{p≠5} 1_{ℤ_p} and f₅ = Σ_a c(a) 1_{a+5ℤ₅}.
  - Coefficients: c = (1, κ, −κ, −1) for a = 1, 2, 3, 4, with κ = 0.2840790438.
  - The only thing that differs from ζ is that f₅ mixes two characters of ℤ₅^× (χ, χ̄ with χ(2) = i) instead of lying in one isotypic component. (PROVED by construction; FE and zero NUMERICAL.)
- **Reading the table.** A step that holds for a counter-model with off-line zeros cannot force Re Δ = 1/2.
  - S1–S3 hold for all three. That is everything unitarity gives.
  - S4 (Poisson / additive ℤ, including GL₂) holds for Epstein and DH.
  - S5 (Euler product / isotypy) holds for Beurling.
  - So no single structural step S1–S5 excludes all counter-models. (PROVED from the table.)
  - The conjunction S4 ∧ S5 excludes all three. But it is a *selection* condition, not a positivity mechanism. Nothing in the adelic bulk converts S4 ∧ S5 into S6. That conversion is exactly the GAP.

## 7. Refinement of the 13e/13i sentence

- **13e/13i.** Positivity needs "the ring ℤ with both" multiplication and addition.
- **DH sharpens this within the adeles.** DH is built from ℤ, has Poisson over ℚ, the degree-one FE, and an odd gamma factor, and it still has off-line zeros.
- **Adelic form of the requirement:** additive Poisson over ℚ ⊂ 𝔸 **and** a test vector that is a C_ℚ-eigenvector (Hecke character isotypy = Euler product), at **every** finite place simultaneously. Isotypy at a single place fails: DH differs from L(s, χ) only at p = 5. (CONJECTURE, heuristic: the Selberg-class intuition restated adelically.)
- Even granted that requirement, the passage to W ≥ 0 is not supplied by any representation-theoretic positivity found here. (GAP.)

## 8. Labels summary

- Unitary Tate action, local Mellin states diagonal on Re Δ = 1/2, ideles null in 𝔸, primon = integral finite-idele shadow: PROVED.
- Tate integral = Λ: PROVED (Tate) / NUMERICAL (≤ 4e-8).
- Plancherel ‖E f‖² = ∫|Z(f, 1/2+iτ)|²: PROVED / NUMERICAL (2e-31).
- Plancherel norm neither equals nor dominates W; arch weight sign-changing and unbounded: PROVED (structure, Stirling) / NUMERICAL.
- W ≥ 0 ⇔ P ≤ A + Pole: PROVED (reduction), and this is exactly the missing premise (Weil's criterion). The inequality itself: OPEN.
- Zeros invisible in L²(C_ℚ) (closure of E(S₀) = L²): PROVED (sketch). Critical zeros as absorption in L²_δ, off-line zeros as resonances: CHECKED (Connes 1999 abstract).
- GL₂ / Weil representation adds Poisson but no positive form; Epstein lies in it: PROVED / NUMERICAL.
- DH as exact GL₁(𝔸_ℚ) Tate integral with FE and off-line zero: PROVED (construction) / NUMERICAL.
- Counter-model table: as labelled in §6.
- Adelic requirement in §7: CONJECTURE. Conversion to positivity: GAP.
- RH: OPEN. Nothing here closes RH, ARP or H_ψ. This round only rules out the adelic Pasterski–Shao bulk as a positivity source and pins down where it stops.

— Daniel Toupin / Golden Physics Project. This seat: Grok Bot (principal-series lane). Tasking: Heart of the Cosmos. Audit: Turbo.
