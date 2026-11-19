# Round 19: "It has to be the natural numbers"

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

**Frame.** d = 1, s = Δ, states live in L²(ℝ>0, dx), and the shadow is Δ ↦ 1 − Δ, implemented by (JF)(x) = x⁻¹F(1/x). **RH is OPEN, and nothing here claims it.** Zeros are never used as inputs. Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN. "Cited" means a published theorem that I state but did not re-derive.

**Daniel's thesis.** ζ counts each n ∈ ℕ exactly once: the coefficients are all ones, which is a completely multiplicative sequence. Theta inversion θ(1/x) = √x θ(x) and prime factorization are one symmetry. The fakes are concocted.

**Short answer.**
- Part 1 is true and precise: "it has to be ℕ" is a theorem (Hamburger; Kaczorowski–Perelli). In d = 1 it says F_ζ is the only shadow-invariant state built from integer-frequency Gaussians.
- Part 2 sorts every fake by the naturalness property it lacks.
- Part 3 is the honest negative. The *same* uniqueness mechanism, run at conductor 5, singles out Davenport–Heilbronn, which has zeros off the line. So uniqueness alone cannot force RH. Whatever proves RH has to use the Euler product essentially, and turn an identity into an inequality (GAP-19, Part 4).

Files in this folder:
- `r19_hamburger.py/.json/.log`: Hamburger rigidity by collocation, for conductor 1 and conductor 5.
- `r19_theta.py/.json/.log`: the theta self-duality test for ℕ, a perturbed prime system, the Round-13e splice, and DMV-type Beurling primes; plus heat-flow exit from the natural class.
- `r19_table.py/.json`: coefficient checks behind the table.

---

## 1. Inversion and factorization as one symmetry on ℕ

### 1.1 Additive side: Poisson for the all-ones sequence on ℤ (PROVED, classical; CHECKED)

The measure δ_ℤ = Σ_{n∈ℤ} δ_n is its own Fourier transform (Poisson summation). Applied to e^{−πt²y}, this gives θ(1/y) = √y θ(y), where θ(y) = Σ_ℤ e^{−πn²y}.

**d = 1 form.** Let F_ζ(x) = Σ_{n≥1} e^{−πn²x²} − 1/(2x), the Round-18 lattice state. It is in L²(dx), its strip Mellin transform is ½Γ_ℝ(s)ζ(s), and it satisfies

  J F_ζ = F_ζ  (shadow invariance).

The proof is short. Write θ₊ = (θ − 1)/2. Then θ(1/x²) = xθ(x²) gives θ₊(1/x) = xθ₊(x) + x/2 − 1/2. Hence x⁻¹F_ζ(1/x) = θ₊(x) + ½ − 1/(2x) − ½ = F_ζ(x).

`r19_theta.json` (NAT row) checks the theta relation with integers up to 3·10⁴ on x ∈ [10⁻⁷, 1]. The maximum residual is 2.2·10⁻¹⁶.

### 1.2 Multiplicative side: unique factorization (PROVED)

ℚ>0 is the free abelian group on the primes, via n ↦ (v_p(n))_p, and ℕ is the free commutative monoid inside it. The all-ones sequence is the indicator function of ℕ ⊂ ℚ>0, and the indicator of a free monoid factors over its generators:

  Σ_{n∈ℕ} n^{−s} = Π_p Σ_{k≥0} p^{−ks} = Π_p (1 − p^{−s})^{−1}.

So "all ones on ℕ" and "completely multiplicative with Euler product" are the same statement.

### 1.3 Why the two are one symmetry (PROVED, Tate 1950; cited)

In Tate's thesis, ℚ sits in the adeles 𝔸_ℚ as a discrete, cocompact, self-dual subgroup, and ℤ = ℚ ∩ Π_p ℤ_p. Take the self-dual test function f = e^{−πx²} ⊗ Π_p 1_{ℤ_p}.
- The archimedean factor of its zeta integral is Γ_ℝ(s). It carries the Gaussian, and so the inversion.
- The p-adic factors are (1 − p^{−s})^{−1}. They carry the factorization.
- Poisson summation for ℚ ⊂ 𝔸 gives the functional equation Λ(s) = Λ(1 − s).

So "inversion and factorization are one symmetry" means precisely this: **ℚ is a self-dual lattice in 𝔸, and its integral points are ℕ ∪ {0} ∪ −ℕ.** Inversion is that self-duality seen at ∞. Factorization is the same lattice seen at the primes.

### 1.4 Rigidity statements (cited, with numerical replication)

**(H) Hamburger (1921), standard form.**
- Hypotheses: f(s) = Σ a_n n^{−s} and g(s) = Σ b_n n^{−s} converge absolutely for Re s > 1; (s − 1)f(s) is entire of finite order; and π^{−s/2}Γ(s/2) f(s) = π^{−(1−s)/2}Γ((1−s)/2) g(1 − s).
- Conclusion: f = g = a₁ ζ.
- No Euler product is assumed. **At conductor 1 with Γ_ℝ(s), the inversion symmetry alone forces the all-ones sequence, and with it the Euler product.** This is the precise sense in which inversion *implies* factorization.

**d = 1 corollary (PROVED, given H).** Let Φ(x) = Σ_{n≥1} a_n e^{−πn²x²} − c/(2x), with polynomially bounded a_n, and suppose JΦ = Φ. Then Φ = c F_ζ.

Proof: JΦ = Φ is equivalent to θ_b(1/y) = √y θ_b(y), where θ_b(y) = c + 2Σ a_n e^{−πn²y} and y = x². Mellin then gives Hamburger's functional equation, so b_n = c for all n. In plain words: **among shadow-invariant states built from Gaussians at integer frequencies, F_ζ is the only one.**

**NUMERICAL replication (`r19_hamburger.json`).** Collocation of θ_b(1/x) = √x θ_b(x) on x ∈ [0.35, 1] at 120 digits, with unknowns a₀…a_N.

| N | smallest log₁₀ σ | next log₁₀ σ | null vector |
|---|---|---|---|
| 6 | −26.2 | −2.2 | (1, 1, 1, 1, 1, 1 − 5·10⁻¹²) |
| 8 | −43.3 | −3.5 | all ones |
| 10 | −65.2 | −5.4 | all ones |
| 12 | −92.2 | −8.4 | all ones |

There is one isolated null direction, and it is the all-ones vector. Its separation from the next singular value grows from 24 orders of magnitude at N = 6 to 84 at N = 12.

Control: the odd class at conductor 1 (Γ_ℝ(s+1), θ(1/x) = x^{3/2}θ(x)) has no isolated null direction, at N = 8 or N = 12. There is no odd character mod 1.

**(KP) Selberg class, degree ≤ 1 (cited).**
- Conrey–Ghosh (1993): S₀ = {1}.
- Kaczorowski–Perelli (Acta Math. 182, 1999): S has no functions of degree 0 < d < 1, and S₁ = {ζ(s + iθ), L(s + iθ, χ) : χ primitive, θ ∈ ℝ}.
- They also describe the extended class S♯₁ (functional equation, no Euler product) as linear combinations of shifted Dirichlet L-functions with Dirichlet-polynomial coefficients; for a fixed conductor this is a finite-dimensional family.
- Conductor 1 together with the pole at s = 1 forces F = ζ.
- The brief also names Bombieri. I do not have a precise Bombieri theorem for this classification, so I do not attach his name to a statement.

**(C) Converse theorems (cited).** In degree 1 the converse theorem *is* Hamburger/Hecke: the functional equation forces modularity, which here means the theta function. In degree 2, Weil (1967) needs the functional equation for all twists. In higher degree there is Cogdell–Piatetski-Shapiro. In general, "FE ⇒ automorphic" needs twists. Degree 1 at conductor 1 is the case where one functional equation is enough.

**(Cór) Self-dual discrete measures (cited; context only).** Córdoba (1989) and Lev–Olevskii (Invent. Math. 2015): a measure on ℝ whose support and spectrum are both uniformly discrete is a finite combination of (translated, modulated) Dirac combs. Kurasov–Sarnak (2020) show that dropping uniform discreteness of the spectrum admits non-periodic positive Fourier quasicrystals. I did not derive "self-dual, unit masses, uniformly discrete ⇒ δ_ℤ" here, so I label that step HEURISTIC. Beurling integers are not uniformly discrete, so this route does not exclude them; Hamburger-type arguments and the numerics below do.

**Summary (PROVED, from cited theorems).** "Natural" = Dirichlet series over ℕ + Γ_ℝ functional equation at conductor 1 (+ pole at 1) ⇒ ζ (Hamburger), and the Euler product comes free. In Selberg-class language: degree 1 + conductor 1 + Euler product + pole ⇒ ζ (KP). **Natural pins ζ.**

---

## 2. Which naturalness property each fake lacks

Sources for the columns:
- "FE" means a Riemann-type functional equation, with its gamma factor and conductor.
- "deg" is the degree.
- "1s on ℕ" means a_n = 1 for all n ∈ ℕ.
- "Ram." is Ramanujan, max_p |a_p| / degree ≤ 1.
- Coefficient columns come from `r19_table.json` (n ≤ 3000); the theta column from `r19_theta.json`; the uniqueness column from `r19_hamburger.json` and theory.

| object | FE | Euler product | deg | 1s on ℕ | Ram. | unique in its own FE class | zeros off the line |
|---|---|---|---|---|---|---|---|
| ζ | yes: Γ_ℝ(s), q = 1 | yes | 1 | **yes** | yes (1.0) | **yes** (Hamburger; NUMERICAL above) | none known (RH OPEN) |
| Epstein x²+5y² (EPS) | yes: Γ_ℂ(s), q = 20 | **no**: a₆ = 2 ≠ a₂a₃ = 0 | **2** | **no** | yes (1.0) | **no**: genus theory, class contains ζ·L(χ₋₂₀) and L(χ₋₄)L(χ₅) | yes (Round 17: ρ = 0.932969697485 + 15.6682495313i) |
| Davenport–Heilbronn (DH) | yes: Γ_ℝ(s+1), q = 5, residual 10⁻¹²¹ | **no**: a₆ = 1 ≠ a₂a₃ = −0.0807 | 1 | **no** (periodic mod 5) | yes (1.0) | **yes**: NUMERICAL, null vector = DH to 10⁻¹⁰²; consistent with KP | yes (Davenport–Heilbronn 1936, cited) |
| F_θ = ζ(s−¼)ζ(s+¼) | yes: Γ_ℝ(s−¼)Γ_ℝ(s+¼), q = 1 | yes (multiplicative, defect 10⁻¹⁴) | **2** | **no**: a_p = p^{¼}+p^{−¼} | **no** (max a_p/2 = 3.77) | — | yes, by construction (ρ ± ¼); also poles at 1 ± ¼ |
| Beurling DMV-type (x₀ = 20) | **no**: theta residual 0.87 | yes | undefined | ones on its *own* integers, which ≠ ℕ (first: 20.83) | — | — | yes, by construction (quartet ½ ± 0.1 ± 28i) |
| Round-13 splice (rational p < 55 + Beurling ≥ 55) | **no**: theta residual 7.4·10⁻⁴ at x = 10⁻⁴, 2.2·10⁻² at x = 10⁻⁷ | yes | undefined | equals ℕ up to 58.94, the first Beurling prime, then not | — | — | yes (quartet), with width-≤ log 55 data identical to ζ's (Round 13e) |
| PERT control (2 → 2.01) | **no**: theta residual 4.9·10⁻³ | yes | undefined | no (first: 2.01) | — | — | not computed |

How to read the table.
- **Each fake fails exactly where its construction cut a corner.** Epstein keeps the additive symmetry (a self-dual lattice in ℝ², degree 2) but loses factorization. DH keeps degree 1 and the functional equation but sits at conductor 5, where the Euler product is not automatic, and it does not have one. F_θ keeps factorization and the functional equation but doubles the degree and breaks Ramanujan. Beurling systems and the splice keep factorization but move the integers off ℕ, and the theta inversion detects that immediately.
- **The splice deserves a note.** Its theta residual is exactly 0 to machine precision for x ≥ 10⁻², because e^{−π·58.94²·0.01} ≈ e^{−109}. The failure appears only once Gaussians wide enough to see n ≈ 59 enter. This is the d = 1 version of the Round-13e statement that the splice agrees with ζ at every width ≤ log 55. (NUMERICAL. The residual at x = 10⁻⁴ involves only integers ≲ 300, all included exactly, so it is not a truncation artifact.)

---

## 3. Does uniqueness give normalizable zeros, or RH?

### 3.1 The decisive counter-model: DH is unique too (NUMERICAL + cited theory)

Run the Hamburger collocation at conductor 5, odd parity, real coefficients, and root number +1: θ(x) = Σ a_n n e^{−πn²x/5}, constrained by θ(1/x) = x^{3/2}θ(x).

| N | smallest log₁₀ σ | next | max deviation of null vector from DH, n ≤ 12 |
|---|---|---|---|
| 13 | −22.7 | −5.9 | 1.6·10⁻⁴ |
| 18 | −41.2 | −8.6 | 4.6·10⁻²⁵ |
| 22 | −60.9 | −11.5 | 1.0·10⁻⁴³ |
| 27 | −91.1 | −16.5 | 1.0·10⁻⁷² |
| 34 | −120.1 | −27.8 | 5.9·10⁻¹⁰² |

There is one isolated null direction, and it is **exactly the Davenport–Heilbronn coefficient vector**. So DH is as rigid in its class as ζ is in its own: unique up to scale. Theory says the same. By KP, S♯₁ at conductor 5 is spanned by L(s, χ) and L(s, χ̄) with χ odd mod 5; the real, root-number-+1 combination is one-dimensional, and that line is DH.

DH has infinitely many zeros off the critical line (cited). Therefore:

> **PROVED (logic), resting on NUMERICAL uniqueness plus the cited KP structure:** "unique object with its inversion symmetry" does **not** imply "zeros on the line."
>
> The rigidity argument that pins ζ at conductor 1 pins DH at conductor 5. The only difference between the two outcomes is the Euler product. It is automatic at q = 1 and absent at q = 5.

Consequence. A proof route through "it has to be ℕ" cannot run through uniqueness. It must use multiplicativity (the Euler product) as an active ingredient. That is exactly GRH territory: the Euler product plus degree 1 gives the Dirichlet L-functions, whose RH is also OPEN.

(Remark. Hecke's theorem that λ > 2 gives infinite-dimensional spaces concerns full Fourier series in e^{2πimz/λ}. These thetas contain only the square exponents m = n², and the numerics confirm that the class here is one-dimensional.)

### 3.2 Variational or extremal characterization (HEURISTIC; no mechanism found)

The natural class H₁ is a single point, so there is nothing to vary inside it. Any extremal principle needs an ambient family. Three were tried:

1. **Heat flow (de Bruijn–Newman).** H_t(z) = ∫Φ(u)e^{tu²}cos(zu) du.
   - N3 (NUMERICAL): H_t/H₀ is not constant for t = ±0.1. At z = 0, 4, 8, 12 the ratio is 0.9954, 0.9986, 1.0106, 1.0573 for t = −0.1, and 1.0047, 1.0014, 0.9892, 0.9430 for t = +0.1. So by Hamburger, H_t leaves the natural class for every t ≠ 0. **The heat line meets "natural" at exactly one point, t = 0.**
   - Separately, all zeros of H_t are real exactly when t ≥ Λ. Rodgers–Tao (2020) prove Λ ≥ 0, and Platt–Trudgian prove Λ ≤ 0.2 (both cited). RH ⇔ Λ = 0.
   - So RH is the statement that **the natural point coincides with the reality threshold**. Rigidity locates the first point; nothing in rigidity locates the second. Λ ≥ 0 says ζ is "at least critical". That is suggestive of a criticality principle (HEURISTIC), but the deformation leaves the class, so rigidity has nothing to say along it.
2. **Conductor family.** At q = 5 the rigid point (DH) is not on the line. So "rigid point ⇒ critical" fails as a general principle (§3.1).
3. **Shadow-symmetric states in d = 1.** F_ζ is the unique J-invariant integer-Gaussian state (§1.4). Its overlaps with the normalizable primaries x^{−Δ}, Re Δ = ½, are ½Γ_ℝ(Δ)ζ(Δ). Uniqueness of the *state* says nothing about where its *overlap function* vanishes.

   The DH analogue makes this concrete. The analogous unique J-invariant state at conductor 5 (odd: Σ a_n n x e^{−πn²x²/5}) has overlaps that vanish off the line, exactly as the Round-17/18 ghosts did. Normalizability still holds only on Re Δ = ½ (Pasterski–Shao (4.7)). Zero-normalizability remains the same OPEN question as in Rounds 17–18.

### 3.3 Positivity: the naturalness positivity is the wrong kind (PROVED as a comparison; the transfer is OPEN)

Naturalness comes with three positivities: δ_ℤ ≥ 0; δ_ℤ is positive-definite, since δ̂_ℤ = δ_ℤ ≥ 0; and Λ(n) ≥ 0. All of them live on the additive group (ℝ, +) or on the prime side alone.

RH is equivalent to Weil positivity: positive-definiteness, on the multiplicative group ℝ>0, of W = pole + arch − prime, with the prime part *subtracted*. The two fakes show these positivities do not transfer:
- Epstein has the additive positivity (δ_L ≥ 0, positive-definite) and fails Weil positivity.
- F_θ has Euler product, functional equation and Λ_F(n) ≥ 0, and fails Weil positivity (degree 2; earlier rounds).

The map x ↦ log x that links the two groups does not preserve positive-definiteness.

### 3.4 Verdict on Task 3

Uniqueness does **not** yield "zeros normalizable" or RH. That is shown by DH (§3.1), not just conjectured. No variational or deformation mechanism found here closes the gap. The Λ = 0 criticality reading is a faithful *restatement* of RH, not a route to it.

---

## 4. GAP-19, in plain words

**GAP-19 (positivity transfer).** Show that the self-duality of ℤ in ℝ (Poisson: δ_ℤ is its own Fourier transform, and is positive and positive-definite), *together with* the fact that ℕ is the free monoid on the primes (the Euler product, used essentially), forces positive-definiteness of the Weil distribution on the multiplicative group ℝ>0.

In adelic words: the single symmetry "ℚ is self-dual in 𝔸" must produce an **inequality** on the idele class group, not only the **identity** it already gives, which is the functional equation.

Plainly:
- "It has to be ℕ" is true, and it is a theorem. It explains exactly why every fake is fake.
- But what it delivers is an equation (inversion) and a product (factorization), and both are identities. RH is an inequality.
- The uniqueness of ℕ cannot by itself be the inequality. At conductor 5, the very same uniqueness picks out Davenport–Heilbronn, which has zeros off the line. The only thing ζ has that DH lacks is the Euler product, so the inequality must come from multiplication, not from uniqueness.
- Over function fields that inequality comes from a surface (C × C) and the Hodge index theorem. For ℕ no such surface is known. This is the same missing object as the Round-15/16 Frobenius gap.

## 5. Status ledger

- PROVED (classical or cited): Poisson ⇔ θ inversion ⇔ J F_ζ = F_ζ; all-ones = Euler product; Tate's one-symmetry statement; Hamburger and the d = 1 corollary; KP degree-1 classification; Rodgers–Tao Λ ≥ 0.
- PROVED (logic, here): uniqueness ⇏ RH, given DH's uniqueness in its class (NUMERICAL here, consistent with KP) and DH's off-line zeros (cited).
- NUMERICAL:
  - Hamburger null space at q = 1 is all ones (gap of 24 to 84 orders);
  - the q = 5 odd null space is DH (to 10⁻¹⁰²);
  - theta residuals of 2·10⁻¹⁶ (ℕ), 4.9·10⁻³ (PERT), 2.2·10⁻² (splice), 0.87 (DMV20);
  - heat-flow exit from the class.
- CHECKED: the coefficient table (r19_table.json).
- HEURISTIC: Λ = 0 as "the natural point is critical"; the Córdoba-type "self-dual ⇒ ℤ" step.
- GAP-19: positivity transfer, as above. OPEN: RH.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
