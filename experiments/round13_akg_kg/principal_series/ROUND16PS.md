# Round 16 PS: is there a representation-theoretic replacement for polarization / Hodge index over ℤ?

Principal-series lane (GPPDiscovery5, branch `team/principal-series`), 2026-11-18 (Toronto). Question from Daniel via Heart of the Cosmos, after Round 15 (`experiments/round13_akg_kg/round15_hasse/ROUND15.md`, main 4dbeb9d).

Location: Heart asked for `principal_series/ROUND16PS.md`. I used the existing lane folder `experiments/round13_akg_kg/principal_series/`, next to ROUND14A/ROUND14F.

**Labels:** PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN.
**Status:** RH, ARP and H_ψ are **OPEN**. Nothing here proves RH. No zeros are used as inputs anywhere: the script uses only the integer coefficients τ(n) of Δ, lattice counts for x² + 5y², and the DH coefficients. Zeros appear only inside explanatory identities.

Notation follows Round 15: F, F̂ = Frobenius and its dual (Rosati adjoint); a = F + F̂, q = F·F̂; ω_i = α_i/√q; the Toeplitz Hasse Gram T[j,k] = Σ_i ω_i^{k−j}; the one-shift Weil block [g, τ_s g].

## Summary

1. **Over F_q, Hodge index is a unitarity statement in a finite-dimensional space.** Rosati positivity makes U = F/√q an element with U·U′ = 1 in an algebra with a positive involution. The Round 15 Toeplitz Gram is exactly the Gram matrix of the powers U^k for the Rosati trace form. PROVED (standard).
2. **The trace = sum, norm = product object exists, and it is the Satake/Frobenius pair.** For a Hecke eigenform of weight k: a_p = α + β and p^{k−1} = αβ. This is the literal analogue of a = F + F̂, q = F·F̂. PROVED (standard).
3. **Unitarity is the wrong positivity; temperedness is the right one.** Every input that only says "this representation is unitary" (Petersson, Hermitian Hecke operators, the local unitary dual) still allows complementary-series parameters ω = ±p^σ, 0 < σ < ½. All of these break the Toeplitz Gram already at size 2. PROVED (elementary); NUMERICAL table below. Unitarity *plus* cuspidal automorphy is not enough either: Saito–Kurokawa / Howe–Piatetski-Shapiro forms are counterexamples (PROVED, literature).
4. **Over ℤ, "a unitary representation of C_ℚ whose character is the explicit formula" is equivalent to RH** (Weil + Bochner–Schwartz + GNS). It is a reformulation, not progress. PROVED (tautology); constructing it is the GAP.
5. **For ζ itself the local Satake parameters are trivially tempered** (all equal to 1). RH for ζ is not a Ramanujan statement. The archimedean pair is (ρ, 1 − ρ): sum 1, product ρ(1 − ρ) = ¼ + γ². A self-adjoint realization of ρ(1 − ρ) would be RH-equivalent. OPEN.
6. **Candidates:** Rosati (no Frobenius endomorphism over ℤ, GAP); a unitary idèle-class representation with a Frobenius-like element (equivalent to RH, OPEN); Petersson/Plancherel (unitarity only, too weak); Connes–Consani/Yoshida archimedean positivity (PROVED but local, so it cannot separate pairs with matching local data).
7. **For the Stepanov port:** the analogue of "N_k is a count" is **c(n) ≥ 0** for the log-derivative coefficients, not a(n) ≥ 0 for the Dirichlet coefficients. EPS has a(n) ≥ 0 for every n, yet c(36) = −7.17. NUMERICAL. Any ported inequality must fail on EPS and DH.

## 1. Over F_q: Hodge index = unitarity of F/√q (PROVED, standard)

Let J be the Jacobian of C/F_q, End⁰(J) = End(J) ⊗ ℚ, and ′ the Rosati involution of the canonical polarization. Weil's positivity says that the trace form ⟨φ, ψ⟩ = Tr(φ ψ′) is positive definite on End⁰(J) ⊗ ℝ, with Tr the trace on V_ℓ(J) (dimension 2g). This is equivalent to the Castelnuovo–Severi / Hodge-index inequality on C×C used in Round 15 §1.1(iii).

- F′ = V (Verschiebung) and F·F′ = q. So U := F/√q in End⁰(J) ⊗ ℝ satisfies **U·U′ = 1**: U is "unitary" for the positive involution.
- ℚ[F] is commutative and stable under ′, because F′ = q/F. A positive involution on a commutative semisimple ℝ-algebra is complex conjugation in every complex embedding, so F′ = F̄ in every embedding. Hence F F̄ = q, that is |α_i| = √q.
- **Exact dictionary:** a = Tr(F + F′) = √q · Tr(U + U^{−1}) and q = F·F′. For g = 1: a = α + ᾱ and q = αᾱ.
- **The Toeplitz Gram is the Rosati Gram of the powers of U.** T[j,k] = Σ_i ω_i^{k−j} = Tr(U^{k−j}) = Tr(U^k (U^j)′), because (U^j)′ = U^{−j}. So T is the Gram matrix of {1, U, U², …} for ⟨·,·⟩, and T ⪰ 0 at every size is Rosati positivity restricted to ℝ[U]. By Herglotz, k ↦ Tr(U^k) is a positive-definite sequence on ℤ, i.e. the Fourier coefficients of the positive measure Σ_i δ_{ω_i} on the unit circle. That is possible only if every |ω_i| = 1.

So the representation-theoretic name of Hodge index is: **Frobenius normalised by √q acts unitarily on the finite-dimensional space H¹ for the Weil-pairing/Rosati inner product.** Two features carry the weight: (i) the inner product is positive *definite*; (ii) the space is *finite-dimensional*, so for unitary U every eigenvalue has |ω| = 1. Over ℤ we have neither.

Your analysis, sharpened: "Toeplitz PSD = matrix coefficients of a unitary operator" is right. More precisely, T is the Gram of {U^k} under the Rosati trace form, so its positivity is Rosati positivity itself, not just a consequence of it.

## 2. Over ℤ: the honest equivalence (PROVED, a tautology, not progress)

Let W be Weil's functional, with spectral side Σ_ρ ĥ(γ_ρ), ρ = ½ + iγ_ρ, on test functions on ℝ₊^× ≅ ℝ. For ζ this means characters of C_ℚ trivial on Ẑ^×.

- **(a) Weil:** W(f * f^*) ≥ 0 for all f ⇔ RH.
- **(b) Bochner–Schwartz:** W positive definite ⇔ Ŵ is a positive tempered measure μ on ℝ. Here μ = Σ m_γ δ_γ, which is a measure on ℝ only if every γ is real.
- **(c) GNS:** given μ ≥ 0, take L²(ℝ, μ) with U(s) = multiplication by e^{iγs}. This is a unitary representation of ℝ (and so of C_ℚ) with cyclic vector 1. Its generator is self-adjoint with spectrum {γ}, and its distributional character is W. Conversely, if a unitary representation has character W, then W(f * f^*) = Tr(U(f) U(f)^*) ≥ 0.

So **"there is a unitary representation of C_ℚ whose character is the explicit formula" ⇔ RH.** Building it from the γ's is exactly what the rules forbid: it uses zeros as inputs. The regular representation on L²(C_ℚ) is unitary, but its spectrum is the whole unitary dual (all |·|^{it}χ), not the zeros. In Connes (1999) the critical zeros appear as an *absorption spectrum*, a cokernel or quotient, and non-critical zeros would appear as resonances. Connes states that the global trace formula, i.e. Weil positivity there, is equivalent to RH for all L-functions with Grössencharakter. Unitarity on that quotient is exactly what is missing. **GAP / OPEN.**

## 3. The key point: UNITARY ≠ TEMPERED

### 3.1 The trace = sum, norm = product object (PROVED, standard)

For a normalised Hecke eigenform f of weight k and level 1, the Euler factor (1 − a_p p^{−s} + p^{k−1−2s})^{−1} = ((1 − αp^{−s})(1 − βp^{−s}))^{−1} gives

  **a_p = α + β (trace = sum), p^{k−1} = α β (norm = product).**

This is the literal analogue of a = F + F̂, q = F·F̂. With ω = α/p^{(k−1)/2}, the normalised pair is {ω, ω^{−1}} and λ_p = a_p/p^{(k−1)/2} = ω + ω^{−1}. The Round 15 Toeplitz Gram of the pair is

  G[j,k] = ω^{k−j} + ω^{−(k−j)} = 2cos((k−j)θ_p) if ω = e^{iθ_p}; or 2cosh((k−j)σ log p) if ω = ±p^σ (sign ignored in the cosh form).

For g = 1, the Round 15 curve form m² − a_p mn + p^{k−1} n² is positive definite ⇔ a_p² < 4p^{k−1} ⇔ |ω| = 1.

### 3.2 What unitarity gives: the trivial strip (PROVED, standard)

The Petersson inner product makes the automorphic representation π unitary, so every π_p is unitary. Equivalently, the Hecke operator T_p is self-adjoint, so λ_p is real. Then λ_p = ω + ω^{−1} real allows two cases:
- (i) |ω| = 1: tempered principal series;
- (ii) ω = ±p^σ real: complementary series, unitary and generic for 0 < σ < ½.

The local unitary dual allows exactly |σ| < ½ for generic spherical representations. In particular |λ_p| < p^{½} + p^{−½}: unitarity gives only this trivial bound. This is the analogue of the trivial strip 0 < Re ρ < 1, and complementary series are the analogue of off-line zeros.

**What is shown, precisely (PROVED, elementary):** for every p and every 0 < σ < ½, the parameters {p^σ, p^{−σ}}
- keep norm = product (ω·ω^{−1} = 1);
- keep λ_p real;
- keep the local representation unitary and generic;
- keep |λ_p| below the unitary-dual bound;

yet the 2×2 Toeplitz Gram [[2, 2cosh x],[2cosh x, 2]], x = σ log p, has eigenvalue 2 − 2cosh x < 0. So **no input that only encodes local unitarity of the Satake pair can force the Toeplitz Gram.** This is a statement about the parameter space. It does not exhibit a cuspidal GL₂/ℚ form with σ > 0. None is known, and Kim–Sarnak prove σ ≤ 7/64 for GL₂ over ℚ.

**Temperedness** (|ω| = 1, the Ramanujan–Petersson conjecture) is the analogue of RH. Equivalently, π_p is weakly contained in the regular representation L²(GL₂(ℚ_p)). The Plancherel measure is supported on the tempered dual, but Plancherel positivity is a statement about the regular representation, not about which π_p occur in L²(GL₂(ℚ)\GL₂(𝔸)). So Petersson and Plancherel positivity are both "unitarity", the wrong, too-weak positivity.

### 3.3 Automorphy + unitarity + cuspidality is still not enough (PROVED, literature)

Howe and Piatetski-Shapiro (1979), and Kurokawa, found cuspidal automorphic representations of Sp(4) and U(2,1) that are non-tempered almost everywhere. The Saito–Kurokawa lift F of f ∈ S_{2k−2}(SL₂(ℤ)) has spinor L-function ζ(s−k+1) ζ(s−k+2) L(s,f). In unitary normalization this is ζ(s−½) ζ(s+½) L(s,f), with a Satake parameter p^{½}. It has an Euler product and a functional equation, and it comes from a unitary cuspidal representation. Yet it has zeros on Re s = 0 and Re s = 1. The corrected conjecture is for *generic* (Whittaker) cuspidal representations, which for GL_n is all cuspidal ones; it is OPEN in general. This sharpens the 14F criterion: "isobaric sum of unitary cuspidal **GL_n**" is doing real work. Non-GL_n unitary cuspidal data can break the line.

### 3.4 Where temperedness actually comes from (correction to my brief)

The brief said "Deligne proved it by returning to F_p and Hodge index/Weil (Kuga–Sato)." Corrected:
- **Weight 2:** Eichler–Shimura relates a_p to Frobenius on the Jacobian of X₀(N) mod p. Weil's curve RH, i.e. the Rosati / Castelnuovo positivity of §1, gives |ω| = 1. Hodge index is the right word here.
- **Weight k ≥ 3 (Δ has k = 12):** Kuga–Sato reduces it to the Weil conjectures, proved by Deligne, *La conjecture de Weil I* (1974). That proof **does not use Hodge index**; the standard conjectures remain open in characteristic p. Its positivity input is the "Rankin trick": tensor powers of a real symplectic local system have log-derivatives with non-negative coefficients, so a pole-location argument gives |α| ≤ q^{w/2}, and Poincaré duality / the symplectic pairing gives equality. It also uses Lefschetz pencils and monodromy.
- **Weight 1:** Deligne–Serre: the Galois representation has finite image, so the Satake parameters are roots of unity. Finite-dimensional plus compact gives tempered.
- **Function fields:** Drinfeld (GL₂) and L. Lafforgue (GL_n) prove Ramanujan for cuspidal GL_n via shtukas and Deligne.
- **Over ℚ, unconditionally for Maass forms:** only partial bounds, e.g. Kim–Sarnak σ ≤ 7/64, via symmetric-power functoriality plus unitarity. The amplification is the representation-theoretic Stepanov/Rankin trick: if Sym^n π were automorphic on GL_{n+1} for all n, then unitarity of Sym^n π_p (exponent nσ < ½, from the Jacquet–Shalika bound) would force σ < 1/(2n) → 0, i.e. Ramanujan (Langlands' argument; CHECKED numerically in §4F as arithmetic). Selberg's ¼ and Ramanujan for Maass forms remain OPEN.

So the honest list of positivity sources that produce temperedness is: (a) a finite-dimensional cohomology with a Frobenius and a positive form (Weil/Rosati), (b) a Rankin/tensor-power trick on an ℓ-adic system (Deligne), and (c) unitarity of *all* tensor powers (functoriality). Plain unitarity is never enough.

### 3.5 ζ itself: the archimedean Satake pair (correction to item 7 of my brief)

For ζ (GL₁), every local Satake parameter is 1, so every local component is tempered. RH for ζ is a global statement about zeros, not a Ramanujan statement. The analogue of (α, β) for the zeros is the spectral pair (ρ, 1 − ρ), as in Selberg's λ = s(1 − s):

  ρ + (1 − ρ) = 1 (sum), ρ(1 − ρ) = ¼ + γ² (product).

- **Self-adjointness gives the unitary dichotomy (PROVED).** If ρ(1 − ρ) were the eigenvalues of a self-adjoint operator (no zeros as inputs), then ¼ + γ² would be real, so γ is real (on the line) or purely imaginary (ρ real in (0,1)). The second case is the "complementary series" case.
- **ζ has no real zeros in (0,1) (PROVED, elementary).** (1 − 2^{1−σ}) ζ(σ) = Σ (−1)^{n−1} n^{−σ} > 0 (alternating, with decreasing terms), and 1 − 2^{1−σ} < 0, so ζ(σ) < 0.

So for ζ such a realization would be RH-equivalent. Over ℤ the obstruction is not "unitary vs tempered at p". It is the absence of the Hilbert space itself: the H¹ of Spec ℤ with a positive form. **OPEN.**

## 4. NUMERICAL / CHECKED (`r16ps_satake.py` → `r16ps_satake.json`)

Inputs: τ(n) from Δ = q ∏(1 − qⁿ)²⁴, exact integers, n ≤ 200; lattice counts of x² + 5y²; DH coefficients (period 5: 0, 1, κ, −κ, −1, κ = 0.284079). No zeros used.

**A. Sanity (CHECKED).**
- τ(2…13) match the known values.
- Hecke multiplicativity τ(mn) = τ(m)τ(n): 100 coprime pairs, 0 failures.
- τ(p²) = τ(p)² − p¹¹ holds for every p ≤ 14.

**B. Trace = sum, norm = product (CHECKED, exact integers).** For all 46 primes p ≤ 199, the integer discriminant τ(p)² − 4p¹¹ < 0. So the curve form m² − τ(p) mn + p¹¹ n² is positive definite, and max |λ_p| = 1.918814 (at p = 103). This is Deligne's theorem seen on data, not a proof.

**C. Δ Toeplitz Gram (NUMERICAL).** G[j,k] = 2cos((k−j)θ_p), worst minimum eigenvalue over all 46 primes:

| size K | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| min eig | 2.000 | 8.12e-02 | −9.8e-16 | −2.0e-15 | −1.5e-15 | −3.1e-15 | −3.1e-15 | −3.9e-15 |

PSD at every size, to float precision: rank 2, so the size-≥3 minima are 0. The size-2 minimum 0.0812 = 2 − max |λ_p|. Growth form: max |Σ_i ω_i^k| over k ≤ 50 and all p is 1.999998 ≤ 2.

**D. Complementary-series fake (NUMERICAL; the size-2 value 2 − 2cosh(σ log p) is exact).** ω ∈ {p^σ, p^{−σ}}, Gram 2cosh((k−j)σ log p):

| p | σ | λ_p | unitary bound p^{½}+p^{−½} | min eig K=2 | min eig K=8 |
|---|---|---|---|---|---|
| 2 | 0.01 | 2.0000 | 2.1213 | −4.80e-05 | −4.04e-03 |
| 2 | 7/64 | 2.0058 | 2.1213 | −5.75e-03 | −4.91e-01 |
| 2 | 0.25 | 2.0301 | 2.1213 | −3.01e-02 | −2.77 |
| 2 | 0.49 | 2.1165 | 2.1213 | −1.16e-01 | −13.8 |
| 199 | 0.01 | 2.0028 | 14.18 | −2.80e-03 | −2.37e-01 |
| 199 | 7/64 | 2.3447 | 14.18 | −3.45e-01 | −75.9 |
| 199 | 0.25 | 4.0221 | 14.18 | −2.02 | −1.13e+04 |
| 199 | 0.49 | 13.454 | 14.18 | −11.45 | −7.72e+07 |
| 199 | 0.50 (non-generic, Saito–Kurokawa-type) | 14.178 | = bound | −12.18 | −1.12e+08 |

Every complementary-series parameter, including tiny σ = 0.01 at p = 2, fails the Gram at size 2 and fails it more strongly at larger sizes. Growth: at p = 2, σ = 7/64, |Σ ω^k| = 2.60 at k = 10 and 44.3 at k = 50 (Round 15 §2.4 growth form). The full grid (p ∈ {2, 3, 5, 7, 11, 101, 199}) is in the JSON.

**F. Sym^n amplification (CHECKED arithmetic).** The largest n for which Sym^n {p^σ, p^{−σ}} stays in the generic unitary range nσ < ½ is 49 for σ = 0.01, 4 for σ = 7/64, and 1 for σ = 0.25 and 0.49. So "unitarity of all Sym^n" kills every σ > 0.

**G. Stepanov-falsifier sign test (NUMERICAL, n ≤ 200).** Dirichlet coefficients b(n) and log-derivative coefficients c(n), with −L′/L = Σ c(n) n^{−s}:

| model | b(n) < 0 (count) | c(n) < 0 (count) | first negative c(n) |
|---|---|---|---|
| ZETA | 0 | 0 | — (c = Λ) |
| EPS (x² + 5y²) | 0 (b counts lattice points) | 6 | c(36) = −7.17, then 54, 84, 126, 189, 196 |
| DH | 80 | 59 | c(3) = −0.312, c(4) = −1.44, … |

EPS is the instructive case. Its coefficients are a genuine count, b(n) ≥ 0 ("N_k ≥ 0" at the level of the Dirichlet series), but its log-derivative coefficients go negative. Over F_q, N_k = #C(F_{q^k}) are the coefficients of log Z. So the correct ℤ-analogue of "positivity is a count" is c(n) ≥ 0 (Euler product + Λ-positivity), not b(n) ≥ 0. Caveat: c ≥ 0 is neither necessary for RH-type statements (L(χ) has signed c) nor sufficient. It is only the ingredient Stepanov's count uses.

## 5. Candidate-by-candidate verdict

| candidate | what it gives | can it be the ℤ-polarization? | label |
|---|---|---|---|
| **Rosati-type positivity** | Exists for abelian varieties and CM motives over number fields: positive involution on End⁰, Hodge–Riemann over ℂ. Gives Frobenius bounds only after reducing mod p (§1). | Over ℤ there is no Frobenius *endomorphism* of a finite-dimensional H¹ for ζ. ζ has only H⁰ and H² of a point, and the zeros have no motive. Hodge–Riemann over ℂ fixes weights, not Frobenius. | GAP |
| **Unitary rep of C_ℚ with a Frobenius-like element** (Connes scaling action, Deninger flow) | The scaling action exists. Connes realizes the critical zeros as an absorption spectrum. | A Hilbert space where the flow is unitary *and* the character is the explicit formula is RH-equivalent (§2). Deninger's foliated-space cohomology with a positive form is a program. | CONJECTURE / OPEN |
| **Petersson / Plancherel positivity** | Unitarity of π and π_p; Plancherel support = tempered dual (regular representation only). | Too weak. Complementary series pass every unitarity constraint and fail the Toeplitz Gram at size 2 (§3.2, §4D). Saito–Kurokawa shows that even cuspidal automorphy does not help without genericity (§3.3). | PROVED (insufficient) |
| **Connes' archimedean Weil positivity** | Connes–Consani (Selecta 27, 2021): for g supported in [2^{−½}, 2^{½}] with ĝ(i/2) = ĝ(0) = 0, W_∞(g * g^*) ≥ Tr(ϑ(g) S ϑ(g)^*) ≥ 0 (Sonin-space compression plus a prolate/Toeplitz correction). This sharpens Yoshida (1992): positivity for small support. | Local: support < log 2 sees no primes, so it is the "genus" term of Round 15 (the diagonal d(g)). It cannot separate ζ from EPS or L(χ₅) from DH, because those pairs have matching local data (same gamma factors and conductors). The semilocal extension to {∞, 2, …, p} is CC's program; the global version is RH. | PROVED (local) / OPEN (global) |

**Natural object with trace = sum, norm = product:** the Satake/Frobenius pair (§3.1). Over F_q, Rosati makes U = F/√q unitary in a finite-dimensional space, which forces |ω| = 1. Over ℤ, the representation theory supplies the pair and its unitarity, but not its temperedness.

## 6. Advice for Heart's Stepanov port (Round 16)

- **What Bombieri–Stepanov uses:** (i) Riemann–Roch on C, to count auxiliary functions of bounded pole order; (ii) the q-power Frobenius as a morphism of degree q that is additive and multiplicative, so the auxiliary function f^q-type construction vanishes to order ≥ q at every F_q-point; (iii) N_k = #C(F_{q^k}) ≥ 0 is a count; (iv) the rationality of Z, to amplify over k.
- **Over ℤ:** there is no Frobenius morphism of ℤ. Auxiliary-function / polynomial methods over ℤ (large sieve, zero-density estimates, positive-proportion results) give density statements, not RH. That is consistent with Round 15 §3: any subexponential uniform bound on W(g, τ_s g) is RH-equivalent.
- **Falsifiers every ported step must pass (team rule):**
  1. **DH:** same functional equation and gamma factor as L(χ mod 5), off-line zeros. Any argument that uses only FE + local data proves a false statement for DH.
  2. **EPS (x² + 5y²):** b(n) ≥ 0 is a genuine count, yet c(36) < 0 and it has off-line zeros. Any argument whose positivity input is "the Dirichlet coefficients count something" proves a false statement for EPS. The step must use c(n) ≥ 0, i.e. the Euler product, which is what N_k ≥ 0 really is.
  3. **F_θ (non-unitary)** and the **complementary-series pair** {p^σ, p^{−σ}} (§4D): any step that uses only "λ_p real / T_p Hermitian" proves a false statement for these.
- **Where to look for the additive+multiplicative coupling:** in the GL_n world, the only known engine that turns unitarity into temperedness is the tensor-power trick (Deligne's Rankin step, Langlands' Sym^n argument), which is multiplicative (⊗) and additive (⊕) at once on representations. The ring structure of the representation ring is the closest representation-theoretic analogue of "x ↦ x^q is a ring map" that I know of. For ζ (GL₁) every tensor power is trivial, so this engine does nothing. HEURISTIC.

## 7. Bottom line

The natural object is the Satake/Frobenius pair (α, β) with a = α + β and q = αβ (or a_p = α + β, p^{k−1} = αβ). Hodge index is precisely "U = F/√q is unitary for a positive-definite form on a finite-dimensional H¹", which in automorphic language is **temperedness**, not mere unitarity. Petersson/Plancherel give only unitarity: complementary-series parameters satisfy every unitarity constraint and fail the Round 15 Toeplitz Gram at size 2 (PROVED; NUMERICAL table §4D), and even cuspidal automorphy does not rescue it without genericity (Saito–Kurokawa). Known sources of temperedness are finite-field geometry (Weil/Rosati for weight 2, Deligne's Rankin trick otherwise), finite Galois image, and functoriality of all tensor powers. For ζ itself the local parameters are trivially tempered. The relevant pair is (ρ, 1 − ρ), with product ¼ + γ², and its unitary/self-adjoint realization is RH-equivalent, as is a unitary C_ℚ-representation whose character is the explicit formula. Connes–Consani/Yoshida positivity is proved but local, so it cannot separate ζ from EPS or DH. **No representation-theoretic replacement for the polarization over ℤ is known; RH is OPEN.**

## Reproduce

`python3 r16ps_satake.py` (about 0.2 s, numpy only) → `r16ps_satake.json`.

## References (checked by web search on 2026-11-18)

- A. Connes, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*, Selecta Math. 5 (1999) 29–106.
- A. Connes, C. Consani, *Weil positivity and trace formula, the archimedean place*, Selecta Math. 27 (2021), no. 4 (arXiv:2006.13771).
- H. Yoshida, *On Hermitian forms attached to zeta functions*, in *Zeta Functions in Geometry*, Adv. Stud. Pure Math. 21 (1992) 281–325.
- P. Deligne, *La conjecture de Weil. I*, Publ. Math. IHÉS 43 (1974) 273–307.
- P. Deligne, J.-P. Serre, *Formes modulaires de poids 1*, Ann. Sci. ÉNS (4) 7 (1974) 507–530.
- H. Kim, *Functoriality for the exterior square of GL₄ and the symmetric fourth of GL₂*, J. AMS 16 (2003) 139–183; Appendix 2 with P. Sarnak (7/64 bound).
- R. Howe, I. I. Piatetski-Shapiro, *A counterexample to the "generalized Ramanujan conjecture" for (quasi-)split groups*, Proc. Sympos. Pure Math. 33, Part 1 (1979) 315–322.
- R. Schmidt, *The Saito–Kurokawa lifting and functoriality* (spinor L-function ζ(s−k+1)ζ(s−k+2)L(s,f); I verified the content but not the journal details).
- E. Bombieri, *Counting points on curves over finite fields (d'après S. A. Stepanov)*, Sém. Bourbaki 430 (1973) (as cited in ROUND15).
- Standard facts used without a specific citation: Weil's Rosati positivity; Herglotz; Bochner–Schwartz; GNS; the unitary dual of GL₂(ℚ_p); the Jacquet–Shalika bound; Drinfeld/Lafforgue.

RH / ARP / H_ψ: **OPEN**.

— Daniel Toupin / Golden Physics Project. Principal-series lane (Grok Bot).
