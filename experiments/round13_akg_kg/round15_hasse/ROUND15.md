# Round 15: the Hasse Gram, testing whether Hasse's two-integer positive form has a shadow over ℤ

> **Erratum (Golden Eye AUDIT15, team/golden-eye 798f3ac), 2026-11-18.**
> 1. The one-shift block follows from W ≥ 0 by Cauchy–Schwarz (W is translation-invariant), so it is a **tautological control**, strictly weaker than W ≥ 0. EPS is indefinite from width 3.0 but fails the block only at s = 4.40 (width 5.0).
> 2. The pole term is 2cosh(s/2)·ĝ(i/2)ĝ(−i/2), which equals |ĝ(i/2)|² only for even g. It is an algebraic identity, so its match is a **code check**.
> 3. CC 2021 needs supp g ⊂ [2^{−1/2}, 2^{1/2}] and ĝ(i/2) = 0.
> 4. "Kronecker ⇒ ρ → 1" is **HEURISTIC**.
> 5. Toeplitz = 2·deg is **algebra-checked only**.
> 6. The scripts now run from any working directory (`_paths.py`).
> Numbers are unchanged, and 'all s ⇔ RH' (PROVED, sketch) stands per AUDIT15.

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Daniel (t417u/t418u): "dig into this, do your thinking magic", about CLUES.md clue 1.

**Labels:** PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN.
**Status:** RH, ARP and H_ψ are all **OPEN**. Nothing here proves RH. No zeta zeros are used as inputs: every matrix is built from integer data c(n), each model's own gamma factor, and its own poles. Zeros appear only in the explanatory explicit-formula identities of §3.

Models (14B/14F):
- ZETA;
- EPS (x² + 5y²);
- DH (Davenport–Heilbronn);
- BEUR2 (spliced Beurling, i.e. BEUR″);
- FTH, the non-unitary control F_θ with θ = 0.2;
- FTAU, the unitary control F_{iτ} with τ = 3.

## Summary

1. **Why the curve form is positive (Part 1).** The positivity of deg(m − nF) = m² − a·mn + q·n² is *inherited*, not computed:
   - degree is a count;
   - it is a quadratic form because the dual isogeny is additive (polarization);
   - in genus g, Hodge index on C×C makes the Toeplitz "Hasse Gram" T[j,k] = Σ ω_i^{k−j} positive.

   a = F + F̂ is an addition and q = F F̂ is a multiplication of the *same* algebraic number. The curve-RH proof that avoids surfaces (Bombieri–Stepanov) runs on one fact: x ↦ x^q is additive and multiplicative **at the same time**. That is the sharpest mathematical form of Daniel's coupling, and it is exactly what ℤ lacks in characteristic 0.
2. **The one-shift Hasse block [g, τ_s g] separates ζ from every fake on our data (Part 2), but only as a tautological control** (AUDIT15). W is translation-invariant, so the block is Cauchy–Schwarz for W. It is implied by W ≥ 0, every failure is an explicit ghost g ± τ_s g of width ≤ l + s, and it is strictly weaker than W ≥ 0: EPS is indefinite from width 3.0, yet fails the block only at s = 4.40 (width 5.0). At l = 0.6 it holds for ζ and F_{iτ} at all 499 shifts s ≤ 8, and fails for F_θ (from p = 3), BEUR″ (p = 37), EPS (p = 149) and DH (p = 547). NUMERICAL. ζ's pass is a code check; F_{iτ}'s pass is expected. ζ nearly saturates it: max ρ = 0.9996.
3. **It is not a new criterion (Part 3).** The family over all s is RH-equivalent (PROVED, sketch), and so is any weak, subexponential version (the Stepanov amplification). The full Hasse Gram reproduces 14B's ghost onsets and adds nothing new.
4. **The GAP** is a zero-free, Stepanov-type uniform bound on W(g, τ_s g). It would need a Frobenius that is additive and multiplicative at once. Such a Frobenius exists only in characteristic 1 (Connes–Consani's arithmetic site). OPEN.

## Part 1 — Why the curve form is positive, and what it maps to over ℤ

### 1.1 Over F_q: three ingredients, each doing a different job

Take a curve C/F_q of genus g, Frobenius F, and α₁…α_{2g} the eigenvalues of F on H¹, so that #C(F_{q^k}) = N_k = q^k + 1 − a_k with a_k = Σ_i α_i^k.

**(i) Positivity is a count (geometric, elementary).** For an elliptic curve, deg(φ) of a non-zero isogeny φ is the degree of a finite map of curves, a positive integer: for separable φ it is #ker φ. So deg(m − nF) ≥ 0 for integers m, n. This part is just counting, which is why it is "physical": a degree is a number of points.

**(ii) Quadraticity is the additive structure (geometric, deep).** deg is a *quadratic form* on End(E) because deg φ = φ∘φ̂ (φ̂ = dual isogeny, the Rosati adjoint) and φ ↦ φ̂ is *additive*. The additivity of the dual is the real theorem: it comes from Pic⁰ / the Weil pairing, i.e. from the polarization of the Jacobian (the theta divisor). Then

  deg(m − nF) = (m − nF)(m − nF̂) = m² − (F + F̂) mn + (F F̂) n² = m² − a·mn + q·n²,

with **F + F̂ = a (an addition)** and **F F̂ = q (a multiplication)**. Positivity of the form ⇔ a² − 4q ≤ 0 ⇔ a and q are the sum and product of a complex-conjugate pair α, ᾱ ⇔ |α| = √q. That is Hasse.

**(iii) For general g, the surface C×C (geometric, deep).** Let H = C×pt, V = pt×C, Γ_j = graph of F^j. Then H² = V² = 0, H·V = 1, Γ_j·H = q^j, Γ_j·V = 1, Γ_j² = q^j(2 − 2g), and Γ_j·Γ_k = q^j N_{k−j} for j ≤ k. With the primitive part Γ′_j = Γ_j − H − q^j V (the part orthogonal to the ample class H + V), a short computation gives (PROVED, standard):

  Γ′_j · Γ′_k = −q^j · a_{k−j},  Γ′_j² = −2g q^j.

The **Hodge index theorem** (NS(C×C)_ℝ has signature (1, ρ−1) because H + V is ample, proved from Riemann–Roch on the surface) says the primitive part is negative definite. Normalising by q^{(j+k)/2} and flipping the sign gives the **finite-field Hasse Gram**

  T[j,k] = −Γ′_j·Γ′_k / q^{(j+k)/2} = Σ_i ω_i^{k−j},  ω_i = α_i/√q,  T ⪰ 0 for every size.

For g = 1 the 2×2 block is [[2, a/√q], [a/√q, 2]] and (m, −n√q) T (m, −n√q)ᵀ = 2(m² − a mn + q n²) = 2 deg(m − nF) (an algebraic identity: CHECKED-as-algebra only, not a test against isogeny degrees). The positive Toeplitz form forces |ω_i| = 1: if a positive metric is preserved by F/√q, F/√q is unitary. **That is literally 14B's holonomy lemma and Daniel's No-Ghosts theorem.** In the curve case the theorem is not the hard part; the hard part is having a positive metric that Frobenius preserves.

**Where add and mult meet.** Correspondences on C×C form a *ring*: + is the sum of divisors, × is composition, and the intersection with Δ is the trace (Lefschetz). Frobenius obeys F∘Fᵗ = q·Δ (F times its transpose is multiplication by q). And in characteristic p, x ↦ x^q is **at the same time multiplicative and additive**: (x + y)^q = x^q + y^q. That single fact is what Bombieri–Stepanov use. Their auxiliary-polynomial proof needs only (a) Riemann–Roch on C, to count auxiliary functions, and (b) additivity of x ↦ x^q, to make the auxiliary function vanish to high order at every rational point. It gives |N_k − q^k − 1| ≤ C q^{k/2} for large k. Then amplification over k (the generating function Σ a_k T^k = Σ_i α_i T/(1 − α_i T) is rational) upgrades it to |α_i| = √q.

### 1.2 Dictionary F_q → ℤ (with what is actually known)

| curve over F_q | over ℤ (log scale u = log x) | status over ℤ |
|---|---|---|
| C, and C×C | adele class space A_ℚ/ℚ^× (Connes 1999); arithmetic site (ℕ^×-topos, ℤ_max) and its **square** with structure sheaf ℤ_min ⊗_B ℤ_min (Newton polygons) (Connes–Consani, Adv. Math. 2016); scaling site, with periodic orbits = circles of length log p (CC 2016–17) | object exists, no intersection theory |
| Frobenius F, F^k | scaling flow x ↦ e^u x; closed orbits of length log p^k. CC: Frobenius correspondences Ψ(λ), λ ∈ ℝ₊^×, with Ψ(λ)∘Ψ(λ′) = Ψ(λλ′) (generic λλ′) | exists (char 1) |
| F is a ring endomorphism (x ↦ x^q additive) | on the arithmetic site, Fr_λ(x) = x^λ is a semifield automorphism of ℝ_max **only in characteristic 1**; on the adeles, scaling is additive but not multiplicative; Borger's λ-rings: Frobenius lifts ψ_p are ring maps, but trivial on ℤ | **missing in char 0** |
| divisor mΔ − nΓ_F | test function f = m·g + n·τ_{log p} g (bump g at the identity plus its shift by log p) | this round |
| intersection D·E | −W(f, h), the Weil bilinear form; explicit formula = Lefschetz trace (Connes trace formula) | PROVED |
| subtracting H, V (primitive part) | the pole terms ĝ(0), ĝ(1). The pole part of the shift block is 2cosh(s/2)·ĝ(i/2)ĝ(−i/2) (= \|ĝ(i/2)\|² only for even g; AUDIT15), the ℤ image of (q + 1)/√q = √q + 1/√q. It is an algebraic identity (c_±(τ_s g) = e^{±s/2} c_±(g)), so the match is a code check | PROVED; code check |
| N_k = #C(F_{q^k}) | the prime sum Σ_{\|log n − s\|<l} Λ(n) n^{−1/2} X(log n − s): ψ in a window at e^s | — |
| genus 2g = −Γ′² | W(g) for supp g of width < log 2 (no primes): the archimedean "genus" of g. It diverges as g → δ, so the genus is infinite (CC 2026, *On the Jacobian of Spec ℤ*: Jacobian = orbit space of the scaling flow, infinite genus) | — |
| local positivity of the genus term | W(g * g*) ≥ 0 for supp g ⊂ [2^{−1/2}, 2^{1/2}] (so supp g * g* ⊂ [1/2, 2]) **and ĝ(i/2) = 0**, so it does not cover the diagonal's pole part: Connes–Consani, *Weil positivity and trace formula, the archimedean place* (Selecta 2021): a Sonin-space trace plus a prolate correction controlled by Toeplitz matrices | PROVED (literature) |
| Riemann–Roch on C | Riemann–Roch for Spec ℤ̄ (CC, Bull. Sci. Math. 2023; for the ring ℤ over the sphere spectrum, CRAS 2024) | partial |
| Hodge index on C×C / Rosati positivity | Weil positivity W(f * f^*) ≥ 0, f^*(x) = x^{−1} conj f(1/x) (Daniel's x → 1/x), equivalent to RH (Weil 1952) | **OPEN = RH** |
| polarization (theta divisor of Jac C) | ? The Jacobi theta Σ e^{−πn²x}, Poisson on ℤ, is the obvious candidate, but nothing makes it a polarization of the CC Jacobian | **missing** |

### 1.3 The missing ingredient, stated exactly

Over F_q the positivity is never computed from ζ. It is **inherited**:
1. from a finite-dimensional H¹ (dim 2g);
2. carrying a positive form, via the polarization / ampleness of H + V;
3. which Frobenius preserves up to the scalar q, because F is a *ring* endomorphism and F Fᵗ = q holds as an identity of correspondences.

Over ℤ, all three fail:
1. The genus is infinite. Measured here: the diagonal W(g) grows as the window shrinks (λ_min(W_l) = 0.0089 at l = 0.6, 0.23 at l = 0.3), and §2.4 shows that already in genus 2 the k = 1 block does not suffice.
2. No polarization is known.
3. No characteristic-0 Frobenius is additive and multiplicative at once.

**Could Daniel's add/mult "same time" coupling be the missing ingredient?** Yes, in one precise sense, which is the Stepanov route rather than the Hodge route. The curve-RH proof that uses no surface at all (Bombieri–Stepanov) needs exactly two things:
- a Riemann–Roch count, which partly exists for Spec ℤ̄ (CC 2023);
- a Frobenius that is **simultaneously additive and multiplicative**, of "degree" e^u.

That second item is Daniel's coupling, made into a precise mathematical statement. Over ℤ the scaling flow is additive only (λ(x + y) = λx + λy, but λ(xy) ≠ λx·λy). The semiring version exists only in characteristic 1 (Connes–Consani's Fr_λ on ℝ_max). In a quick search I found no published Stepanov-type argument on the arithmetic site. Label: **HEURISTIC / OPEN**. It is a lead, not a result.

---

## Part 2 — The Hasse Gram for ℤ

### 2.1 Construction (`hasse_core.py`)

**Basis and shifts.**
- Basis: Dirichlet sines φ_i on a base window of width l < log 2. With that width no prime sees g alone, so the diagonal is pure "genus" (arch + pole).
- Shift: τ_s φ(u) = φ(u − s), the ℤ-analogue of Γ_{F^k}, with s = k log p.

**The bilinear block.** B_s[i,j] = W(φ_i, τ_s φ_j). It comes from the 14B linear functional Q applied to the symmetrized cross-correlation (X_ij(t − s) + X_ij(−t − s))/2, where X_ij is in closed form (14B `Aclosed`).
- Translation invariance gives G[(a,i),(b,j)] = B_{s_b − s_a}[i,j] and B_{−s} = B_sᵀ.
- Inputs are integer and own-gamma-factor data only. **No zeros are used.**

**Checks.**
- B_0 equals the 14B Weil matrix to 6e-14, 1.4e-13, 7e-14 and 1.3e-13 for ZETA, EPS, DH and FTH (`check_core.json`).
- An independent brute-force grid quadrature of W(g + τ_s h) agrees to 1e-7 to 1e-6 relative, at s = 0.45, log 2, log 3 and 1.7, for ZETA, EPS, FTH and FTAU (`check_shift.json`). That is the grid's own accuracy.
- Doubling the panels at the tightest shifts changes μ by less than 7e-14 (`check_quad.json`).

**The add+mult dictionary, measured.** For g = first sine at l = 0.6 (ζ), the pole part of b_s, divided by 2cosh(s/2), equals 0.488418 = ĝ(i/2)ĝ(−i/2) at every s from 0.6 to 7.8. That is an algebraic identity (the ℤ image of (q + 1)/√q), so this is a code check, not a measurement. The prime-window sum tracks it (17.691 vs 17.711 at s = 7.18), and the remainder b_s stays within ±0.008, while d = W(g) = 0.0118. This is N_k ≈ q^k + 1 with the error bounded by the genus, played out on ℤ (CHECKED).

### 2.2 (a)+(c) One-shift Hasse block [g, τ_s g] (`hasse_blocks.py`)

For a 2×2 block "PSD for every g in the window", the condition is W_l ± Bsym_s ⪰ 0. Write
- μ = min over the two signs of λ_min(W_l ± Bsym_s), and
- the Hasse ratio ρ(s) = ‖W_l^{−1/2} Bsym_s W_l^{−1/2}‖.

The Hasse inequality |b_s(g)| ≤ d(g) for all g is the same as ρ ≤ 1. The scan covers s ∈ [l, s_max] at a 0.02 step plus every prime shift log p (all of them up to p ≈ 300, then thinned), with s_max = log(n_max) − l − 0.05 ≈ 8.

| model | λ_min(W_l) l=0.6 | l=0.6: violating shifts / scanned | first violating s (first prime p) | min μ | max ρ | l=0.3: violations, first s | l=0.3 max ρ |
|---|---|---|---|---|---|---|---|
| ZETA | 0.00892 | 0/499 | none | 2.12e-05 | 0.9996 | none | 0.9560 |
| FTAU | 0.0737 | 0/499 | none | 6.24e-04 | 0.9963 | none | 0.8685 |
| EPS | 2.21 | 220/498 | 4.40 (p = 149) | -1.69e+01 | 4.2148 | 146, 5.64 | 2.5805 |
| DH | 0.786 | 75/499 | 6.30 (p = 547) | -2.11e+00 | 1.4992 | none | 0.9993 |
| BEUR2 | 0.00892 | 118/384 | 3.60 (p = 37) | -1.14e-01 | 9.6598 | none | 0.9560 |
| FTH | 0.0176 | 377/499 | 0.66 (p = 3) | -4.49e+00 | 2.5384 | 205, 1.78 | 2.4975 |

**Readout.**
- At **l = 0.6** the one-shift block separates ζ from **every** fake. This is a tautological control, a weaker ghost detector than 14B's W ≥ 0 (Cauchy–Schwarz; AUDIT15 §1): each failure exhibits a ghost g ± τ_s g of width ≤ l + s, and EPS (indefinite from width 3.0) first fails at width 5.0. ζ and the unitary control F_{iτ} pass at every one of 499 shifts. Each fake fails beyond a model-dependent onset: F_θ at s = 0.66 (p = 3), BEUR″ at 3.60 (p = 37, the first window that reaches its Beurling primes ≥ 55), EPS at 4.40 (p = 149), DH at 6.30 (p = 547). NUMERICAL.
- At **l = 0.3** (narrower g, larger "genus" d) only EPS (from s = 5.64) and F_θ (from 1.78) fail before the data cutoff. DH (max ρ 0.9993) and BEUR″ (0.956) pass. **The onset moves out as the genus grows.** That is the expected e^{δs}-growth picture: violation needs |G(ρ)|² e^{δs} ≳ d.
- **ζ saturates the bound:** max ρ = 0.99959 at s = 5.78 (l = 0.6). It is the ℤ-shadow of the Weil bound's sharpness: under RH, ρ(s) = sup_g |Σ_γ |G(γ)|² e^{−iγs}| / Σ_γ |G(γ)|², and Kronecker makes it approach 1 (**HEURISTIC**; F_{iτ} nearly saturates too, at 0.9963). The analogy to s = Δ / de Bruijn–Newman Λ = 0 is also HEURISTIC. The ρ values are NUMERICAL.

### 2.3 (b) Full Hasse Gram (`hasse_gram.py`)

The generalized eigenproblem G v = λ N v, with N-null directions dropped and ghosts counted below −1e-7·max|G| (the 14B convention). Families:
- P_p = {k log p}, block-Toeplitz: the literal ℤ[F_p];
- INT = {log n, n ≤ K}: the multiplicative monoid ℕ;
- PP = {0} ∪ {log p^k ≤ log K}.

| model | P2: width 8.2 λ_min / ghosts / onset width | P3: 8.3 | P5: 7.0 | INT (K=40, width 4.29) | PP (K=128, width 5.45) | 14B ghost onset L |
|---|---|---|---|---|---|---|
| ZETA | -7.6e-14 / 0 / — | -5.1e-14 / 0 / — | -9.1e-15 / 0 / — | -7.8e-15 / 0 / — | -1.6e-14 / 0 / — | none to 8 |
| FTAU | -3.4e-10 / 0 / — | -3.4e-10 / 0 / — | 7.3e-07 / 0 / — | -3.3e-10 / 0 / — | -3.4e-10 / 0 / — | — (14F: 0 ghosts L=2,4) |
| EPS | -2.8e+01 / 24 / 3.37 | -2.1e+01 / 18 / 3.90 | -7.3e+00 / 9 / 5.43 | -2.1e+00 / 4 / 3.08 | -5.8e+00 / 16 / 3.16 | 3.0 |
| DH | -1.0e+01 / 8 / 4.07 | -6.2e+00 / 6 / 4.99 | -2.8e+00 / 2 / 5.43 | -1.0e+00 / 4 / 3.64 | -2.5e+00 / 8 / 3.74 | 4.0 |
| BEUR2 | -5.2e+00 / 68 / 4.76 | -3.4e+00 / 35 / 4.99 | -1.7e+00 / 19 / 5.43 | -5.1e-01 / 8 / 4.10 | -1.6e+00 / 43 / 4.21 | 4.5 |
| FTH | -6.9e+00 / 59 / 1.29 | -6.0e+00 / 52 / 1.70 | -2.4e+00 / 33 / 2.21 | -1.3e+00 / 53 / 1.29 | -2.6e+00 / 76 / 1.29 | 1.0 |

**Readout.**
- **ζ and F_{iτ} have no ghosts in any family**, at every width up to 8.3. ζ's λ_min sits at the float floor, −1e-14 to −8e-14, as in 14B. F_{iτ} sits at −3.4e-10, its own quadrature floor, and is well inside threshold.
- **Every fake grows ghosts**, and the onset widths reproduce 14B's ghost onsets: EPS ≈ 3.1–3.4 (14B: L = 3), DH ≈ 3.6–4.1 (L = 4), BEUR″ ≈ 4.1–4.8 (L = 4.5), F_θ ≈ 1.3 (L = 1).
- So the full Hasse Gram is simply the Weil form on a lattice of translates. It finds nothing 14B didn't. Its ghost counts can exceed 14B's at the same width (e.g. EPS P2 at 8.2: 24 vs 10) because the shifted windows carry kinks, i.e. higher frequencies than 14B's basis. Both are lower bounds on the inertia of W on width-L functions (the basis dependence Golden Eye flagged in AUDIT14).
- **The 2×2 one-shift block is much weaker but still decisive:** it needs s ≈ 4.4 (EPS) and 6.3 (DH) where the Gram needs width ≈ 3.1 and 3.6. This is the Hasse/Stepanov pattern: a weak, one-Frobenius bound, amplified over the powers k.

### 2.4 (d) Finite-field side by side (`ff_hasse.py` → `ff_hasse.json`)

| object | 2×2 Hasse block [[2g, c₁],[c₁, 2g]] | Toeplitz T up to K = 8–10 | note |
|---|---|---|---|
| elliptic curves, probe B (p = 5…19) | PSD (e.g. p = 5: eigenvalues 0.658, 3.342) | rank 2, 0 ghosts | (m, −n√q)T(…)ᵀ = 2 deg(m − nF) for \|m\|, \|n\| ≤ 4 (algebra check only: compared to m² − a·mn + q·n², not to isogeny degrees); N_k > 0 |
| class-0 partial zeta Z₀ (probe B "Epstein", no Euler product) | c₁ = √q > 2: **1 ghost** for every q ≥ 5 (p = 5: −0.236) | 1 ghost | the EPS failure exactly |
| genus-2 curves y² = f(x), p = 5, 7, 11, 13 (brute-force N₁, N₂) | PSD | rank 4, 0 ghosts, all \|ω\| = 1 | Weil |
| synthetic genus-2 "fake" {e^{±iφ}, r, 1/r} (FE, off-circle pair) | **passes** (c₁ = 1.32, 0.45, 3.11) | first ghost at K = 3 for all three, i.e. size 2g = 4 | the one-shift block [Δ, Γ_{F^k}] first fails at k = 3, 4, 6 |

Lesson (PROVED, elementary): in genus g, the Toeplitz block of size 2g decides RH (Carathéodory–Toeplitz). The single block at k = 1 decides it only when g = 1. A family of one-shift blocks over all k decides it by growth: |Σ ω^k| ≤ 2g for all k ⇔ |ω| = 1. Over ℤ the genus is infinite, so no finite Gram suffices. Only the growth form survives, and §2.2 sees exactly that.

---

## Part 3 — Verdict

**1. A ζ-vs-all-fakes separator exists on our data: the one-shift Hasse block.**

  d(g) := W(g) ≥ |b_s(g)| := |W(g, τ_s g)| for every g supported in a window of width l < log 2, for every shift s.

At l = 0.6 it holds for ζ and F_{iτ} on all 499 shifts up to s ≈ 8 (p ≲ 2800), and it fails for EPS, DH, BEUR″ and F_θ (NUMERICAL). It is implied by W ≥ 0 (Cauchy–Schwarz) and strictly weaker, so it is a tautological control (AUDIT15). Its shape is genuinely Hasse's |a| ≤ 2g√q:
- the **diagonal is prime-free**: archimedean + pole, the "genus", positive by Connes–Consani 2021 only for supp g ⊂ [2^{−1/2}, 2^{1/2}] with ĝ(i/2) = 0, so the pole part is not covered;
- the **off-diagonal is one Frobenius window**: the pole main term 2cosh(s/2)·ĝ(i/2)ĝ(−i/2), the analogue of q^k + 1, minus the prime count near e^s, the analogue of N_k.

So it reads: *prime-count fluctuation at e^s, tested by g ⊗ g, is bounded by the archimedean genus of g.*

**2. But it is tautological as a criterion: RH-equivalent.** PROVED (sketch, standard). By the explicit formula, b_s(g) = Σ_ρ G(γ_ρ) G(−γ_ρ) e^{−iγ_ρ s}, with G(z) = ∫ g e^{izu} and ρ = ½ + iγ_ρ.
- **RH ⇒ the bound.** Under RH, b_s = Σ|G(γ)|² e^{−iγs}, so |b_s| ≤ Σ|G|² = d. This is Cauchy–Schwarz / Bochner.
- **An off-line zero ⇒ failure.** A zero with Re ρ = ½ + δ, δ > 0, contributes e^{δs}·G(γ)G(−γ). That coefficient is non-zero for generic g, since G is entire and non-zero. The Laplace transform ∫₀^∞ b_s e^{−zs} ds would then have a pole at Re z = δ > 0, impossible if b_s is bounded. So boundedness fails.
- **The weak version is no easier.** The same argument shows that **any** bound |b_s(g)| ≤ C(g) e^{εs} for all ε > 0 already implies RH. That is the ℤ-version of Stepanov's amplification. Over ℤ such a weak bound is equivalent to the smoothed ψ(x) = x + O(x^{1/2+ε}), i.e. RH itself.
- **Finite checks prove nothing.** Any finite range of s is only a necessary condition. A zero with tiny δ or great height evades it, as DH and BEUR″ do at l = 0.3.

**3. The theorem that would turn it into a proof (the GAP).** A uniform-in-s bound |W(g, τ_s g)| ≤ C(g)·e^{o(s)} for one fixed g of width < log 2 with G(γ) ≠ 0 off the line, proved **without zeros**. Over F_q this is exactly what Bombieri–Stepanov prove, from (i) Riemann–Roch on C and (ii) additivity of x ↦ x^q. Over ℤ:
- (i) partly exists (CC Riemann–Roch for Spec ℤ̄);
- (ii) exists only in characteristic 1 (Fr_λ on ℤ_max / ℝ_max).

The GAP is a Stepanov-type auxiliary-function count on the arithmetic or scaling site. **OPEN.**

**4. Where the analogy breaks** (each item measured or proved above):
- (i) **Infinite genus.** In genus g a 2g-block decides (fake genus-2: first ghost at K = 3 every time). Over ℤ the diagonal "genus" W(g) depends on the window, and no finite Gram decides.
- (ii) **Incommensurate clocks.** Over F_q all orbit lengths are multiples of log q, so one F generates everything and Γ_F counts all degree-1 points at once. Over ℤ the log p are ℚ-independent, so there is no single generator. The P_p Grams see only one prime's powers, while ghosts are global (14D/14F).
- (iii) **No Hodge index.** Over F_q the Gram is positive because it is an intersection form on a surface with an ample class. Over ℤ the only route to the Gram's positivity is RH itself. The Hasse Gram (b) is just W on a lattice of translates and reproduces 14B's ghost onsets exactly.
- (iv) **No additive Frobenius in characteristic 0.** This is the precise place Daniel's add/mult coupling would have to enter.

**Labels.** The F_q facts are PROVED. The engine and dictionary are CHECKED. The ζ/fake block and Gram results are NUMERICAL. The Hasse-block ⇔ RH equivalence and the weak-bound amplification are PROVED (sketch). The Stepanov-on-the-arithmetic-site route is HEURISTIC. RH, ARP and H_ψ are **OPEN**.

## Best next step

**Round 16: Stepanov for ℤ, finite version.** Bombieri's proof of curve RH is a dimension count: Riemann–Roch gives enough functions of bounded pole order, and Frobenius additivity makes them vanish to order ≥ q at the rational points, so #points ≤ (degree)/q. Build the literal ℤ-toy:
- the space of trigonometric/Gaussian test functions on the window, as the "functions of bounded degree";
- the CC tropical Frobenius Fr_λ (x ↦ λx on max-plus), which is additive and multiplicative at once in characteristic 1;
- then ask whether a Stepanov-type count bounds the prime window sum in b_s by C·e^{s/2}·(genus) with **any** constant. By §3.2 any such bound with a subexponential factor would imply RH, so the test must be done in the finite toy first: F_q, probe B. There the Stepanov count is checkable line by line, and we can see exactly which step uses additivity. Then port that step.

## Files (this folder)

| file | role |
|---|---|
| `hasse_core.py` | shifted bilinear Weil block B_s, models (+ FTAU), translation-invariant Gram blocks |
| `check_core.py` → `check_core.json` | B_0 = 14B W (≤ 1.4e-13) |
| `check_shift.py` → `check_shift.json` | independent grid quadrature of W(g + τ_s h) |
| `check_quad.py` → `check_quad.json` | panel convergence at the tightest shifts |
| `hasse_blocks.py` → `blocks_<M>_l{0.6,0.3}.json`, `blocks*.log` | parts (a) and (c) |
| `hasse_gram.py` → `gram_<M>_<FAM>_l0.6.json`, `gram*.log` | part (b) |
| `ff_hasse.py` → `ff_hasse.json` | part (d) |
| `summarize15.py` → `round15.json`, `summary15.md` | tables |

Reproduce (from any working directory: every script imports `_paths.py`, which sets sys.path and chdirs to this folder; for a scratch copy elsewhere set `GPP_R13=/path/to/round13_akg_kg`):
1. `python3 check_core.py; python3 check_shift.py; python3 check_quad.py`
2. `for m in ZETA EPS DH BEUR2 FTH FTAU; do python3 hasse_blocks.py $m 0.6; python3 hasse_blocks.py $m 0.3; for f in P2 P3 P5; do python3 hasse_gram.py $m $f 14; done; python3 hasse_gram.py $m INT 40; python3 hasse_gram.py $m PP 128; done`
3. `python3 ff_hasse.py; python3 summarize15.py`

Data cutoff: NX = e^{8.6} (14B `models14b.NX`). Bases: fmax 45, and 125 for DH.

## References

- A. Weil, *Sur les "formules explicites" de la théorie des nombres premiers* (1952).
- E. Bombieri, *Counting points on curves over finite fields (d'après S. A. Stepanov)*, Sém. Bourbaki 430 (1973).
- E. Bombieri, *Remarks on Weil's quadratic functional in the theory of prime numbers* (2000).
- A. Connes, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*, Selecta Math. 5 (1999).
- A. Connes and C. Consani:
  - *Geometry of the arithmetic site*, Adv. Math. 291 (2016);
  - *The scaling site*, CRAS 354 (2016), and *Geometry of the scaling site*, Selecta Math. 23 (2017);
  - *Weil positivity and trace formula, the archimedean place*, Selecta Math. 27 (2021);
  - *Riemann–Roch for Spec ℤ̄*, Bull. Sci. Math. 187 (2023), and *Riemann–Roch for the ring ℤ*, CRAS 362 (2024);
  - *Spectral triples and ζ-cycles*, Enseign. Math. 69 (2023);
  - *On the Jacobian of Spec ℤ*, J. Noncommut. Geom. (forthcoming, 2026).
- J. Borger, *Λ-rings and the field with one element* (2009).
- C. Deninger, *Analogies between analysis on foliated spaces and arithmetic geometry* (2008).
- S. Haran, *Non-additive geometry*, Compositio Math. (2007).

RH / ARP / H_ψ: **OPEN**.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
