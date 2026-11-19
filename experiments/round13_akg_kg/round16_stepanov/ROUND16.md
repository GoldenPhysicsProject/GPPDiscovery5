# Round 16: Stepanov for ℤ, finite version

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Daniel (t419u): "Yes and anything else that you think you should", about Round 15's next step.

**Labels:** PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN.
**Status:** RH, ARP and H_ψ are **OPEN**. Nothing here proves RH. No zeta zeros are proof inputs. Every ℤ-side number comes from integer data c(n), each model's own gamma factor and poles (Round 14B/15 engine), or from mpmath values of the model's own Dirichlet series on the real axis. The F_q side is exact finite-field arithmetic.

Models: ZETA; EPS (x² + 5y²); DH (Davenport–Heilbronn); BEUR2 (spliced Beurling, BEUR″); FTH (non-unitary control F_θ, θ = 0.2); FTAU (unitary control F_{iτ}, τ = 3).

## Summary

1. **Bombieri–Stepanov executed line by line (Part 1, CHECKED).** On 4 probe-B elliptic curves and 2 genus-2 curves over F_q (q = p², p = 5…13) I built the auxiliary function F = Σ c_uv u^p v^q explicitly. It is non-zero. It is an exact p-th power and vanishes at every rational point. It gives N ≤ 1 + pole(F)/p (e.g. q = 25: N = 27 ≤ 39; Bombieri's formula gives 41). The twist gives the lower bound, and the extensions plus the functional equation give |ω| = 1.
2. **The add+mult step is S4+S5.** Multiplicativity of x ↦ x^q, together with the fact that rational points are its fixed points, makes v^q equal to v at *every* rational point at once (S4). Additivity of x ↦ x^p turns the whole linear combination into one perfect p-th power, so it vanishes to order p = √q everywhere at once (S5). Ablations (CHECKED): without additivity the bound collapses to 1779 against q = 121; with a Frobenius that does not fix the points, F vanishes at only 13 of 129 points.
3. **The port map (Part 2).** Riemann–Roch, Siegel's lemma, zeros = poles (the product formula), the functional equation, and even the lower-bound step (via Landau's theorem plus ζ(σ) < 0 on (0,1)) all have real ℤ-analogues. The **additive-and-multiplicative Frobenius (S5) has none in characteristic 0**. It exists mod p, one prime at a time (Fermat), and globally only in characteristic 1 (Connes–Consani Fr_λ). The extension step (S8) has only a formal analogue: τ_{ns} = τ_s^n holds exactly, but there is no constant-field extension, so no uniform constant.
4. **Best provable θ (Part 3).** θ = 1 in ψ-normalisation (θ_W = ½ in Hasse-block normalisation), i.e. |b_s − pole| ≤ C·e^{s/2}. Chebyshev's binomial auxiliary integer proves it with C explicit (PROVED, zero-free). Diamond–Erdős improve the constant to 1 + ε, and Vinogradov–Korobov give a subexponential *saving* exp(−c s^{3/5}(log s)^{−1/5}) using a zero-free region. Any θ < 1 would be a zero-free strip Re s > θ, which is OPEN. RH needs θ = ½ (θ_W = 0).
5. **What fails on the fakes (NUMERICAL).** Chebyshev's add+mult identity Σ c(n)⌊x/n⌋ = log⌊x⌋! is exact only for ζ. Its *inequality* consequences (θ = 1) also hold for EPS, DH and BEUR″. Only FTH breaks the Chebyshev bound, and only EPS and FTH break the large-sieve bound. A genuine Stepanov argument must use the exact identity, not its θ = 1 shadow.
6. **The tensor-power analogue (Part 3).** For every fake, the Hasse bound holds at s₀ and then fails at n·s₀ (EPS from n = 4–8, DH 5–11, BEUR″ 3–6, FTH 2–5). So a bound at s does **not** propagate to ns without structure. Over F_q the structure is finite genus (Hankel rank 2g). Over ℤ the Hankel rank is essentially full (177 of 184 for ζ).
7. **The GAP** is one lemma: an integer Stepanov lemma, i.e. auxiliary integers with square-root-sharp size against divisibility (§4.2). It is RH-strength. The ℤ object it needs is exactly the missing S5.

---

## Part 1 — Bombieri's proof executed on concrete curves (`ffq.py`, `stepanov_ff.py`, `ablation_ff.py`)

### 1.1 Setting

- C: y² = f(x), deg f = 2g + 1, coefficients in F_p. Base-changed to F_q with q = p². Bombieri needs q to be a square and q > (g + 1)⁴; the base change is itself step (iv).
- O is the single point at infinity, with ord_O x = −2 and ord_O y = −(2g + 1).
- Curves:
  - elliptic (probe B): y² = x³ + x + 1 (p = 5, 11), x³ + x + 3 (p = 7), x³ + x + 6 (p = 13);
  - genus 2 (Round 15): x⁵ + 2x³ + x + 1 (p = 11), x⁵ + x⁴ + 5 (p = 13).
- Each curve is paired with its quadratic twist y² = d·f(x), d a non-square in F_q.

### 1.2 The steps, with ingredients and checks

Ingredient tags: **(i)** additivity of x ↦ x^Q, **(ii)** multiplicativity of x ↦ x^Q / rational points are Frobenius-fixed, **(iii)** Riemann–Roch / dimension count, **(iv)** extensions F_{q^n} (tensor-power trick), **(v)** functional equation.

| step | statement | ingredient | check (all 6 curves + 6 twists) |
|---|---|---|---|
| S1 | L(mO) = span{xⁱyʲ : 2i + j(2g+1) ≤ m}, dim = m + 1 − g for m ≥ 2g − 1; gaps {1} (g=1), {1,3} (g=2) | (iii) | CHECKED |
| S2 | u ⊗ v ↦ u^p v^q is injective on L(l) ⊗ L(m) when l < p: pole orders p·a + q·b are all distinct | (iii) + base-p digits | CHECKED |
| S3 | δ: Σ c u^p v^q ↦ Σ c u^p v ∈ L(lp + m) has a kernel when dim L(l)·dim L(m) > dim L(lp + m) | (iii) | kernel found by Gauss elimination over F_q, δ(F) = 0 exactly |
| S4 | v^q(P) = v(P)^q = v(P) at every F_q-point, so F(P) = δ(F)(P) = 0 | **(ii)** | F vanishes at all N − 1 affine points |
| S5 | F = H^p with H = Σ c^p u v^p, since (Σ)^p = Σ(·)^p; hence ord_P F ≥ p | **(i)** | F equals H^p computed as a literal p-fold ring product, exact; monomial Frobenius formula = ring power |
| S6 | zeros = poles: p(N − 1) ≤ pole_O(F) ≤ lp + mq, so N ≤ 1 + l + mp | deg div = 0 | CHECKED |
| S7 | N(C) + N(C^d) = 2(q + 1) (count over the x-line); an upper bound on C^d is a lower bound on C | Galois cover / twist | CHECKED |
| S8 | N_n = q^n + 1 − Σα_iⁿ (rationality from RR), and S1–S7 over every F_{q^n} with the **same** constant 2g + 1. Then max\|ω\| ≤ (2g+1)^{1/2n} → 1, and the FE α ↔ q/α gives \|ω\| ≥ 1 | (iv), (iii), (v) | N₂ from the α's equals the brute count; \|N_{2n} − q^n − 1\| ≤ (2g+1)q^{n/2} for n ≤ 15; \|ω\| = 1; FE pairing exact |

Numbers (`stepanov_ff.json`). (l, m) is the cheapest choice satisfying S2+S3.

| curve | q | (l, m) | dims L(l), L(m), L(lp+m) | kernel dim | pole(F) ≤ lp+mq | N | upper bound | Bombieri q+1+(2g+1)√q | Weil | twist N′ | lower bound | Bombieri lower |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E p=5 | 25 | (4, 7) | 4, 7, 27 | 1 | 190 ≤ 195 | 27 | 39 | 41 | 36 | 25 | 13 | 11 |
| E p=7 | 49 | (5, 9) | 5, 9, 44 | 1 | 469 ≤ 476 | 60 | 68 | 71 | 64 | 40 | 32 | 29 |
| E p=11 | 121 | (7, 13) | 7, 13, 90 | 1 | 1639 ≤ 1650 | 140 | 150 | 155 | 144 | 104 | 94 | 89 |
| E p=13 | 169 | (8, 15) | 8, 15, 119 | 1 | 2626 ≤ 2639 | 195 | 203 | 209 | 196 | 145 | 137 | 131 |
| G2 p=11 | 121 | (10, 15) | 9, 14, 124 | 3 | 1892 ≤ 1925 | 130 | 173 | 177 | 166 | 114 | 71 | 67 |
| G2 p=13 | 169 | (11, 17) | 10, 16, 159 | 1 | 2990 ≤ 3016 | 207 | 231 | 235 | 222 | 133 | 109 | 105 |

Every bound holds and every intermediate identity is exact. Over the extensions, max over n ≤ 15 of |N_{2n} − q^n − 1|/q^{n/2} is 2.0–2.2 for g = 1 and 3.3–3.6 for g = 2. That is ≤ 2g, and in particular below Bombieri's 2g + 1.

**Fakes on the F_q side.** None of these is a curve, so S1–S6 cannot even be set up. Their "N_n" break every C·q^{n/2}:
- the class-0 partial zeta Z₀ of probe B has |ω| = 1.618 / 0.618 at q = 5, and |Σωⁿ| > 3 from n = 3 on;
- the synthetic genus-2 fake (FE, one off-circle pair) passes > 5 at n = 2, 6 or 13.

**Ablation, i.e. which ingredient carries the weight (`ablation_ff.json`, CHECKED):**

| curve | full proof | A: u^p → u^{p+1} (multiplicative, not additive) | B: v^q → v^p (additive + multiplicative, but does not fix F_q-points) | C: m − 2 (RR count too small) |
|---|---|---|---|---|
| E p=7 | vanishes at 59/59, N ≤ 68 | no kernel at all | vanishes at 7/59 | no kernel |
| E p=11 | 139/139, N ≤ 150 | no kernel | 15/139 | no kernel |
| G2 p=11 | 129/129, N ≤ 173 | kernel, but multiplicity only ≥ 1: N ≤ 1779 (> 14q) | 13/129 | no kernel |

So the square-root gain is the vanishing multiplicity p = √q from S5 (additivity). The "everywhere at once" comes from S4: one Frobenius fixes all rational points. **That is Daniel's coupling made exact: one map that is additive and multiplicative at once, and whose fixed points are exactly the objects being counted.**

---

## Part 2 — Port map F_q → ℤ

| step | F_q ingredient | ℤ candidate | real analogue? |
|---|---|---|---|
| S1 RR count | dim L(mO) = m + 1 − g | (a) Paley–Wiener / prolate count: window width L, band T, dim ≈ LT/π (Landau–Pollak–Slepian); (b) **Siegel's lemma / Minkowski**, the RR of geometry of numbers; (c) Connes–Consani RR for Spec ℤ̄ (2023) and for ℤ over 𝕊 (2024); (d) CC RR on the periodic orbit C_p = ℝ₊^×/p^ℤ: Dim_ℝ H⁰(D) − Dim_ℝ H⁰(−D) = deg D (genus-1 shape, real-valued) | **yes** (several) |
| S2 injectivity | distinct pole orders p a + q b | uniqueness of base-p expansions; CC ℤ_min ⊗_𝔹 ℤ_min (Newton polygons) | partial |
| S3 kernel | linear algebra over F_q | Siegel's lemma: more unknowns than equations ⇒ a small non-zero integer solution | **yes** |
| S4 v^q(P) = v(P) at all points | Frobenius multiplicative, points = Frobenius fixed points | Fermat a^p ≡ a (mod p), **one p at a time**. Globally: primes = periodic orbits of the scaling flow (Connes 1999 trace formula = Lefschetz), but no ring in which "evaluation at all orbits at once" is a ring map | **only per prime** / char 1 |
| S5 F = H^p | x ↦ x^p **additive** in char p | mod p: (a+b)^p ≡ a^p + b^p. Borger λ-rings: Frobenius lifts ψ_p are ring maps, but the identity on ℤ. Char 1: Fr_λ(x) = x^λ on ℝ_max is additive (max) and multiplicative and scales divisors by λ (CC 2016–17) | **no in char 0**; exists in char 1 |
| S6 zeros = poles | deg div f = 0 | **product formula** \|N\|_∞ ∏_p \|N\|_p = 1 for N ∈ ℚ^× ("divisible by many p^k ⇒ archimedean size large"); also deg div = 0 on C_p | **yes** |
| S7 twist / lower bound | N(C) + N(C^d) = 2(q + 1) | **Landau's theorem + no real zero**. A one-sided bound ψ_k(x) ≤ x∫k + O(x^θ) plus ζ(σ) ≠ 0 on (θ,1) gives no zeros with Re ρ > θ (Landau; Ingham 1932, ch. V). ζ(σ) < 0 on (0,1) is PROVED (η-series), so the lower bound is free over ℤ. Characters (Σ_χ orthogonality) are the other twist analogue. The sieve parity barrier blocks the sieve-style lower bound | **yes** (Landau route); measured in §3.4 |
| S8 extensions + tensor power | C/F_{q^n}: same genus, constant 2g+1 uniform in n | formal: τ_{ns} = τ_sⁿ exactly (x ↦ xⁿ), and the Laplace / Landau amplification works (Round 15 §3.2). But there is no constant-field extension of ℚ. The cyclotomic tower ℚ(μ_n) ("F_{1ⁿ}") has genus ∝ log disc ≈ φ(n) log n → ∞. Scaling L → nL is not arithmetic: it sends log p to log p / n | **no** (formal only) |
| S9 FE | α ↔ q/α | ζ(s) = ζ(1 − s) completed, x → 1/x (Daniel). Lets a half-plane bound give the line | **yes** |
| — | sieve / dimension upper bounds | Beurling–Selberg / Vaaler extremal majorants (one-sided, band-limited "functions of bounded degree") + large sieve (Montgomery–Vaughan: N − 1 + Q²): Brun–Titchmarsh π(x+y) − π(x) ≤ 2y/log y. This is S1+S3+S6 **without S4/S5**, and it loses a factor 2 (parity) instead of gaining √ | yes, but too weak |

**Read-off.** Everything except S5 (and the uniform S8) has a genuine ℤ shadow. S5 is the one place where "add and multiply at the same time" does real work, and it is the one place ℤ has nothing in characteristic 0. In characteristic 1 Fr_λ exists, and C_p even has a genus-1-shaped RR. What is missing there is S4 for *all primes at once*: each C_p is its own circle of length log p, with incommensurate clocks (Round 15 §4(ii)). HEURISTIC: no published Stepanov-type argument on the arithmetic or scaling site was found.

---

## Part 3 — Trying it on ℤ (`zport.py` → `zport.json`)

### 3.1 Normalisations and the best provable θ

For g in the l = 0.6 window, A = autocorrelation of g, x = e^s, and k(y) = A(log y) y^{−1/2}:
b_s(g) = W(g, τ_s g) = pole_s − prime_s + arch_s, with prime_s = Σ Λ(n) n^{−1/2} A(log n − s) = x^{−1/2} ψ_k(x).

So |b_s − pole_s| ≲ C e^{θ_W s} ⇔ ψ_k(x) − x∫k = O(x^θ) with θ = ½ + θ_W. RH ⇔ θ = ½ (Round 15 §3.2). Trivial: θ = 1.

| method (zero-free, Stepanov-type) | θ proved | label |
|---|---|---|
| Chebyshev binomial auxiliary integer C(2x, x): ψ(2x) − ψ(x) ≤ 2x log 2. In b_s terms: \|prime_s\| ≤ 2 log 2·e^{3l/2}·e^{s/2}‖g‖² | **θ = 1** | PROVED (here, elementary) |
| Diamond–Erdős (1980) factorial-ratio auxiliary integers: constant 1 + ε, never o(x) | θ = 1 | literature (recalled) |
| integer Chebyshev / Gelfond–Nair lower bounds ∫P^{2k}: capped at ψ(x) ≥ 0.8643x, since t_ℤ[0,1] > 0.4213 (Montgomery, *Ten lectures*, ch. 10). The method cannot reach PNT | θ = 1 (lower) | literature |
| large sieve / Selberg (dimension count + Beurling–Selberg majorants): window upper bound with constant 2 | θ = 1 | literature |
| zero-free region (not Stepanov-type): Vinogradov–Korobov saving exp(−c s^{3/5}(log s)^{−1/5}) on e^{s/2} | θ = 1, subexponential saving | literature |
| any θ < 1 | would be a zero-free strip | **OPEN** |

**Best provable θ = 1** (θ_W = ½), i.e. exactly the trivial square-root *loss*. In ζ's measured data (l = 0.6, sine1) the proven Chebyshev bound exceeds the true |b_s| by a factor of 600 (s = 1) to 11,000 (s = 7) relative to d. The true |b_s|/d stays ≤ 0.89. The tensor-power trick amplifies **constants**, not exponents: e^{θ_W s} for all s with θ_W > 0 yields only the strip Re ρ ≤ ½ + θ_W. So the whole gain e^{−s/2} must come from S5-type high-order vanishing. Over F_q that is exactly the multiplicity √q.

### 3.2 Measured exponents (Round-15 blocks, ρ(s) = sup over g of \|b_s\|/d)

| model | ρ-envelope slope θ_W (s ∈ [2, 8]) | max ρ | sine1 θ_W | what it lacks |
|---|---|---|---|---|
| ZETA | **0.000** | 0.9996 | 0.009 | — |
| FTAU | 0.002 | 0.996 | 0.004 | (unitary control, no off-line zeros) |
| EPS | 0.363 | 4.21 | −0.08 | Euler product |
| DH | 0.123 | 1.50 | 0.04 | Euler product |
| BEUR2 | (jump at the splice, slope −0.17) | 9.66 | −0.86 | Poisson summation on ℤ (Beurling integers ≠ ℤ) |
| FTH | 0.157 (true δ = 0.2) | 2.54 | 0.149 | unitarity |

NUMERICAL. ζ is flat: θ_W = 0 empirically, which is the RH value, but measured, not proved. F_θ's slope recovers its built-in δ = θ = 0.2 to within the fit.

### 3.3 Do the Stepanov-type ingredients fail on the fakes?

| model | Chebyshev identity Σ c(n)⌊x/n⌋ − log⌊x⌋!, /x at x = 2700 | binomial identity residual /x | (ψ(2x) − ψ(x))/(2x log 2) ≤ 1? | MV max ratio y = √x / x^{3/4} / Hasse width | Landau: real zero on (½, 1)? |
|---|---|---|---|---|---|
| ZETA | **0 (exact)** | **0** | 0.714 ✓ | 0.46 / 0.46 / 0.48 ✓ | none (ζ(σ) ∈ [−199, −1.48]) |
| FTAU | −7.89 | −0.068 | 1.149 ✗ | 0.75 / 0.82 / 0.72 ✓ | none (\|ζ(σ+3i)\|² > 0) |
| EPS | 0.224 | −0.087 | 0.524 ✓ | **1.82 / 1.11** / 0.61 ✗ | none (Z_Q(σ) ∈ [−280, −2.3]) |
| DH | −7.09 | −1.37 | 0.014 ✓ | 0.32 / 0.19 / 0.08 ✓ | none (∈ [0.83, 0.92]) |
| BEUR2 | −0.016 | −0.003 | 0.719 ✓ | 0.40 / 0.43 / 0.48 ✓ | N.A. |
| FTH | 13.4 | 5.62 | **3.88 ✗** | **1.99 / 1.93 / 2.30 ✗** | real **pole** at 1 − θ = 0.8 (main term), no zero |

EPS check: the formula Z_{x²+5y²} = ζ·L(χ₋₂₀) + L(χ₋₄)L(χ₅) matches the lattice sum at s = 2 (2.50242 vs 2.50216 truncated).

**Readout (NUMERICAL, the identities CHECKED).**
- The **exact add+mult identity** Λ * 1 = log (unique factorisation in ℤ: an additive count, log⌊x⌋! by Stirling, equals a multiplicative factorisation) holds **only for ζ**. It fails for every fake, *and* for the unitary control FTAU. It is the ℤ-shadow of S4+S5: one identity valid at all primes at once.
- Its **inequality** consequences do not separate:
  - Chebyshev's θ = 1 bound holds for EPS, DH and BEUR″;
  - MV / large sieve fails only for EPS (weight 2 log p on split primes) and FTH;
  - DH and BEUR″ pass both.
  
  So every *proven* Stepanov-type step is blind to DH and BEUR″. A genuine port must fail on BEUR″ through Poisson on ℤ (BEUR″ has an Euler product but its integers are not ℤ) and on DH/EPS through the Euler product. The θ = 1 arguments use neither in a way that survives to the √ scale.
- The **Landau step** (the S7 analogue) is free for ζ but does not separate either: EPS and DH have no real zeros on (½, 1) either.

### 3.4 Tensor-power analogue: does a bound at s imply one at ns?

The sequence is ρ(n s₀) along progressions s₀ ∈ {0.6, 0.8, 1.0, 1.2, 1.4}; ρ ≤ 1 is the Hasse bound.

| model | holds at n = 1 (all s₀) | first n with ρ(n s₀) > 1 | Hankel numerical rank of n ↦ b_{0.6 + 0.02n} (K = 184) |
|---|---|---|---|
| ZETA | yes | never (n s₀ ≤ 8) | 177 |
| FTAU | yes | never | 165 |
| EPS | yes | 8, 6, 5, 4, 4 | 149 |
| DH | yes | 11, —, 7, —, 5 | 141 |
| BEUR2 | yes | 6, 5, 4, 3, 3 | 182 |
| FTH | yes | 3, 5, 2, 3, 3 | 180 |
| F_q curve, genus g | — | never (Weil) | **2g** exactly (power sum of 2g terms) |

NUMERICAL.
- Every fake satisfies the bound at s₀ and violates it later along the same progression. So "bound at s ⇒ bound at ns" is **false** without extra structure.
- Over F_q the structure is finite genus: a_n is a power sum of 2g terms, and the bound for every n is re-proved over F_{q^n} with the same constant.
- Over ℤ, n ↦ b_{ns}(g) is also formally a power sum, b_{ns} = Σ_ρ G(γ_ρ)G(−γ_ρ) e^{−iγ_ρ ns}, so the amplification half of the trick survives. But it has infinitely many terms (Hankel rank essentially full), and there is no "ℤ over a bigger constant field" on which to re-prove the input with a uniform constant.
- A matrix-pencil growth estimate was also run. It is ill-conditioned (FTAU shows a spurious δ), so it is recorded in the JSON and **not** used as evidence.

---

## Part 4 — Verdict

### 4.1 Where the port holds and where it breaks

- **Holds (real ℤ analogues):**
  - S1/S3 dimension count: Siegel's lemma; prolate/Paley–Wiener; CC RR for Spec ℤ̄ and on C_p;
  - S6 zeros = poles: the product formula;
  - S7 lower bound: Landau + ζ(σ) < 0 on (0,1), PROVED and easier than over F_q;
  - S9 the FE (x → 1/x);
  - and the amplification half of S8 (Round 15 §3.2; Landau).
- **Breaks:**
  - **S5**: no characteristic-0 map that is additive and multiplicative and makes an auxiliary object vanish to high order at all primes simultaneously. Mod p it exists one prime at a time, and in characteristic 1 (Fr_λ) it exists only per orbit C_p.
  - **S8-input**: no constant-field extension, so no uniform constant.
- **Consequence:** every proven Stepanov-type bound over ℤ stops at θ = 1 (the e^{s/2} loss). The missing e^{−s/2} is exactly the multiplicity √q that S5 supplies over F_q.

### 4.2 The GAP, as one lemma

**Lemma S_ℤ (integer Stepanov lemma; GAP, RH-strength).** For every ε > 0 there is C_ε such that for every x ≥ 2 there are an integer N_x ≠ 0 and an integer μ_x ≥ 1 with
1. (vanishing) p^{μ_x} | N_x for every prime p ∈ (x e^{−l}, x e^{l}] (the S4+S5 analogue: high-order vanishing at all window primes at once);
2. (size) log|N_x| ≤ μ_x·( x(e^{l} − e^{−l}) + C_ε x^{1/2+ε} ) (the RR/pole budget, with square-root-sharp slack);

constructed without any information on zeros of ζ.

**Why it is the right lemma.**
- (1) + (2) + the product formula (S6) give Σ_{window} log p ≤ log|N_x|/μ_x = x(e^{l} − e^{−l}) + O(x^{1/2+ε}), one-sided. That is Bombieri's N ≤ q + O(√q) for ℤ.
- Landau plus ζ(σ) < 0 on (0,1) (S7) upgrades the one-sided bound to "no zeros with Re ρ > ½ + ε", i.e. RH (PROVED sketch). In the Hasse-block language this is b_s ≥ −C_ε e^{εs}.
- Status of the known constructions:
  - Chebyshev's N_x = C(2x, x) satisfies the lemma with μ = 1, but its slack is ≈ 0.39x (from 2 log 2 against 1) instead of x^{1/2+ε};
  - Diamond–Erdős factorial ratios bring the slack to εx;
  - nothing known reaches o(x) by this route.
- Over F_q the lemma is true with μ = √q, from S5. Over ℤ it would need an "additive Frobenius" that raises multiplicity at all window primes at once while paying archimedean size only μ·x. **OPEN.**

### 4.3 One best next step

**Round 17: the characteristic-1 Bombieri on the scaling site, one prime at a time, then glue.**
1. Run S1–S6 literally on Connes–Consani's periodic orbit C_p, with real-valued Dim_ℝ RR (genus-1 shape, like probe-B elliptic curves) and Fr_λ as the additive+multiplicative Frobenius. Build the tropical auxiliary function Σ c ⊗ u^{λ} v^{λ'} and see what it bounds. This is a finite, checkable piecewise-linear computation.
2. Then test whether two orbits C_p, C_{p′} can share one auxiliary function. That is the incommensurate-clock obstruction, and exactly where S4 for "all primes at once" must appear.

If even the single-orbit version gives nothing, then S5 has no shadow in characteristic 1 either, and the Stepanov route over ℤ is closed.

## Files (this folder)

| file | role |
|---|---|
| `ffq.py` | F_{p²} and hyperelliptic function-field arithmetic (exact) |
| `stepanov_ff.py` → `stepanov_ff.json`, `.log` | Part 1, S1–S8 on 6 curves + 6 twists, F_q fakes |
| `ablation_ff.py` → `ablation_ff.json`, `.log` | Part 1 ablations A/B/C |
| `zport.py` → `zport.json`, `.log` | Part 3 (a)–(f): envelopes, Chebyshev/binomial identities, MV, Landau real axis, progressions/Hankel, proven bound |

Reproduce: `python3 stepanov_ff.py; python3 ablation_ff.py; python3 zport.py` (zport reads `../round15_hasse/blocks_*_l0.6.json`).

## References

- E. Bombieri, *Counting points on curves over finite fields (d'après S. A. Stepanov)*, Sém. Bourbaki 430 (1973).
- S. A. Stepanov, *On the number of points of a hyperelliptic curve over a finite prime field*, Izv. (1969); W. M. Schmidt, *Equations over finite fields* (1976).
- H. Iwaniec and E. Kowalski, *Analytic Number Theory*, ch. 11.
- T. Tao, *The Bombieri–Stepanov proof of the Hasse–Weil bound* (blog, 2014); K. Kedlaya, *Weil cohomology* notes, ch. 5.
- P. L. Chebyshev (1852); H. G. Diamond and P. Erdős, *On sharp elementary prime number estimates*, Enseign. Math. (1980).
- H. L. Montgomery, *Ten lectures on the interface between analytic number theory and harmonic analysis* (1994), ch. 10 (integer Chebyshev).
- H. L. Montgomery and R. C. Vaughan, *The large sieve*, Mathematika 20 (1973).
- J. D. Vaaler, *Some extremal functions in Fourier analysis*, BAMS (1985).
- E. Landau (1905), A. E. Ingham, *The distribution of prime numbers* (1932), ch. V: one-sided bounds and real singularities.
- A. Connes and C. Consani:
  - *The scaling site*, CRAS 354 (2016), with RR on C_p;
  - *Geometry of the arithmetic site* (2016);
  - *Riemann–Roch for Spec ℤ̄* (2023), and *for the ring ℤ* (2024).
- J. Borger, *Λ-rings and the field with one element* (2009).

RH / ARP / H_ψ: **OPEN**.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
