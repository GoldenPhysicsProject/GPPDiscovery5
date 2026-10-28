# PROLATE_BOUNDS — archimedean prolate/Sonin spectral bounds for Λ ≥ 1 (inputs for Queen's M1)

Seat: Spectrum & Measure. **Everything here is archimedean, local and RH-silent** (no prime, no ℚ*-sum, no E). **RH is OPEN.**
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Date: 2026-10-28.
Scope: eigenvalue asymptotics, plunge counts, explicit nonasymptotic tail/plunge/count bounds with constants, and explicit HS/trace
norms of the local (1 − P_Λ)-pieces, in Queen's normalization (team/queen/ARP_queen_Q8.md: B^P_{Λ,0} = span of even prolates
ψ_n ⊗ 1_Ẑ, n ≤ 4Λ², Connes (31)–(33)). **Queen does the step from these bounds to the leak (1 − S_Λ)Q′_Λ after E; not done here.**
Labels: PROVED (verified citation, or elementary and proved here) / via restatement / CONJECTURE / GAP. [Λ-track] throughout.
R-line (whole file): R1 every Λ ≥ 1 | R2 not engaged (archimedean only) | R3 — | R4 — | R5 fake-blind: every object with ζ's archimedean place has the same numbers.

## 1. Normalization. PROVED (definitions; Connes 1999 verified)
- Connes (31): "P_Λ = {ξ ∈ L²(ℝ); ξ(x) = 0 ∀|x| > Λ}", P̂_Λ = F P_Λ F⁻¹ with α(x) = e^{−2πix}. A_Λ := P_Λ P̂_Λ P_Λ has kernel
  K_Λ(x, y) = sin(2πΛ(x − y))/(π(x − y)) on [−Λ, Λ]².
- Karnik–Romberg–Davenport (KRD, arXiv 2006.00427, ACHA 2021; verified from the arXiv source): "(T_T B_Ω T_T y)(t) = ∫_{−T/2}^{T/2}
  sin[Ω(t−t′)]/(π(t−t′)) y(t′) dt′ … the eigenvalues … only depend on the product ΩT … we use the notation λ̃_k(c) … ΩT/2 = c."
  Matching: T = 2Λ, Ω = 2πΛ, so **c = 2πΛ², 2c/π = 4Λ², 100c/π = 200Λ²**. Write λ_n = λ_n(Λ) = λ̃_n(2πΛ²), decreasing, 0 < λ_n < 1.
- Rescaling (LAMBDA_L_dictionary (a2); Principal series P.9.2(a)): only c = Λ² · 2π matters; (P_Λ, P̂_Λ) ≅ (P₁, P̂_{Λ²}).
- Connes' transition-band statement (verified, math/9811068 after (32)): "eigenvalues λn of the operator PΛ P̂Λ PΛ are decreasing very
  slowly from λ0 ≃ 1 until the value n ≃ 4Λ² of the index n, they then decrease from ≃ 1 to ≃ 0 in an interval of length ≃ log(Λ) and
  then stay close to 0." §§2–4 make this explicit.
- Even sector (Queen uses even prolates): eigenvalues of A_Λ on L²(ℝ)ev are the λ_n with n even. Every even-sector sum below is
  bounded by the full sum (subset of nonnegative terms).

## 2. Exact traces. PROVED (elementary)
- Tr A_Λ = ∫_{−Λ}^{Λ} K_Λ(x, x) dx = 2Λ · 2Λ = **4Λ²** (= 2c/π; Connes' 4Λ²).
- Even sector: Tr(A_Λ|_ev) = ∫_{−Λ}^{Λ} ½[K_Λ(x,x) + K_Λ(x,−x)] dx = **2Λ² + Si(4πΛ²)/(2π)**. At Λ = 1: 2.23748…, matching CC 2021
  footnote 7 ("the even λ(n)² sum to 2.237…", as quoted in Principal series P.9.1 (D4); CC's λ(n)² are our even λ_n at Λ = 1).

## 3. Exact plunge mass with explicit constant (Queen's ask). PROVED
**Identity.** Σ_n λ_n(1 − λ_n) = Tr(A − A²) = ‖P_Λ P̂_Λ (1 − P_Λ)‖²_HS, because A − A² = P P̂ (1 − P) P̂ P = X X* with X = P P̂ (1 − P).
**Closed form.** With x := 8πΛ² (= 4c) and si(x) := ∫_x^∞ (sin v / v) dv,
  D(Λ) := Σ λ_n(1 − λ_n) = (1/π²) [ log x + 1 + γ − Ci(x) − cos x + x·si(x) ].
Derivation: D = ∫_{|x|≤Λ}∫_{|y|>Λ} K_Λ² = (2/π²) ∫_0^∞ min(v, 2Λ) sin²(2πΛv)/v² dv = (1/π²)[Cin(x) + x∫_x^∞ (1 − cos w)/w² dw],
Cin(x) = γ + log x − Ci(x), and x∫_x^∞(1 − cos w)/w² dw = 1 − cos x + x·si(x) (integration by parts).
(Identity check only, not exploration: the closed form matches Σλ(1−λ) from a 900-node Gauss–Legendre discretization to 10⁻¹² at
Λ = 0.7, 1, 1.5.)
**Explicit remainder.** R(x) := −Ci(x) − cos x + x·si(x) satisfies **|R(x)| ≤ 2/x² + 4/x³** for x > 0.
Proof: with the auxiliary functions f(x) = ∫_0^∞ e^{−xt}/(1+t²)dt, g(x) = ∫_0^∞ t e^{−xt}/(1+t²)dt one has Ci = f sin − g cos and
si = f cos + g sin (DLMF 6.2.17–6.2.18 / 6.7.13–6.7.14 type representations; standard), so R = cos x·(xf − 1 + g) + sin x·(xg − f).
Here xf − 1 = −∫ x e^{−xt} t²/(1+t²) dt ∈ [−2/x², 0], g ∈ [0, 1/x²], and (integrating xg by parts) xg − f = −∫ e^{−xt} 2t²/(1+t²)² dt
∈ [−4/x³, 0]. ∎ (DLMF numbers given as the standard source of the f, g representation; the bounds are proved here.)
So, for all Λ > 0:
  **D(Λ) = (1/π²)(log(8πΛ²) + 1 + γ) + θ·(2/x² + 4/x³)/π², |θ| ≤ 1, x = 8πΛ².**
At Λ = 1: D = 0.48633 (x = 25.13, remainder ≤ 3.5·10⁻⁴). Leading term (1/π²) log c agrees with Landau–Widom (§5): Landau's log c growth
with all constants explicit.
**HS form for Queen.** ‖P_Λ P̂_Λ (1 − P_Λ)‖²_HS = ‖(1 − P_Λ) P̂_Λ P_Λ‖²_HS = D(Λ); by F-conjugation also ‖P̂_Λ P_Λ(1 − P̂_Λ)‖²_HS = D(Λ).

## 4. Nonasymptotic consequences of §3. PROVED (elementary)
Let N_½ := #{n : λ_n > ½}, M := 4Λ².
- Tail trace: Σ_{n ≥ N_½} λ_n ≤ 2D (λ ≤ 2λ(1−λ) for λ ≤ ½); head defect Σ_{n < N_½}(1 − λ_n) ≤ 2D; jointly the two sums are ≤ 2D.
- Index location: |N_½ − M| ≤ 2D (N_½ − M = Σ_{n<N_½}(1−λ_n) − Σ_{n≥N_½}λ_n).
- Count (Chebyshev): #{n : ε ≤ λ_n ≤ 1 − ε} ≤ D/(ε(1 − ε)). (Weak in ε; §5 has the log(1/ε) form.)
- With §3: all three are ≤ (2/π²)(log(8πΛ²) + 1 + γ) + 2(2/x² + 4/x³)/π², e.g. ≤ 0.973 at Λ = 1.

## 5. Literature bounds (verbatim)
**KRD (arXiv 2006.00427), verified from the LaTeX source** (their numbering in the arXiv version):
- Theorem (PSWF eigenvalue gap, `thm:PSWFEigGap`, Theorem 3 in the paper): "For any c > 0 and ε ∈ (0, ½),
  #{k : ε < λ̃_k < 1 − ε} ≤ (2/π²) log(100c/π + 25) log(5/(ε(1 − ε))) + 7."
- Corollary (`cor:PSWFEigBounds`, Cor 3): "For any c > 0 and ε ∈ (0, ½), λ̃_k ≥ 1 − 10 exp[ −(⌊2c/π⌋ − k − 7) / ((2/π²) log(100c/π + 25)) ]
  for 0 ≤ k ≤ ⌊2c/π⌋ − 1, and λ̃_k ≤ 10 exp[ −(k − ⌈2c/π⌉ − 6) / ((2/π²) log(100c/π + 25)) ] for k ≥ ⌈2c/π⌉."
- Corollary (`cor:PSWFEigSumBounds`, Cor 4): "For any c > 0, Σ_{k=0}^{K−1} (1 − λ̃_k) ≤ (20/π²) log(100c/π + 25) exp[ −(⌊2c/π⌋ − K − 7) /
  ((2/π²) log(100c/π + 25)) ] for 1 ≤ K ≤ ⌊2c/π⌋, and Σ_{k=K}^∞ λ̃_k ≤ (20/π²) log(100c/π + 25) exp[ −(K − ⌈2c/π⌉ − 7) /
  ((2/π²) log(100c/π + 25)) ] for K ≥ ⌈2c/π⌉."
- KRD remark (verbatim): the bound "scales like O(log(c) log(1/ε)) and has the correct leading constant 2/π²".
**In Queen's normalization** (c = 2πΛ², M = 4Λ², β_Λ := (2/π²) log(200Λ² + 25)):
- count: #{n : ε < λ_n < 1 − ε} ≤ β_Λ log(5/(ε(1−ε))) + 7;
- tail trace: Σ_{n ≥ K} λ_n ≤ 10 β_Λ exp(−(K − ⌈M⌉ − 7)/β_Λ) for K ≥ ⌈M⌉;
- head defect: Σ_{n < K}(1 − λ_n) ≤ 10 β_Λ exp(−(⌊M⌋ − K − 7)/β_Λ) for 1 ≤ K ≤ ⌊M⌋;
- termwise: 1 − λ_k ≤ 10 exp(−(⌊M⌋ − k − 7)/β_Λ) (k < ⌊M⌋), λ_k ≤ 10 exp(−(k − ⌈M⌉ − 6)/β_Λ) (k ≥ ⌈M⌉).
- Example Λ = 1: β_1 = (2/π²) log 225 = 1.0975; Σ_{n≥15} λ_n ≤ 10.975·e^{−4/1.0975} ≈ 0.29; Σ_{n≥20} λ_n ≤ 3.0·10⁻³.
These are the "explicit nonasymptotic bounds with constants" Queen asked for. PROVED (KRD, verified). Even-sector sums ≤ these.
**Bonami–Jaming–Karoui (arXiv 1804.01257, J. Math. Phys. 62 (2021) 033511), verified**, for Q_c f(x) = ∫_{−1}^{1} sin c(y−x)/(π(y−x)) f(y) dy
(same λ_n(c)):
- (1.2) "For any c > 0 and any 0 ≤ n < 2c/2.7, 1 − (7/√c)((2c)^n/n!) e^{−c} ≤ λn(c) < 1."
- (1.3) "For any c ≥ 22 and η = 0.069, we have, for 2c/π + log c + 6 ≤ n ≤ c, λn(c) ≤ exp(−(η/2)·(n − 2c/π)/(log c + 5))."
- Restated there (so **via restatement**): Widom (2.7) "λn(c) ∼ (ec/(4(n + ½)))^{2n+1}"; Landau–Widom (2.8), Λ_ε = {k : λk(c) ≥ ε},
  "|Λε| = 2c/π + (1/π²) log((1−ε)/ε) log(c) + o(log c)"; Osipov "|Λε| ≤ 2c/π + K(log c)²" (K unspecified); Israel's bound with K_η.
- Slepian / Landau–Pollak (1960–78): via restatement (KRD introduction); not quoted further.

## 5′. Local (1 − P_Λ)-pieces on the prolate span (inputs for Queen's M1). PROVED (elementary + KRD)
Let ψ_n be the L²(ℝ)-normalized prolates (bandlimited to |ξ| ≤ Λ, eigenfunctions of A_Λ after restriction), Π_N the projection onto
span{ψ_n : n < N} (Queen's B^P with N ≈ 4Λ², even n). Double orthogonality ⟨P_Λψ_n, P_Λψ_m⟩ = λ_n δ_nm gives
⟨(1−P_Λ)ψ_n, (1−P_Λ)ψ_m⟩ = (1 − λ_n) δ_nm, so (1 − P_Λ)Π_N has singular values √(1 − λ_n), n < N:
- ‖(1 − P_Λ)Π_N‖²_HS = Σ_{n<N} (1 − λ_n) ≤ 10 β_Λ exp(−(⌊M⌋ − N − 7)/β_Λ) (KRD Cor 4), and ≤ 2D + #{n < N : λ_n ≤ ½} in general.
- ‖(1 − P_Λ)Π_N‖₁ = Σ_{n<N} √(1 − λ_n) ≤ √10 · exp(−(⌊M⌋ − N − 7)/(2β_Λ)) / (e^{1/(2β_Λ)} − 1) (KRD Cor 3 termwise, geometric sum),
  for N ≤ ⌊M⌋; and by Cauchy–Schwarz ≤ √(N · Σ_{n<N}(1 − λ_n)).
- ‖(1 − P_Λ)Π_N‖ = √(1 − λ_{N−1}).
- If N is pushed to the plunge (N ≈ M), these are O(log Λ)-size: Σ_{n<N}(1−λ_n) ≤ 2D + (number of n < N with λ_n ≤ ½) — the trace-norm leak
  of the full prolate span is **not** small; it is ~ (2/π²)log Λ in HS² and needs the choice N ≤ M − C β_Λ log(1/δ) to be ≤ δ.
  That trade-off is the quantitative content Queen needs; the passage to (1 − S_Λ)Q′_Λ after E is hers.
- The P̂-placement (Queen's (P̂)) has the same numbers by F-conjugation.

## 6. Hooks
- **Sonin side (Λ = 1, CC 2021):** S(1,1) is the eigenvalue-1 eigenspace of P P̂ P (CC (81), cited via Principal series P.9.1 (D3)).
  On the generic part the angle operator is shared by (P₁, P̂₁) and (P, P̂) (P.9.1 (D1)), so the numbers above (at c = 2π) are the
  Sonin-side spectral data too.
- **Connes–Moscovici, arXiv 2112.05500 ("Prolate spheroidal operator and zeta"), verified:** Thm 2.6: "(ii) Wsa commutes with the
  projections Pλ and P̂λ … (iv) The spectrum of Wsa is discrete and unbounded on both sides." Abstract: the restriction of W_sa to the
  complement of J has "negative eigenvalues whose ultraviolet behavior reproduce that of the squares of zeros of the Riemann zeta
  function. Furthermore, their corresponding eigenfunctions belong to the Sonin space." Hook only: a UV counting statement, RH-silent.
- **Principal series (P9.4) at Λ:** the prolate data at bandwidth 2πΛ² (ε^{(Λ)}) use exactly these λ_n(Λ); the K_I^{(Λ)} GAP
  (LAMBDA_L_dictionary (e5)) would use §§3–5′.
- **Celestial §C.13:** CELESTIAL_INPUTS_C13.md uses §5 at c = 2π.

## Status
| item | label |
|---|---|
| c = 2πΛ², Tr = 4Λ², Tr_ev = 2Λ² + Si(4πΛ²)/(2π) | PROVED |
| Σλ(1−λ) = Tr(A − A²) = ‖PP̂(1−P)‖²_HS, closed form, |R| ≤ 2/x² + 4/x³ | PROVED |
| tail/head ≤ 2D, |N_½ − 4Λ²| ≤ 2D, count ≤ D/(ε(1−ε)) | PROVED |
| KRD Thm 3, Cor 3, Cor 4 (explicit constants) | PROVED (verified) |
| BJK (1.2), (1.3) | PROVED (verified) |
| Widom, Landau–Widom, Osipov, Slepian | via restatement |
| (1 − P_Λ)Π_N HS / trace / operator norms | PROVED |
| passage to Queen's leak (1 − S_Λ)Q′_Λ after E | not done (Queen) |
| RH | OPEN |
