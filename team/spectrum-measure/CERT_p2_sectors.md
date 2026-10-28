# CERT_p2_sectors: certified sector lower bounds for Q_L at the p = 2 entry. Odd past ℓ = 0.028 (to ℓ = 0.21); even by-products

Spectrum & Measure seat, team/spectrum-measure. Date: 2026-10-28. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
**RH is OPEN.** Numerics are evidence, not proof. This file is archimedean plus the single prime 2 (2L ≤ log 3). It is RH-silent.
Script: cert/cert_galerkin_Q.py (python-flint 0.9.0, Arb balls). Log: cert/cert_galerkin_Q.log (prec 100). A 128-bit replay of the ℓ = 0.2 job follows in the next commit.

Notation as in CERT_p2_constants.md (7906287): a = log 2, c₂ = a/√2, 2L = a + ℓ, ‖ξ‖²_H = ∫(ξ′² + ¼ξ²) = (1/2π)∫(¼ + t²)|ξ̂|²,
Ω(t) = Re ψ(¼ + it/2) − log π, R(ξ) = (1/2π)∫(¼ + t²)Ω|ξ̂|², Q_L = R − 2c₂B, B(ξ) = ∫[ξ′(x)ξ′(x+a) + ¼ξ(x)ξ(x+a)]dx.

## 0. Results (R-lines)

| # | Statement | Label |
|---|---|---|
| R1 | **Odd (K2), direct:** Q_L(ξ) ≥ 0.008‖ξ‖²_H for every real odd ξ ∈ H¹₀(−L, L), all 2L ≤ log 2 + 0.21 | **CERTIFIED** |
| R2 | Odd, larger margins on smaller ranges: ≥ 0.012 for ℓ ≤ 0.2; ≥ 0.06 for ℓ ≤ 0.15 | **CERTIFIED** |
| R3 | By-product, even: m₀^e(L₀) := inf_{even ξ ∈ H¹₀(−L₀,L₀)} R(ξ)/‖ξ‖²_H ≥ 0.05 | **CERTIFIED** |
| R4 | By-product, even (K3): Q_L(ξ) ≥ 0.009‖ξ‖²_H for real even ξ, all 2L ≤ log 2 + 0.06 (superseded by the C.15 certificate in the next commit) | **CERTIFIED** |
| R5 | Monotonicity in ℓ (Lemma 1) | PROVED |
| R6 | Method (Thm S): capped-symbol Galerkin with Schur complement and the Lemma A tail | PROVED |

R1 supersedes CERT_p2_constants R5 (ℓ ≤ 0.028, margin 0.0126). The certificates use only the classical form Q_L = R − 2c₂B (Weil explicit formula,
as recorded in CERT_p2_constants §"Source of Q_L"), Prime's (Ω2) (Ω increasing on [0, ∞)), and our Lemma A (7906287). No CC numerics are used, and no Mori.

Cross-checks (NUMERICAL, not used): Prime cf4ee37 (sN/n14_KL_odd_gap) Ritz upper bounds m₀^odd ≈ 0.548, m₀^even ≈ 0.072. Our capped Galerkin
gives 0.545–0.547 (odd) and 0.0696–0.0712 (even) at L₀. Our float capped-Galerkin λ_min for odd Q_L: 0.137 (ℓ = 0.15), 0.072 (0.2), 0.034 (0.25), 0.0137 (0.3).
This is consistent with Golden Eye's odd margin at ℓ = 0.2–0.3.

## 1. Lemma 1 (monotonicity in ℓ; PROVED)

For ξ ∈ H¹₀(−L, L) with 2L ≤ log 3, ξ(x)ξ(x + a) vanishes off [−L, L − a] = S₋. So B(ξ) = ∫_ℝ[…] = (1/2π)∫(¼ + t²)cos(at)|ξ̂|² (Plancherel,
real ξ), and Q_L(ξ) = (1/2π)∫(¼ + t²)[Ω(t) − 2c₂cos(at)]|ξ̂|² ∈ (−∞, +∞] has an L-independent formula. Extension by zero embeds
H¹₀,par(−L, L) into H¹₀,par(−L′, L′) for L < L′, preserving ξ̂, Q and ‖·‖_H. So λ_par(L) := inf Q_L/‖·‖²_H is nonincreasing in L.
**One certificate at ℓ_max covers all 0 ≤ ℓ ≤ ℓ_max.** (Celestial Prop C.15.B, f13efdd, states the same; verified independently here.)
Remark: at 2L = log 3 exactly, the n = 3 Weil term sits at the boundary of supp(ξ⋆ξ*); all certificates here have 2L ≤ log 2 + 0.21 < log 3.

## 2. Theorem S (capped-symbol Galerkin/Schur certificate; PROVED)

Fix a parity, N, T₂ > 0, r > 0 and η ∈ (0, 1). Put Ω* = Ω(T₂). Let V_N be spanned by the first N H-orthonormal Dirichlet modes:
- even: f_m = cos(kx)/√(L(k² + ¼)), k = (2m+1)π/(2L);
- odd: f_m = sin(kx)/√(L(k² + ¼)), k = (m+1)π/L;
- |f̂_m(t)| = L|sinc((t−k)L) ± sinc((t+k)L)|/√(L(k² + ¼)), with + for even and − for odd.

Write ξ = u + w with u ∈ V_N and w ⊥_H V_N. Define:
- N_N = [(1/π)∫_0^{T₂}(Ω* − Ω)(¼ + t²)f̂_i f̂_j]: the capped low-band defect on V_N;
- n_⊥ = ∫_0^{T₂}(Ω* − Ω)·dens_N, with dens_N = dens_par − (1/π)(¼ + t²)Σ_{m<N}|f̂_m|² (Lemma A, projected);
- B_N = [B(f_i, f_j)] and G_E = [2∫_{S₋}(f_i′f_j′ + ¼f_i f_j)];
- β₀ = Ω* − c₂ − (1 + 1/r)n_⊥, β = (1 − η)β₀, and M = Ω*I − (1+r)N_N − 2c₂B_N − (c₂²/β)G_E.

If M ≥ λ₀I and β₀ > 0, then Q_L(ξ) ≥ min(λ₀, ηβ₀)‖ξ‖²_H on the sector.

*Proof.*
1. Ω is increasing (Ω2), so R ≥ Ω*‖ξ‖²_H − 𝒩(ξ), where 𝒩(ξ) = (1/π)∫_0^{T₂}(Ω* − Ω)(¼ + t²)|ξ̂|² ≥ 0.
2. 𝒩(ξ) ≤ (1+r)𝒩(u) + (1 + 1/r)𝒩(w).
3. Projected Lemma A: |ŵ(t)|² ≤ K_N(t)‖w‖²_H with K_N = K − Σ_{m<N}|f̂_m|², since f̂_m(t) is the H-inner product of f_m with the Green function. Hence 𝒩(w) ≤ n_⊥‖w‖².
4. S1 by Cauchy–Schwarz on the two strips S₋ and S₋ + a = −S₋ (parity): |B(u, w)| ≤ ½‖u‖_{H(E)}‖w‖_{H(E)}, with ‖u‖²_{H(E)} = uᵀG_Eu, ‖w‖_{H(E)} ≤ ‖w‖_H, and −2c₂B(w) ≥ −c₂‖w‖².
5. The cross term: 4c₂|B(u, w)| ≤ (c₂²/β)uᵀG_Eu + β‖w‖².
6. Collect, using ‖ξ‖² = ‖u‖² + ‖w‖² (H-orthogonality): Q_L(ξ) ≥ uᵀMu + ηβ₀‖w‖². ∎

This is a **direct operator lower bound**, not a Rayleigh–Ritz value. So no Temple / Lehmann–Goerisch a-priori spectral-gap input is needed. The truncation is paid for rigorously by step 3.

## 3. Arb implementation (cert/cert_galerkin_Q.py)

- Every quantity is an arb/acb ball. Ω(z) = ½[ψ(¼ + iz/2) + ψ(¼ − iz/2)] − log π is holomorphic near ℝ (poles at ±i(2n + ½)).
- N_N and n_⊥ come from acb.integral over unit-free panels (width 10 and 5) on [0, T₂], with entire sinc forms.
- B_N and G_E come from acb.integral over S₋.
- M − λ₀I is checked by ball LDLᵀ; every pivot must be certified > 0. Pivots enclose those of the true matrix, so positivity is rigorous.
- Any non-finite ball aborts the run.
- Usage: `python cert_galerkin_Q.py PREC [par ell N T2 r eta lam0]`.

## 4. Certificates (from cert/cert_galerkin_Q.log, prec 100 bits)

| sector | ℓ | N | T₂ | r | η | λ₀ | β₀ (ball) | ηβ₀ | min pivot of M − λ₀I | margin | time |
|---|---|---|---|---|---|---|---|---|---|---|---|
| odd | 0.15 | 40 | 100 | 0.003 | 0.1 | 0.06 | 2.18140… | 0.218 | > 0 | **0.06** | ~5 min |
| odd | 0.20 | 80 | 200 | 0.002 | 0.1 | 0.012 | 2.79515… | 0.280 | 0.0970 | **0.012** | 994 s |
| odd | 0.21 | 100 | 200 | 0.002 | 0.05 | 0.008 | 2.87904… | 0.144 | 0.1033 | **0.008** | 1421 s |
| even | 0 (L₀; R only, B = 0) | 40 | 100 | 0.003 | 0.5 | 0.05 | 2.22283… | 1.111 | 0.0937 | **0.05** (R3) | 199 s |
| even | 0.06 | 40 | 100 | 0.003 | 0.1 | 0.009 | 2.20701… | 0.221 | 0.163 | **0.009** (R4) | 274 s |

Ω(100) = 2.76728895…, Ω(200) = 3.46043925… (balls in the log). n_⊥ ≤ 3.5·10⁻⁴ in every run.

## 5. What limits the odd range; GAP beyond ℓ = 0.21

- The loss is the Schur cross term (c₂²/β)G_E. G_E grows with ℓ, because the H(E) mass of low modes on the strips grows.
- Float scans of this method: N = 120, T₂ = 300 gives a certificate value of +0.026 at ℓ = 0.2, and −0.011 at ℓ = 0.25.
- The sharper route is Celestial's Thm C.15.A (f13efdd). Its p = 2 term is a band-limited Fourier multiplier inside a_T, so there is no G_E cross term. Its proof is parity-free.
- An odd C.15 run at ℓ = 0.25–0.3 is the next odd step. It is paused per the lead's 2:07 order, behind the even certificate.
- **GAP:** certified odd (K2) for ℓ ∈ (0.21, log(3/2)].
- Heuristic (Golden Eye): Q_L = Σ_γ(¼ + γ²)|ξ̂(γ)|² under RH, so margins must shrink as the window grows. This is not used.

## 6. Status of the even items (pointer)

- (P) / GAP b-1 and m₀^e are off the critical path (lead, 1:41 and 2:07).
- Celestial C.14.K is **false**: NUMERICAL in C.14.4c (d447151), and Arb-refuted at (a/6, 14) in §C.15 (f13efdd, c14K_refute_arb.*). It is not used here.
- The even certificate on (0, 0.1] (target 0.2) by Thm C.15.A is the next commit.

R-line: R1 every ℓ ≤ ℓ_max via Lemma 1, constants explicit | R2 only p = 2 (through c₂cos(at)); Ω is FE/archimedean | R3 the p = 2 term is signed (odd charged) | R4 ĝ(±i/2) = 0 in Q-coordinates, so there is no pole entry | R5 fakes with ζ's W∞ and the same p = 2 term pass identically | R6 Arb balls, LDLᵀ certified.
