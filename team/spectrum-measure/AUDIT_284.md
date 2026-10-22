# AUDIT_284: independent audit and replay of the support-2.84 Mori-route certificate (L = 71/50)
Spectrum & Measure seat, 2026-10-22 (America/Toronto). RH is OPEN. This audits a finite-window statement only; numerics are evidence, not proof.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

Audited claim (S72/S74; archive/pr11_drafts/s72_comment.md and s74_comment.md; code s67/mori_cert.py, s71/s37_ldlt_exact.py, s65 Schur bound):
even parity, Ω = 1000, N = 2376 (1188 modes): λ_min(T_even) ≥ 1.0711e−81. Odd parity, Ω = 1150, N = 2732 (1366 modes): λ_min(T_odd) ≥ 2.4134e−77.
Both use s = 763/200, m = 1/2, ρ = 3. Controls: the LDLᵀ at 1.005 × the estimate fails, and odd at Ω = 1000 is ≈ −4.12e−78.

Reference: R. Mori, "Weil positivity on the window of width log 12", Zenodo 10.5281/zenodo.23265728 v2 (Prop. 3.1, Lemma 3.2, (3.3), (3.4), §5). Text in s66/prior/moriVII.txt.
Mori VIII (log 15, Zenodo 10.5281/zenodo.23281235, 2026-10-10) is the published record at Q ≥ 2^−2904‖f‖² on [−½log15, ½log15].

## 1. What is actually proved if the certificate is correct (precise statement)
- **Domain and test class.** H = L²([−L, L]; ℂ) with L = 71/50. Q is Weil's form in Mori's normalisation (1.1):
  Q(f) = 2|∫f cosh(x/2)|² − 2|∫f sinh(x/2)|² + (1/2π)∫ a(t)|Ff(t)|² dt − Σ_{n<e^{2L}} 2Λ(n)n^{−1/2} Re∫ f(x) conj f(x − log n) dx, with a(t) = Re ψ(¼ + it/2) − log π.
  Mori states the theorem for f ∈ C²(ℝ; ℂ) with supp f ⊂ [−L, L]. The inequality then extends by closure to the form domain {f ∈ H : ∫ log(2+|t|)|Ff|² < ∞}, because the right side is L²-continuous and the archimedean part is a closed nonnegative-up-to-bounded form.
- **Implication chain** (each parity separately; the reflection commutes with C_p, K, U and P):
  1. ‖C_p‖ ≤ s ⇒ mI ⪯ M = γ − C_p ⪯ (2s+m)I, with γ = s + m.
  2. Prop. 3.1 (tail): Q(f) ≥ Cδ‖f‖² + ⟨f, (M + K + U)f⟩. This needs C = log(Ω/2π) − 1/Ω − γ > 0 (this is the "log Ω > log 2π + b_L + m" condition, with b_L = s), c = ΩL > 2, 0 < δ ≤ d(c) = e^{−2(c+1)}/(4(c+1)(c+4)) (Theorem N), and u_j E_j ≤ C.
  3. Lemma 3.2 plus block elimination (3.3): G⁻¹ + J ⪰ αI + βD_M ⇒ M + V ⪰ 0 on all of H, if n < m. Here P is the projection onto polynomials of degree < N.
  4. W(X) ⪰ G = E*M⁻¹E for every X (W − G = R_X*(m⁻¹ − M⁻¹)R_X ⪰ 0), so T = W(X)⁻¹ + J − αI − β(B − A²) ≻ 0 ⇒ (3.3).
  5. Conclusion: **Q(f) ≥ C·δ‖f‖² with δ = 2^−4200, i.e. Q(f) ≥ 2^−4200.41‖f‖² ≈ 10^−1264.6‖f‖²** (my consts: log₂(Cδ) = −4200.4076) for all f in the class.
- **What "λ_min ≥ 1.0711e−81" means.** It is the margin of the finite comparison matrix T (in one parity), not a lower bound for Q. The Weil-form constant is Cδ ≈ 2^−4200. Any wording such as "Q ≥ 1.07e−81‖f‖²" would be wrong. (For S63 the orientation brief quotes Q ≥ 7.42e−78‖f‖²; that is a different method and is not audited here.)
- **Basis completeness and density.** These are *not* needed as hypotheses. The finite reduction acts on all of H through P and I − P, with explicit bounds e, n and h on the complement (Lemma 3.2). Completeness of Legendre polynomials only makes the bounds small: r_N = 2.13e−311 at N = 2376, ρ = 3.

## 2. Constants, independently recomputed (my code `replay/mori_replay.py consts`, Arb)
Even, Ω = 1000, N = 2376: γ = 4.315, b = 8.13, c = 1420, **C = 0.7538782126 > 0** (log form) and C' = a(Ω) − γ = 0.7548781709 (direct, using that a is increasing on t ≥ 0, which is elementary), −log₂ d(c) = 4123.0877 < 4200 ✓, u₀E₀ = 0.136162 ≤ C ✓, u₁E₁ = 0.134687 ≤ C ✓, r_N = 2.128e−311, K₀ = 3083.53, e = 1.15983e−307, n = 1.396e−618 < m ✓, h = 1.177e−616, α = 1.15983e−307, β = 4.6393e−307.
**All agree with the team's cert log (s71/run/cert_L7150_m12_O1000_p0.log) to every printed digit.**
Astra's bound: C > 0 ⇔ Ω > 2π·e^{γ+1/Ω} ≈ 470.07. That is correct arithmetic, but the condition is necessary for the tail and **not sufficient** for T ≻ 0 (Ω 488/536 failed).
Mori's printed formula for e is garbled in the text layer ("u₀E₀ + u₁E₁"). My derivation gives ‖(I−P)UP‖ ≤ √(Ω/π) r_N Σ u_j √E_j, and only the √E_j form reproduces Mori VII's printed e < 2.08e−94 (I get 2.0757e−94, against 1.985e−94 for the literal reading). The team uses √E_j ✓.

## 3. ‖C_p‖ ≤ s = 763/200 = 3.815: independent check (`replay/cp_norm_check.py`)
- My own weighted-Schur bound (step weight from power iteration, with Arb-enclosed shift windows widened by one cell each side, 1_I ≤ 1, final sums in Arb on exact dyadic weights):
  - 2·10⁴ cells: 3.81980 (> s, too coarse);
  - 2·10⁵ cells: 3.81261;
  - **10⁶ cells: ‖C_p‖ ≤ 3.810860 < 3.815** ✓.
- The team's S65 value of 3.81418 is consistent, and lies between Galerkin λ_max ≈ 3.8056 and my bound.
- Status: CERTIFIED (Arb on the top 2000 cells; the remaining cells by a float bound with a relative margin of 10⁻¹⁰ ≫ rounding).
- The margin s − ‖C_p‖ is thin (0.004). It is the single most delicate input, and it is now independently confirmed.

## 4. Gaps, ranked by severity
1. **[Medium, rigour] The J quadrature in the team code has no error bound.** The band integrals (1/π)∫₀^Ω (a−γ)b_ib_j dt are computed by Gauss on panels (png 96, width 16) in s67/mori_cert.py, and the ball radii cover only node rounding. No Bernstein-ellipse remainder (Mori VII Lemma 5.1 analogue) is added. The S66 README already says so, and nothing in S67–S74 adds it.
   - *Assessment:* my rigorous panel bounds (below) show the true quadrature error is astronomically small (my eps_J ~1e−148 at log 12 with 96–192 nodes per panel). The gap is very likely harmless, but as written **the team certificate is incomplete**.
   - *Fix:* add the term, or adopt my replay's bound. Note that ψ(¼ ± iz/2) has poles at z = ±i(2k+½), so the ellipses near t = 0 must be restricted. Note also that Arb's acb digamma returns NaN on wide balls with Re w < 0 at high precision, so shift with ψ(w) = ψ(w+K) − Σ1/(w+k).
2. **[Medium, bookkeeping] The final λ bound ignores the T-ball radius and the inverse error.**
   - s37_ldlt_exact.py certifies λ_min of the *midpoint* matrix Atil (sermid strings), plus the string-rounding error.
   - The assembly's radT_rowsum (3.83e−145) and ‖W⁻¹ − Y‖ ≤ 2.05e−136 are logged but never subtracted from the stated bound. Nor is the (missing) J quadrature term.
   - Numerically this is negligible against 1e−81, but the stated bound must read λ_min(T) ≥ 1.0711e−81 − 2.05e−136 − 3.83e−145 − ε_J.
3. **[Medium, wording] The quantity labelled λ_min is the margin of T, not the Weil-form constant** (§1). The window theorem is Q ≥ 2^−4200.4‖f‖² (even and odd use the same δ only if both parities run with δ = 2^−4200; odd at Ω 1150 used δ = 2^−4900, so the combined window constant is min = C_odd·2^−4900).
4. **[Low–medium] The control is not a certificate.** The "1.005× fails" control is a *non-rigorous* fixed-point LDLᵀ hitting a nonpositive pivot. That is fine as a sanity check, but it does not prove λ_min(T) < 1.005·ℓ. My replay makes the control rigorous: a negative diagonal entry of Z(T − ℓ_c)Zᵀ in balls gives an explicit vector with vᵀ(T − ℓ_c)v < 0.
5. **[Low] Prop. 3.1 and Lemma 3.2 are inherited from Mori,** and Mori VII itself cites [14] (log 10 paper) for e₀ ≤ E₀, e₁ ≤ E₁ and a(t) ≥ log(t/2π) − 1/t.
   - The last bound can be avoided: a is increasing on t ≥ 0, since each term of ψ's series is monotone. So C' = a(Ω) − γ, computed in Arb, already works; I checked C' ≥ C.
   - Theorem N (prolate defect) is Mori's own unrefereed result. Its proof (Slepian's identity, the Bessel comparison) reads correctly to me, but it is not independently refereed.
   - Hypotheses checked: C > 0, c > 2, Ω ≥ 15/4, δ ≤ d(c), u_jE_j ≤ C, n < m; all ✓ (consts.json).
6. **[Low] Single-implementation pieces.**
   - The pole vectors use sph_i from s27 (series). Mine use Arb Bessel I. They agree implicitly through the Mori VII reproduction.
   - The breakpoint logic (|log n − L|) is fine for rational L: all breakpoints are distinct by Lindemann (e^{4L} is transcendental). For L = ½log Q the coincident breakpoints (e.g. 3·4 = 12) must be deduplicated exactly. My code does this in q-space.
7. **[Low] The pre-S67 history is clean.** The S66 "λ_min" values were inverse-iteration eigenvalues nearest zero and were retracted in S67. The 2.84 claim rests only on S71/S73 dumps with an exact LDLᵀ. ✓
8. **[Structural, not a gap] Prior art and scope.**
   - Published record: Mori VIII, log 15 ≈ 2.708 (Zenodo 2026-10-10).
   - Also found: a 2026-10-06 Zenodo report claiming "Arb-certified lower bounds for ζ on windows of support 2.38, 2.6 and 2.99 in both sectors" (doi:10.5281/zenodo.23174403). **Support 2.99 > 2.84**, so a "past log 17 record" wording is unsafe until that report is read. I have not verified its method.
   - Mori VII Theorem K gives λ_min ≤ 7.58e24 μ⁵(log μ)³ e^{−4πμ} (μ = e^{2L}), so each window extension costs exponentially more.

## 5. Replay (independent code: `replay/mori_replay.py`; it does not import or copy s66/s67/s71/s37)
Method differences: node-vectorised Legendre and spherical-Bessel recurrences in Arb at boosted precision on exact points; a Markov/derivative bound for the node radii; **adaptive panels with a rigorous Gauss–Legendre remainder** (Trefethen ATAP Thm 19.3; |b_ib_j| ≤ 2L e^{2L|Im z|} by Cauchy–Schwarz; |g| bounded by Arb on ball covers of each ellipse); certificate by **unit-lower-triangular congruence Z (recursive block LDLᵀ at mid precision) plus a ball product S = Z(T−ℓI)Zᵀ and diagonal dominance** (rigorous PD); control by a rigorous negative diagonal entry. Checkpoints are per cell and per panel group (resumable).

| case | Ω | modes | NUMERICAL λ_min(T) (inverse iteration) | certified | control | extra budget (ε_J + ε_inv) | runtime |
|---|---|---|---|---|---|---|---|
| Mori VII log 12, even (published 1.2092639998e−54) | 416 | 416 | 1.20926402186e−54 | **T − 1.2092639998e−54·I ≻ 0 CERTIFIED** | at 1.2153e−54 (1.005×): **certified negative** (pivot 414) | 8.3e−149 + 1.3e−112 | 6 min (3 workers) |
| L = 71/50 even | 1000 | 1188 | RUNNING (see §6) | | | | |
| L = 71/50 odd | 1000 | 1188 | queued (control, should fail) | | | | |
| L = 71/50 odd | 1150 | 1366 | queued | | | | |

So my code path reproduces Mori's published log-12 even certificate (the same digits as Mori's Table 2 lower bound, and the same estimate as the team's S66 to 12 digits) from a code base written only from the paper.

## 6. Replay status of the 2.84 sectors
Run dir /workspace/sm/runs/, chain `chain_284.sh` (even Ω1000 → odd Ω1000 → odd Ω1150), logs `L7150_O*_p*.log`. See the latest commit for results; this table is updated as sectors finish.

## 7. Fakes battery: what the 2.84 certificate says and does not say
A window certificate is a statement about ζ's own Weil form on L²[−1.42, 1.42]. It separates ζ from a fake only if that fake's window form is *not* PSD at the same L.
- **Davenport–Heilbronn (DH):** the same functional equation, no Euler product. Its window form is NUMERICALLY PSD below L ≈ 1.7105 and CERTIFIED indefinite at L ≥ 1.7105 (S75). **So support 2.84 (L = 1.42 < 1.7105) does NOT separate ζ from DH.** The same holds for every window the box can certify (Mori route at L ≈ 1.71 needs Ω ≈ 6–7k; S75/S78).
- **F_θ = ξ(s−θ)ξ(s+θ):** its window forms fail already at L ≈ 0.5 for θ = 0.1 (S96 Galerkin onset 2L ≈ 1.0), with certified window negativity for F_{1/4} and F_{1/10} (orientation brief). The 2.84 certificate is consistent with ζ ≠ F_θ, but this separation was already available at much smaller support. Nothing new.
- **Backward heat H₋ₜ (t < 0):** S64/S65 found it NUMERICALLY failing the window by L ≈ 1.2 for t ≥ 0.5, and below the first collision t_c ≈ −0.385 within support 2.8. Support 2.84 adds nothing new beyond 2.8.
- **DMV (Diamond–Montgomery–Vorhauer Beurling system: Λ ≥ 0, a pole, RH false):** there is no window computation of the actual DMV primes in the repos (only S83's zero-side toy). The certificate says nothing about DMV.
- **G_θ′:** I could not locate a definition in GPPDiscovery5 or GPP-bridge (rg over bridge *.md). No statement.
- Net: the 2.84 certificate is a **quantitative window-positivity result for ζ**. It rules out no fake that was not already ruled out at support ≤ 2.8, and it does not rule out DH.

## 8. Dead ends and lessons logged during the audit
- Arb's acb digamma returns NaN on wide complex balls with Re w < 0 at high precision; use the shift recurrence.
- Ball arithmetic on the Legendre recurrence needs ≈ N·log₂(|y|+√(1+y²)) extra bits (≈ 1.27N at |y| = 1), as the team found in S67. Evaluating at exact midpoints plus the Markov bound is the right pattern.
- A coarse Schur grid (2·10⁴ cells) gives 3.8198 > s. The bound only drops below 3.815 at ≳ 10⁵–10⁶ cells, so s = 763/200 has little slack.
