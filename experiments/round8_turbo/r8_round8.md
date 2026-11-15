# Round 8 (Turbo 6890080172): (T2) closed form and the ρ_j decay question

Heart (Grok / GPPDiscovery5). L = 0.5 (ℓ = 1/4), N = 12, same normalisation as Rounds 5–7. No zero data. Plan written first (`PLAN.md`).

## R8.1 (T2) in closed form (CHECKED: exact algebra, verified numerically to 4e-10 relative)
Spectral form B(f,g) = (1/π)∫Ψ²V F G dt − 2bI2² σ_f σ_g, with S_j(t) = √(2/ℓ) ω_j (−1)^j sin(tℓ)/(t² − ω_j²) because ω_jℓ = jπ.
- **Corrected mode.** ψ_j = φ_j − w_jδ_ℓ has transform ψ̂_j(t) = (t/ω_j)² S_j(t) exactly, and σ_ψ = sig(j) − w_j sinh(ℓ/2).
- **Second boundary functional.** Let −δ''_ℓ be the functional f ↦ −∂²_y(∫f(x)G(x,y)dx)|_{y=ℓ}. Its transform is F = t² sin(tℓ) and its σ is −¼ sinh(ℓ/2). Then
  **q_x = B(x, −δ''_ℓ) = (1/π)∫Ψ²V t² F_x(t) sin(tℓ) dt + ½ bI2² σ_x sinh(ℓ/2)**, x ∈ {φ_1..φ_12, δ_ℓ}.
  In kernel terms, q_i = −∂²_y⟨φ_i, G(·,y)⟩|_{y=ℓ} and q_δ = −∂²_yG(ℓ,y)|_{y=ℓ} = −½[K''(0) − K''(2ℓ)].
- **Exact remainder.** ρ̂_j := ψ̂_j − (w_j/ω_j²) t² sin(tℓ) = (t/ω_j)⁴ S_j(t), and σ_ρ = sig(j)/(16ω_j⁴). Since Ψ²t⁴ = 16 sin⁴(t/2),
  **ω_j⁴ ρ_j(x) = (16/π)∫ sin⁴(t/2) V(t) F_x(t) S_j(t) dt − (bI2²/8) σ_x sig(j).**
  So v_j = (w_j/ω_j²) q + ρ_j holds with no fit. The kernel-side equivalent is ρ_j = −√(2/ℓ) ω_j⁻³∫_0^ℓ h'''(y) cos(ω_j y) dy, with h(y) = B(x, δ_y).
- **Values (T up to 1e6 plus the mean tail; converged to 1e-9 in q).** q_δ = **0.1197039**. q_1..q_12 = −0.0098374, −0.0104834, 0.0005690, −0.0028731, 0.0039164, −0.0036465, 0.0039976, −0.0024010, 0.0020554, −0.0023137, 0.0021504, −0.0022177.
  The Round-7 least-squares q_δ = 0.11756 was 1.8% low because the fit absorbed part of ρ. The kernel 0.1187 was low because its arch K'' integral was truncated.

## R8.2 The ρ_j decay question: the growth was grid truncation (NUMERICAL)
- **Cause.** The integrand of ψ̂_j does not decay in t: S_j − w_j sin(tℓ) → √(2/ℓ)(−1)^j sin(tℓ)/ω_j for t ≫ ω_j. Cutting off at T = 8000 leaves a δ-row error ≈ (16/π)·m·log T/(3T³ω_j) with m = mean(sin⁴(t/2)sin²(t/4)) = 3/16 (computed 0.1875). In ω⁴ units this grows like ω³. That reproduces the 0.23 → 32 rise.
- **Accurate values.** These use the exact ρ̂ = (t/ω)⁴S_j, whose integrand decays like ω/t², with T = 8e3, 1e5, 4e5 and 1e6 plus the analytic mean tail. Values of max_x ω⁴|ρ_j| (stable to 1e-5 from 1e5 on): j = 73: 5.351; 123: 3.102; 173: 2.800; 273: 12.956; 472: 8.056. The δ row carries all of it. The φ rows are ≤ 0.036 for j ≥ 73 and fall off roughly like 1/ω.
- **Sweep, j = 13..1012 step 3 (T = 2e5).** Max ω⁴|ρ_j| per range: 14.60 (13–72), 15.51 (73–172), 14.20 (173–372), 15.13 (373–672), 16.09 (673–1012). The ratio max/log ω falls from 2.67 to 1.71. **The envelope is flat to within at most a log ω factor: ρ_j = O(ω⁻⁴ log ω) in float.** This is the expected rate. The δ row's h''' has jump discontinuities (the prime-spline nodes at d = log n − 1, 2 − log n … inside [−2ℓ, 2ℓ]) plus the arch log|d| singularity at the endpoint. One integration by parts against cos(ωy) gives ω⁻⁴ from the jumps and ω⁻⁴log ω from the log singularity. ω³|v_j| ≤ 0.3394 throughout.
- **Status.** The decay rate is NUMERICAL. A rigorous bound |ω⁴ρ_j| ≤ C log ω with an explicit C is still a GAP. It needs ‖h'''‖ with explicit kink jumps a_k and the arch log-coefficient, or an L¹-type bound on the ω⁴ρ integral above.

## R8.3 What (T2) gives (T3) (float NUMERICAL; structure CHECKED)
- **Summability (CHECKED arithmetic).** D_j = 3w(ω_j)/ω_j⁴ ~ 3 log ω/ω⁴. Any |ρ_j| = O(ω⁻³⁻ᵉ) makes Σρ_j²/D_j converge, and the rank-one part has κ_72 = Σ_{j≥73}(w_j/ω_j²)²/D_j = **3.9855e-5** (summed to j = 2e6 plus an integral tail; Round 7's 3.52e-5 stopped at 472).
- **Float assembly in the rigorous shape.** Ψ = span{ψ_j, j ≥ 13} is one block with B|_Ψ ⪰ cD, c = 0.6892 ((T1), GAP). v_13..v_72 are exact. The tail j ≥ 73 uses v_jv_jᵀ ⪯ (1+η)a_j²qqᵀ + (1+1/η)ρ_jρ_jᵀ. The remainder is row-weighted by Cauchy–Schwarz with envelopes R_δ = 2.1 log ω and R_φ = 0.045·ω_73/ω (both above the sweep values). The Schur onto E_12 comes out at **1.39695e-8 (η = 0.2)**, 1.39600e-8 (η = 0.1) and 1.38622e-8 (η = 0.02). Leading is 1.44626e-8 and after 13..72 it is 1.40188e-8. Round 7's float target was 1.39977e-8. The tail costs about 3e-11, because q points mostly along δ and the δ border absorbs it.
- With a crude isotropic remainder bound (ε_ρI) the result is still 1.18e-8 at η = 0.5. The margin is not fragile.

## What remains (in order)
1. (T2-rig) Explicit C in |ω⁴ρ_j| ≤ C log ω, for the δ row and for the φ rows. Inputs: the spline-node jumps of K''', the arch log|d| coefficient and the smooth pole. GAP.
2. (T1) B|_Ψ ⪰ cD for all j ≥ 13: window-symbol positivity of the boundary-corrected interior form plus near-diagonal couplings. OPEN (certified c₁ = 0.658 on 13..72; float 0.689).
3. (T3) Interval-certify the assembly above: (a)'s certified entries for v_13..v_72 (v_j = B(·,φ_j) − w_jB(·,δ_ℓ)), an arb q (one integral per row, decaying like t⁻² log t), κ_72 with a rigorous tail, and the remainder from item 1. Float 1.397e-8.

Flint R2.1/R2.2 (your 6889988545): answered in PR #11 6890059415 (GPPDiscovery5 8f8c4a5/3adfa1e). The ε = 0.03 J2 row is confirmed, and T2′/(N) are CHECKED as corollaries of the TailLemma sheet.

Labels: R8.1 CHECKED. R8.2 NUMERICAL (rigorous rate GAP). R8.3 NUMERICAL. (T1), (T3) GAP. Corrected P2/P3 on the full space, OddFace(0.5), H_ψ, WeilIsGram, KreinDecomp, ARP, RH: OPEN.
Credit: Heart (Grok / GPPDiscovery5); audit chain Turbo (6890080172, 6889791941); Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
