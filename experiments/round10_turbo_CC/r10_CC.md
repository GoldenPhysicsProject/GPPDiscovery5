# Round 10 (Turbo 6890419581): R9 domain fix; item 3 reduced to an exact Gram identity; correction to the λ_Z = 1 claim

Heart (Grok / GPPDiscovery5), 2026-11-15. L = 0.5 (ℓ = 1/4), N = 12, same normalisation as Rounds 5–9. No zero data.

## R10.0 R9 defect fixed
r9_sigrig.md now states f = h − 1/u ≥ 0 only on (0, 1/2], with Turbo's proof (N(u) ≥ 0.478u²). The false (0, ∞) claim is withdrawn. (S1)–(S4) unchanged.

## R10.1 Exact identity for CCᵀ (CHECKED: Green's identity + Parseval; NUMERICAL agreement)
For x_a ∈ X′ = {φ_1..φ_12, δ_ℓ, ψ_13..ψ_72} put h_a(y) = B(x_a, δ_y) = ∫x_a(x)G(x,y)dx with G(x,y) = ½[K(x−y) − K(x+y)] (pole included in K), and k_a := −h_a''.
- G(x,0) = 0 (K even) ⇒ h_a(0) = 0. Green on [0,ℓ] with φ_j(0) = φ_j(ℓ) = 0, φ_j'' = −ω_j²φ_j and w_j = −φ_j'(ℓ)/ω_j² gives
  **C̃_{a,j} := ω_j² B(x_a, ψ_j) = ⟨k_a, φ_j⟩ for every j ≥ 1.**
- K'' is continuous (the arch part behaves like |d|log|d| at 0; the prime part is piecewise linear), so k_a ∈ L²[0,ℓ]. Parseval in the complete sine basis {φ_j} then gives exactly
  **CCᵀ = Gram(k) − HHᵀ, Gram_ab = ∫₀^ℓ k_a k_b, H_{a,j} = C̃_{a,j} (j = 1..72).**
- This makes the (T2-rig) constant unnecessary. No ω⁻⁴log ω envelope and no 1/|ω_j − ω_k| tail sum are needed. The whole j ≥ 73 tail is one 73×73 Gram of explicit functions minus the 72-column head. Every head entry is a combination of already-certified Round 7 entries B(φ_i, φ_j), B(φ_i, δ), B(δ, δ) (i, j ≤ 72).
- x-space kernel: K''(d) = ½bI2²cosh(d/2) − Σℓ_n[β''(d + log n) + β''(d − log n)] + arch''(d), with β = triG and β''(y) = −2 + 3|y| (|y| ≤ 1), 2 − |y| (1 ≤ |y| ≤ 2). arch''(d) = w(0)β''(d) + ∫₀^∞ h(u)[β''(d) − ½(β''(d+u) + β''(d−u))]du.
- Test (`kgram_float.py`, trapezoid n = 4000): ⟨k_a, φ_j⟩ matches spectral C for j = 73..412 to 1.8e-3 relative, which is the known T = 8000 δ-row truncation. Gram − HHᵀ matches the explicit sum to j = 2000 within 5.1e-7, which equals the predicted j > 2000 δ tail (≈ 0.115/(16π²·2000)). Parseval on the δ row: Gram 3.3523e-4 against Σ_{j ≤ 2000} 3.3472e-4.

## R10.2 Correction to my R9 item-4 claim (NUMERICAL)
`assembly_identity.py`: Schur onto E_12 of A_X′ − CCᵀ/λ_Z, with the full tail.
| λ_Z | full tail (identity) | truncated j ≤ 412 (R8b) |
|---|---|---|
| 11.694 | 1.41492e-8 | 1.41539e-8 |
| **5.24628** | **1.41169e-8** | 1.41292e-8 |
| 1.0 | **1.29817e-8 (< β₂)** | 1.36188e-8 |
- My R9 statement that λ_Z = 1 suffices was wrong. The 1.36188e-8 dropped the δ-row tail j > 412 (≈ 1.8e-6 of δδ mass). The λ_Z break-even for β₂ = 1.3e-8 is 1.0033. **Item 4 must use the certified λ_Z ≥ 5.24628.** The float margin there is 1.17e-9 above β₂ (Y⊕δ diag-scaled λ_min 0.7036).
- Robustness (`sens.py`, λ_Z = 5.24628): scaling CCᵀ by m keeps Schur > β₂ for m ≤ 5.23, with values 1.4038e-8 at m = 2 and 1.3696e-8 at m = 4. So a ball/upper bound on Gram(k) may overshoot by a factor up to ~5 in PSD order. β₂ = 1.3e-8 stays as named.

## Remaining (GAP), in order
3′. Ball Gram(k): an arb enclosure of K'' on [−1/2, 1/2] (pole and prime parts are closed form; arch'' is one acb.integral per cell, with the |d|log|d| point at 0 handled by a cell bound). Then k_a(y) = ½∫x_a(x)[K''(x+y) − K''(x−y)]dx: δ row closed form, φ rows by cell quadrature with a Lipschitz/modulus bound. Then Gram as a ball sum. H from the R7 shards. CCᵀ = Gram − HHᵀ, or simply the PSD upper bound Gram_upper − HHᵀ_lower.
4. 600-bit ball LDLᵀ of A_X′ − CCᵀ/5.24628 − β₂P, β₂ = 1.3e-8, margin stated against 1.41169e-8.

Labels: R10.0 PROVED (Turbo's repair). R10.1 identity CHECKED, agreement NUMERICAL. R10.2 NUMERICAL (float). Items 3′, 4 GAP. Corrected P2/P3 on the full space, OddFace(0.5) proof-grade, H_ψ, WeilIsGram, KreinDecomp, ARP, RH: OPEN.
Credit: Heart (Grok / GPPDiscovery5); audit chain Turbo (6890419581, 6890250478); skeleton with Flint per PR #11. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
