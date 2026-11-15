# Round 6 (Turbo 6889647810): P2/P3 at L = 0.5, N = 12, β-normalisation (L² on [0, 1/4])

Heart (Grok / GPPDiscovery5). Scripts: `p2_boundary_float.py`, `p2_schur_sweep.py`, `p2_remainder_probe.py` (float64 spectral quadrature, same exact reduction as Round 5; it reproduces the certified Round-5 block to 2.8e-11 relative), and `p2_bordered_arb.py` (arb, 160-bit, interval-certified). No zero data. Outcome: **the P2/P3 statements in my Round-4 skeleton are wrong as stated. I give the corrected statements, a certified leading-order Schur condition, and a named remainder. P2 and P3 stay OPEN.**

## R6.1 Exact boundary expansion (CHECKED by hand)
G(s,u) = ½[K_S(s−u) − K_S(s+u)], with K_S ∈ C² and K_S''' locally integrable (arch has a |d|³log|d| singularity, prime is triG piecewise cubic, pole is smooth). Integrate by parts three times in s, using G(0,u) = 0, ∂²_sG(0,u) = 0 (K even) and sin(ω_j ℓ) = 0:
  (Gφ_j)(u) = √(2/ℓ)(−1)^{j+1}[G(ℓ,u)/ω_j − ∂²_sG(ℓ,u)/ω_j³] − √(2/ℓ)ω_j^{−3}∫₀^ℓ ∂³_sG(s,u)cos(ω_j s)ds.
Hence, with w_j := √(2/ℓ)(−1)^{j+1}/ω_j, g := G(ℓ,ℓ) = ½[K_S(0) − K_S(L)] and u_i := B(φ_i, δ_ℓ) = ∫φ_i G(·,ℓ):
  M_jk = g·w_j w_k + O(ω^{−3}·ω^{−1}) + (interior diagonal ≈ 3w(ω_j)/ω_j⁴),   M_ij = u_i w_j + O(ω_j^{−3})   (i ≤ 12 < j, k).
Spectrally, δ_ℓ has transform sin(tℓ) and pole coefficient sinh(ℓ/2). So g = (1/π)∫Ψ²V sin²(tℓ) − 2Ψ(i/2)²sinh²(ℓ/2), and similarly for u_i.

## R6.2 What this does to P2/P3 as I stated them (CHECKED + NUMERICAL)
- **g = 1.23787e-4 > 0** (NUMERICAL; the spectral and kernel-domain values agree to 2.6e-15; arb ball [1.2379e-4 ± 3.1e-9]).
- **P2 as framed is wrong.** E_12^⊥ is not "frequency above T_*". Every high mode has boundary content g·w_j² = (2/ℓ)g/ω_j², which dominates its window-symbol energy ≈ 3w(ω_j)/ω_j⁴ by a factor of ~ω_j². NUMERICAL: M_jk/(g w_j w_k) = 1.0038 at (40,41), 1.0024 at (60,90), 1.014 at (90,90). The diagonal excess is ≈ 2.1–2.7·log ω/ω⁴, which matches the window-symbol order. Most of the high block is the rank-one boundary term, so "window symbol above T_*" alone cannot be P2.
- **P3 against β‖e‖² is false.** λ_max(M12 M22^{−1} M21) = 1.98e-5 (K = 400 high modes), versus β = 1.66823e-8. The cross-term loss is ≈ uuᵀ/g (|u|²/g = 1.95e-5), and it lives in directions where M11 is also large. The scalar form |B(e,r)|² ≤ (1−η)β‖e‖²κE(r) fails by a factor of ~10³.
- **Correct P3 (strengthened Cauchy–Schwarz):** |B(e,r)|² ≤ γ² B(e,e) B(r,r) with γ < 1. NUMERICAL: γ² = λ_max(M11^{−1}M12M22^{−1}M21) = 0.9667 (K = 400); the leading order alone gives 0.9644. Equivalently, **the Schur complement M11 − M12M22^{−1}M21 ≻ 0**. P1 ∧ P2(M22 ≻ 0) ∧ P3(γ < 1) ⇒ B ≥ (1−γ)(B(e,e) + B(r,r)) ≥ 0 (CHECKED). β still enters, because P1 certifies M11 and the Schur margin is measured in the same L² normalisation.

## R6.3 Certified leading-order Schur condition (NUMERICAL, interval-certified; not Lean)
Bordered Gram of the form on {φ_1..φ_12, δ_ℓ}. **g > 0 and M11 − uuᵀ/g ⪰ β₀ I, with β₀ = 1.30163e-8**: this is the LDLᵀ of M11 − uuᵀ/g − β₀I in ball arithmetic, with all 12 pivot lower bounds positive (1.256e-5 down to 1.496e-9). M11 is the banked Round-5 ball block. g and u_i come from acb.integral on [0, 4000] plus a hand tail bound (CHECKED: Ψ² ≤ 16/t⁴, |V| ≤ log t + C, |S_i| ≤ (4/3)ω_i/t²). This is the exact limit of the Schur loss along the boundary direction (wᵀ(g wwᵀ + H)^{−1}w ↑ 1/g).

## R6.4 Full truncated Schur complement (NUMERICAL, float)
λ_min(M11 − M12M22^{−1}M21) with modes 13..12+K: 1.4710e-8 (K=20), 1.4341e-8 (50), 1.4234e-8 (100), 1.4190e-8 (200), 1.4171e-8 (400). The increments shrink by a factor of ~0.4 per doubling, so the limit is ≈ 1.416e-8, above β₀. The beyond-leading-order loss is ≈ 3e-10, which is ~2% of the certified margin β₀.

## R6.5 KILLED: subtract-the-rank-one route
H := M22 − g wwᵀ is **indefinite**: λ_min(D^{−1/2}HD^{−1/2}) = −0.53, −4.5, −17, −56 at K = 50, 100, 200, 400, with D_j = 3w(ω_j)/ω_j⁴. The cause is the second boundary term −∂²_sG(ℓ,·)/ω³. In D-units it couples mode j to mode k with weight ~ω_k/ω_j. So "M22 = g wwᵀ + (positive remainder) and bound E by H^{−1}" is wrong. The boundary block cannot be removed as a rank-one term with a coercive remainder.

## Named remainder (the exact next work; GAP)
Let Y = E_{12+K} for a fixed K (e.g. K = 60) and Z = its complement.
- **(a) Finite part:** certify by arb LDLᵀ that the bordered Gram on Y ⊕ {δ_ℓ} has Schur complement onto E_12 ⪰ β₁ I (same β-normalisation).
- **(b) Tail part:** for j > 12+K, bound the residual coupling ‖(Gφ_j) − w_j G(·,ℓ)‖ by the explicit ω^{−3} terms of R6.1, namely the sup of ∂²_sG(ℓ,·), and ∫|∂³_sG| (log-integrable arch part plus kink jumps a_k) of R6.1. Then show that the extra Schur loss from Z is < β₁. The float data predict that the extra loss from modes > 12+K is ≲ 1e-10 for K ≥ 50, against a margin of ~1.3e-8.
- Gate discipline: no pointwise symbol (OB1), no prime-sum J-Gram (J5). EvenFace stays deferred.

Labels: R6.1 CHECKED. g > 0, ratios, γ² and sweep: NUMERICAL. R6.3 leading-order bordered Schur (β₀ = 1.30163e-8): NUMERICAL, interval-certified, not Lean. Scalar P3 against β: KILLED (fails by ~10³). Rank-one-subtraction route: KILLED. Corrected P2 (M22 ≻ 0) and P3 (γ < 1 / Schur ≻ 0) with explicit tail: OPEN (GAP named above). OddFace(0.5), H_ψ, WeilIsGram, ARP, RH: OPEN.
Credit: Heart (Grok / GPPDiscovery5); audit chain Turbo (6889647810); skeleton with Flint per PR #11.
