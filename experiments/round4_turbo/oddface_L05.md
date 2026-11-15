# OddFace at the single width L = 0.5 (Turbo 6889159916 item 3)

Heart (Grok / GPPDiscovery5), Round 4. Kernels exactly as in-tree (KreinHPsi.lean, WeilGramSemilocal.lean, WeilGramParitySplit.lean): archTerm, primeTermS(S_L), bigPsiI2 = Ψ(i/2) = 16 sinh²(1/4). Script `oddface_L05.py`, results `oddface_L05_results.json` (archTerm by 96-node Gauss–Legendre per 2π period, t ≤ 3200π, dps 30). No zero data anywhere. Outcome: **a proof skeleton with named premises plus a named obstruction. OddFace(0.5) itself stays OPEN.**

Notation: ℓ = L/2 = 1/4, K_S = poleTerm − primeTerm_{S_L} + archTerm, w(t) = Re ψ_Γ(1/4 + it/2) − log π, Ψ(t)² = (sin(t/2)/(t/2))⁴ (the in-tree `bigPsi t ^ 2`), P_S(t) = Σ_{n: all prime factors in S} Λ(n) n^{−1/2} cos(t log n), V_S = w − 2P_S.

## OF0. Exact reformulations (CHECKED by hand)
- kform(K_S, reflX x, oddC c) = Σ_{i,j} c_i c_j G(x_i, x_j), with G(s,u) := ½[K_S(s−u) − K_S(s+u)]. Replacing (x_i, c_i) by (−x_i, −c_i) leaves the odd measure unchanged, so **OddFace L ⟺ G is positive semidefinite on [0, L/2]** (every finite point set). (C1's 6-point D-eigenvalues are 2 × the eigenvalues of G; the run reproduces C1's 6.527e-7 as 2 × 3.2636e-7.)
- Spectral form: g = ψ⋆ψ has ĝ = Ψ², so primeTerm_S(d) = (2/π)∫₀^∞ Ψ² cos(td) P_S(t) dt (the extra prime powers of S beyond e^{2+L} contribute 0 at |d| ≤ L). With S_c(t) := Σ c_i sin(t x_i):
  **OddFace L ⟺ (1/π)∫₀^∞ Ψ(t)² V_{S_L}(t) S_c(t)² dt ≥ 2·|Ψ(i/2)·S_c(i/2)|²** for every odd exponential sum with frequencies x_i ∈ [0, L/2] (Σ c_i sinh(x_i/2) = −i S_c(i/2)). The pole burden is the evaluation of the same sine sum at t = i/2.

## OF1. Small-scale structure; what C1's margin is (CHECKED by hand; NUMERICAL)
- Near d = 0, primeTerm_S is a cubic polynomial (nearest kink |d| = 2 − log 7 = 0.0541) and poleTerm is entire. archTerm has a |d|³ log|d| singularity at d = 0 (from Ψ² ≍ 6 t⁻⁴ on average times w(t) ≍ log t), and it is the only non-smooth term at 0. It is positive at high frequency.
- Hence for points at scale h: c^T G c = a·m₁² + O(h³ log(1/h))·|c|², m₁ = Σ c_i x_i, a := −K_S''(0).
- NUMERICAL: a = −K_S''(0) = **0.24538** = 0.76660 (arch − prime) − 0.52123 (pole Ψ(i/2)²/2). The pole burden takes 68% of the m₁ coefficient. [archTerm''(0) = 1.52192, primeTerm''(0) = 2.28852.]
- NUMERICAL scaling of C1's configuration x = (k − 2.5)h: λ_min(G)/(h³ log(1/h)) = 0.0735, 0.0709, 0.0672, 0.0710 at h = 0.005, 0.01, 0.02, 0.04 (0.0356 at h = 0.08, which is outside the asymptotic regime). The largest eigenvalue tends to a·Σy_i² (the m₁ direction).
- **Reading:** C1's 6.53e-7 is the h³ log(1/h) scale of a clustered configuration (h = 0.01), not a near-violation. The smallest direction is a higher-moment direction controlled by the archimedean singularity, not the pole direction (consistent with C1's 8e-9 pole shift).
- Spread grids on [0, 1/4] (x_i = i/(4N)): all eigenvalues of G are positive; λ_min = 4.01e-6 (N = 8) and 1.88e-6 (N = 12). This is consistent with OddFace(0.5). It is not a proof: G is compact, so λ_min → 0 as N grows.

## OB1. Named obstruction for a pointwise-symbol proof (CHECKED)
The naive symbol Ψ²V_S is pointwise negative far out. P_S(0) = Σ_{p≤11} log p/(√p − 1) = 6.6937, and by Kronecker (1, log 2, …, log 11 are Q-independent) P_S(t) returns arbitrarily close to 6.6937 at arbitrarily large t. So V_S < 0 recurs until w(t) > 13.387, i.e. t ≈ e^{15.23} ≈ 4·10⁶. No proof can bound V_S below pointwise in t. It must use that only the restriction of the kernel to |d| ≤ L matters.

## OF2. Window-restricted symbol and crossover (CHECKED by hand; NUMERICAL constants)
On |d| ≤ 0.5, primeTerm_S is piecewise cubic with 16 kinks d_k = ±(log n − 1) (n = 2, 3, 4; coefficient −(1/3)Λ(n)/√n) and ±(log n − 2) (n = 5, 7, 8, 9, 11; coefficient +(1/12)Λ(n)/√n). These come from triG's third-derivative jumps at |y| = 1, 2. The total |coefficient| is Σ|a_k| = 1.44558. Using |x|³ ↔ 12 t⁻⁴, the window kernel's high-frequency symbol is (1/(π t⁴))·[6 w(t) − 12 Σ_k a_k cos(t d_k)] plus faster-decaying terms (smooth parts and the singularities of archTerm at d = ±1, ±2, which lie outside the window). This is guaranteed positive once w(t) > 2Σ|a_k| = 2.891, i.e. **t > T_* = 113.19**. That is about ℓT_*/π ≈ 9 odd modes on [0, 1/4] below the crossover (compared with 4·10⁶ for the naive symbol).

## Proof skeleton for OddFace(0.5) (logic CHECKED; premises OPEN)
Split ν = e + r, with e in the span E_N of the first N ≈ 12 odd modes on [0, 1/4] and r in its complement.
- **P1 (finite block, OPEN):** the form on E_N is ≥ β‖e‖², with β > 0 certified in interval arithmetic (a #597-style LDLᵀ certificate on the 12-dimensional block).
- **P2 (high-mode coercivity, OPEN):** on E_N^⊥, the form is ≥ κ·E(r), where E is the window-symbol energy above T_*. This includes explicit bounds on the faster-decaying remainder symbol between the modes of E_N^⊥ and T_*.
- **P3 (cross term, OPEN):** |B(e, r)|² ≤ (1 − η)·β‖e‖²·κE(r) for some η > 0.
- P1 ∧ P2 ∧ P3 ⇒ OddFace(0.5) by the 2 × 2 Schur argument (CHECKED).

## Scope (CHECKED / CONJECTURE as marked)
- OddFace(0.5) follows from RH (Weil positivity), so I expect no counterexample. The obstruction is to methods, not to truth.
- **Uniformity in L (CHECKED as a method bound):** the kink mass in |d| ≤ L is ≍ Σ_{n ≤ e^{2+L}} Λ(n)/√n ≍ 2e^{1+L/2} (PNT). So this method's crossover is T_*(L) ≈ 2π·exp(c·e^{L/2}), and the P1 block size grows like L·T_*(L). A single-L certificate does not extend to ∀L by this route. ∀L needs a structural argument, not more certificates.

Labels: OF0 CHECKED. OF1 CHECKED (structure) + NUMERICAL (a, scaling, grids). OB1 CHECKED. OF2 CHECKED + NUMERICAL. P1–P3 OPEN. OddFace(0.5), OddFace ∀L, EvenFace, H_ψ, WeilIsGram, ARP, RH: OPEN.
