
## 2. Fast rigorous PD pipeline (`s28_onestroke_fast.py`)
The mathematics and error budget are the same as `s27_onestroke_arb.py`. The engineering is new:
1. Assembly runs in parallel over interleaved panels (multiprocessing, fork). The Bessel ladder j_0..j_{2K} is computed once per node and shared by both parities. Partial Gram sums come back as Arb balls printed with radius (the printed ball contains the original).
2. The shift is chosen by inverse power iteration on the midpoint matrix. This is a non-rigorous guess, so no bisection is needed.
3. **One rigorous PD check by blocked Schur complements in ball arithmetic** (block 48): ball-Cholesky of A₁₁, then S = A₂₂ − A₂₁·A₁₁⁻¹A₁₂ by `arb_mat.solve` and `arb_mat` multiply (C level), then recurse. It returns True only if every symmetric matrix in the ball is PD (A ≻ 0 ⇔ A₁₁ ≻ 0 and S ≻ 0).
Speed: L = 1, T = 300, 262 modes takes 59 s wall in total (S27: about 10 min per parity). Sanity checks: L = 0.8, T = 80 lifted gives λ_est 7.78558e−18 vs S27 fail-shift 7.78558e−18. L = 1 gives 4.142300e−30 vs 4.142302e−30, and the certified values above.

## 3. CERTIFIED (Arb, computer-assisted, Lean-open; single implementation): ζ Weil positivity on support 2.38 (L = 1.19)
Setup: per-prime lifted one-stroke (`../s27/lift_L1.19.json`, comb n = 2, 3, 4, 5, 7, 8, 9). Arb κ bounds on a 40000-cell θ grid: κ₂ ≤ 1.25372, κ₃ ≤ 1.13797, κ₅ ≤ 0.75727, κ₇ ≤ 0.77381, so A ≤ 3.922743 (raw A_L = 7.075; unlifted T₁ = 7.4e3; lifted T₁ = 317).
Run parameters: T = 700, β = 0.789032, 730 Legendre modes per sector (degrees ≤ 1458/1459), 1400 panels × 60 Gauss nodes, 384-bit Arb. M_S ≤ 123.12, M₀ ≤ 21.893. Quadrature error per entry ≤ 9.30e−57; ‖A − Ã‖ ≤ 6.79e−54; ε_B ≤ 1.1e−143. Wall time 375 s on 8 cores.

| sector | λ estimate | Schur-certified shift | certified λ ≥ (Arb tails, `s28_tails_arb.py`) |
|---|---|---|---|
| even | 6.888847e−48 | 6.882e−48 | **6.8819e−48** |
| odd | 3.918448e−44 | 3.915e−44 | **3.9149e−44** |

⇒ **Q_ζ(f) ≥ 6.88·10⁻⁴⁸‖f‖² for every complex f supported in [−1.19, 1.19]**, i.e. Weil positivity for all g = f⋆f̃ with supp g ⊆ [−2.38, 2.38]. This is the support Zhu claimed and then retracted. Zhu's certified range is 1.6; S27 reached 2.0.
- Stability: the same run at 40 nodes/panel gives λ_est 6.888842e−48 (vs 6.888847e−48 at 60), but there its quadrature radius (4.6e−35/entry) swamps λ, so it is not a certificate (logged).
- Caveats: computer-assisted, single implementation (needs an independent rerun), not Lean. The cost grows doubly exponentially (§1). Finite windows do not imply RH. **RH OPEN.**
- φ does not appear anywhere in S28.
