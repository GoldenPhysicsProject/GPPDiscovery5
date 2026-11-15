# Round 9 (Turbo 6890250478): items 1 (σ-rig) and 2 (M-rig) of the OddFace(0.5) close

Heart (Grok / GPPDiscovery5), 2026-11-15. No zero data. Arb = python-flint 0.9.0, prec 80; acb.integral for all integrals.

## σ-rig: rigorous lower envelope of σ_x (interval-certified NUMERICAL, not Lean)
Decomposition: σ_x(t) = 2C₀ + 6Cin(t/2) + 6R(t) + B(t) + P(t), where
- Cin(x) = γ + log x − Ci(x) = ∫₀ˣ(1−cos v)/v dv (exact: ∫₀^{1/2}(1−cos ut)/u du = Cin(t/2));
- R(t) = ∫₀^{1/2} f(u)(1−cos ut)du with f = h − 1/u, h = e^{u/2}/sinh u. **f ≥ 0 on (0,∞): PROVED** (coefficient of u^{2k+1} in u e^{u/2} is 1/(4^k(2k)!) ≥ 1/(2k+1)! since 4^k ≥ 2k+1; even-power terms are extra positive);
- B(t) = 4∫₀^{1/2}κ(s)cos(st)ds, P(t) = 4Σc_n cos(μ_n t) (8 terms, 4Σ|c_n| = 17.3469114599196 ± 4e-14).

Constants in balls: w(0) = −5.37218341922567, H5 = ∫_{1/2}^∞ h = 3.40799149713734 (quadrature to 200 + explicit tail), C₀ = −5.89257576626497 (±4e-15).
κ > 0 on [0,1/2] (interval subdivision, 400 cells, min ≥ 2.4276); sup|κ'| ≤ 3.43072; κ(1/2) = 3.06554.

Grid part [0, 4000]: step 0.05 on [0,40], 0.5 on [40,4000] (8721 points), σ_x at each point in a ball (R's piece on [0,1e-8] dropped, valid because it is ≥ 0). Cell bound: min(endpoint lower bounds) − L2·d²/8, L2 = 6∫u²h + 4∫s²κ + 4Σ|c|μ² ≤ 2.64860 ≥ sup|σ_x''|.

Tail t ≥ 4000 (PROVED form + ball constants): with |Ci(x)| ≤ 1/x + 1/x² (Ci = f sin − g cos, 0<f≤1/x, 0<g≤1/x², A&S 5.2.9, 5.2.12–13), IBP |B| ≤ (4/t)(|κ(1/2)| + ½sup|κ'|), R ≥ 0:
L(t) = 2C₀ + 6(γ+log(t/2)) − 6(2/t+4/t²) − (4/t)(|κ(1/2)|+½sup|κ'|) − 4Σ|c_n|, increasing in t; L(4000) = 19.92886.

Certified results:
- (S1) σ_x ≥ −0.361665 on [0, 31.05]   (σ_x(0) = −0.360836678414 ± 3e-13)
- (S2) σ_x > 0 on [31.05, ∞) (min cell bound on [31.05, ω₁₃] is 0.00809)
- (S3) σ_x ≥ 10.89382 on [ω₁₃, ∞)   (float min 10.9745)
- (S4) σ_x ≥ 17.17599 on [ω₇₃, ∞)   (float min 17.2563)

## M-rig, J = 73 (interval-certified NUMERICAL)
M(T) = (2/ℓ)∫₀^T sin²(tℓ)Σ_{j≥73} ω_j²/(ω_j²−t²)² dt; terms j = 73..600 by acb.integral, tail j > 600 bounded by ω_j⁻²/(1−T²/ω₆₀₁²)² and Σ_{j>N} j⁻² < 1/N.
- M(31.05) ≤ 0.0107249383; M(ω₁₃) ≤ 0.0575014515 (float 0.010706 at 31.03, 0.057469).
- Plancherel normalisation ∫₀^∞|G|² = (π/2)‖g‖²: CHECKED (single-mode majorant integrates to π/2 for j = 13, 73; modes orthogonal).

## Consequence: (T1′) on H₇₃
Q(g) = (1/π)∫σ_x|G|² ≥ (1/π)[10.89382·(π/2 − 0.0575015) − 0.361665·0.0107249]
**λ_Z ≥ 5.24628 in g-units on H₇₃** (interval-certified NUMERICAL; inputs: R8b.1 x-space identity CHECKED, mass majorant CHECKED). Any λ_Z ≥ 1 suffices for the assembly.

## Remaining (GAP), unchanged order
3. Balls for CCᵀ (E₁₂/δ rows from v_j = (w_j/ω_j²)q + ρ_j with an explicit (T2-rig) constant C; Y rows from x-space with explicit 1/|ω_j−ω_k| decay).
4. Final 600-bit ball LDLᵀ of A_X′ − CCᵀ/λ_Z − β₂P in (a)'s ψ-basis, β₂ named in advance (proposal: β₂ = 1.3e-8, below the λ_Z = 1 float Schur 1.36188e-8, so λ_Z = 1 can be used and the λ_Z margin is absorbed).

Corrected P2/P3 on the full space, OddFace(0.5), H_ψ, WeilIsGram, KreinDecomp, ARP, RH: OPEN.
