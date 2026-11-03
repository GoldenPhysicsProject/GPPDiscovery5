# C.16 odd (K2) certificate: matrix spec for the Arb run (S&M)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.
Theorem: ARP §C.16, Thm C.16.A (odd analogue of Thm C.15.A) + Prop C.15.B (monotone in L). Float values are NUMERICAL
(arp_c/c16_schur_general.out, c16_odd_log3.out). Callbacks: arp_c/c16_arb_callbacks.py (class Callbacks, parity='odd').

## Recommended single run (covers ALL of 2L ∈ (log 2, log 3], i.e. ℓ ∈ (0, log 1.5])
| run | 2L | prime set S | T | N | majorants | float μ_A | float s_V | float ε^abs | float ε^neg | float λ_cert |
|---|---|---|---|---|---|---|---|---|---|---|
| **O3** (recommended) | log 3 | {2} | 60 | 128 | analytic | 8.050e-4 | 3.3147 | 1.639e-5 | 4.092e-6 | **+7.63e-4** |
| O3' (sharper) | log 3 | {2} | 100 | 256 | analytic | 8.85e-4 | 3.8267 | 9.38e-6 | 2.39e-6 | +8.65e-4 |
| O1 (as asked) | log 2 + 0.1 | {2} | 60 | 64 | analytic | 0.217254 | 2.7336 | 5.24e-5 | 1.22e-5 | +0.217119 |
| O2 (as asked) | log 2 + 0.2 | {2} | 60 | 64 | analytic | 0.063857 | 2.7596 | 7.06e-5 | 1.75e-5 | +0.063696 |

O3 alone implies O1, O2 by Prop C.15.B (odd version). At 2L = log 3 the n = 3 line contributes exactly 0 (g(±2L) = 0 for
ξ ∈ H¹₀), so S = {2} there. Margin check: λ_cert(O3)/Galerkin = 7.63e-4/9.52e-4; Schur loss 4.3e-5, comparison loss 1.46e-4.

## Objects (odd sector)
- L = (log 2 + ℓ)/2 (O3: L = ½ log 3). c₂ = log 2/√2. a = log 2.
- Modes e_k(x) = L^{-1/2} sin(ω_k x), ω_k = kπ/L, k = 1..N (code index k−1). f_k = e_k/n_k, n_k = (ω_k² + ¼)^{1/2}; H-orthonormal.
- ê_k(t) = √L[sinc((ω_k − t)L) − sinc((ω_k + t)L)] (real form; |ξ̂|² uses ê²). Entire.
- Ω(t) = ½[ψ(¼ + it/2) + ψ(¼ − it/2)] − log π (meromorphic, poles at |Im t| ≥ ½). q(t) = Ω(t) − 2c₂cos(at). q_T = Ω(T) − 2c₂.
- A_jk = q_T δ_jk + (1/π)∫_0^T (¼ + t²)(q(t) − q_T) ê_j ê_k dt/(n_j n_k).
- M_jk = (1/π)∫_0^T (¼ + t²) m_abs(t) ê_j ê_k dt/(n_j n_k), m_abs = Ω(T) − Ω(t) + 2c₂(1 − cos at).
- ε^g = (1/π)∫_0^T (¼ + t²) g(t) K_W(t) dt, g ∈ {m_abs, m_neg = Ω(T) − Ω(t)}, K_W = K_o − Σ_{k≤N} ê_k²/n_k²,
  K_o(t) = [L − L·sinc(2tL) − 2 sin(tL)(½coth(L/2)sin(tL) − t cos(tL))/(t² + ¼)]/(¼ + t²) (Lemma C.16.K).
  Alternative (no cancellation): mode sum k = N+1..K plus remainder ≤ sup g·(T/4 + T³/3)/π·(4/L)·(L/π)/(3(ω_K − π/L − T)³).
- λ_cert = ½[μ + c − √((μ − c)² + 4 s ε^abs)], c = q_T − ε^neg (upper bound ε^neg), μ ≤ λ_min(A), s ≥ λ_max(M).

## Required enclosures (REMARK: hypotheses of Thm C.16.D)
1. Every A_jk, M_jk as an arb ball (acb.integral on [0, T], panels of width 5; prec ≥ 128). Symmetry: compute j ≤ k only.
2. μ: interval Cholesky (or LDLᵀ) of A − μI ≻ 0 with all entries balls; suggested μ = 0.999·(float μ_A) = 8.04e-4 (O3).
   An enclosed float eigenvalue is NOT acceptable.
3. s: Cholesky of sI − M ≻ 0 (suggested s = 3.32 for O3), or s := enclosed ‖M‖_F (float 6.100; gives λ_cert = 7.27e-4, still > 0).
4. ε^abs (upper), ε^neg (upper), q_T (enclosed). Then λ_cert as an arb ball; success iff its lower endpoint > 0.
5. Nothing is integrated on |t| ≥ T (step 1 uses q − q_T ≥ 0 there: Ω increasing, PROVED).
Cost estimate: O3 has 8256 entries per matrix, ~1–2 s each at 128 bits (self-test timing), i.e. ~5–9 CPU-hours per matrix; embarrassingly parallel.
