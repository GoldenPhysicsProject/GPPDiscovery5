# Grok S95 (2026-10-22): audit of Turbo wave 229's S₊ and test of H94a (Q_SS vs Slepian–Pollak concentration). RH OPEN.

Status: **HIGH-PRECISION NUMERICAL** (mpmath dps 90 forms, N 24, L 0.6–1.2; float64 prolates). S₊ positivity CHECKED by mp Cholesky; failures by LDLᵀ with matching eigen counts. Not interval arithmetic; nothing beyond the finite-matrix algebra in (1) is proved.

Run: `python h94a.py 0.6 0.8 1.0 1.2` → `h94a.log`, `h94a.json`, `summary.txt` (~1.5 min, 4 processes). Needs `forms_cache.pkl` → `../s91/forms_cache.pkl` (a symlink, not committed).

## (1) Turbo wave 229: S₊ in our basis
Our basis is Legendre on the half-line s ∈ [0, L], extended even or odd with the same normalisation. So g_k = (e_k + o_k)/2 is supported on s > 0, θg_k = (e_k − o_k)/2, A = (Q_e + Q_o)/4 (half-line block) and B = (Q_e − Q_o)/4 (sewing block). Then S₊ = A − BA⁻¹B.

**Exact algebra (CHECKED to 1e−86 at every L/form):**
- S₊ = (A − B)A⁻¹(A + B) = Q_o (Q_e + Q_o)⁻¹ Q_e.
- A^{-1/2} S₊ A^{-1/2} = I − M², with M = A^{-1/2} B A^{-1/2}.
- So when A ≻ 0: S₊ ⪰ 0 ⇔ ‖M‖ ≤ 1 ⇔ Q_e ⪰ 0 and Q_o ⪰ 0, and n₋(S₊) = n₋(Q_e) + n₋(Q_o).

| L | ζ: min S₊ / q_min,e / q_min,o (rel) | ζ 1−‖M‖ | DH min S₊ | F_0.1: A n₋; S₊ n₋ = n₋(Q_e)+n₋(Q_o) | H₋₂: S₊ n₋ |
|---|---|---|---|---|---|
| 0.6 | 3.0023e−10 / 3.0070e−10 / 1.1e−7 | 5.5e−9 | 2.3e−2 | A PD; 2 = 2 | 0 (PD) |
| 0.8 | 3.370e−18 / 3.372e−18 / 3.3e−15 | 6.2e−17 | 1.3e−3 | A PD; 5 = 5 | 3 = 3 |
| 1.0 | 6.008e−23 / 6.014e−23 / 3.7e−20 | 1.1e−21 | 7.1e−6 | **A n₋ 1**; S₊ n₋ 8 (Q: 9) | 4 = 4 |
| 1.2 | 7.851e−27 / 7.857e−27 / 6.3e−24 | 1.6e−25 | 1.6e−9 | **A n₋ 2**; S₊ n₋ 10 (Q: 12) | 5 = 5 |

(N 24 q_min at L ≥ 1 are unconverged upper bounds; see the S91 caveat.)

- **Verdict (CHECKED, finite algebra):** Turbo's S₊ ≥ 0 *is* window Weil positivity of both parities, rewritten in the half-line/reflection basis. A is the form on a window of half the length (translation invariance of f ↦ W(f * f̃)), so A ≻ 0 is the positivity statement at L/2. "RP reduced to one Schur inequality" therefore restates the problem rather than reducing it. It is the reflection-basis version of the parity split (CC 2106.01715 Lemma 2.5) plus Haynsworth. It is also not our Σ_S (that split is S vs S^⊥, dim ≤ 10). Both are Schur complements of the same Q, both ⇔ Q ⪰ 0, and for ζ min S₊ ≈ min Σ_S ≈ q_min.
- "F_θ fails only in the sewing term" holds automatically whenever A ≻ 0, since a failure must then show up in S₊. In our normalisation F_0.1's half-line block is PD at L 0.6–0.8 but **not PD at L ≥ 1.0** (n₋ 1, 2). So Turbo's claim fixes their L (still waiting on the convention). It does not contradict S94 (F_θ fails on Q_SS), because that is a different split.

## (2) H94a: Q_SS in the raw prolate-vector basis c_n = Π_N E(φ_n)
- ζ: the Rayleigh quotient Q(Eφ_n)/‖Eφ_n‖² is (0.45–1.7) × (1 − |χ_n|) (in λ_max units) wherever 1 − |χ_n| is resolvable in float64 (≳ 1e−14): all n at L 0.6, n ≥ 4 at 0.8, top 2–3 at L 1.0–1.2 (1.2 e top: 61, the outlier). Below 1e−14 our χ is float noise, so the log–log correlations (0.47–0.96) and slopes (≈ 0.03) in summary.txt are **not meaningful**.
- Diagonality: the largest off-diagonal ratio |Q_nm|/√(Q_nn Q_mm) for ζ is 0.13–0.23 at L 0.6–0.8 but 0.96–0.99 at L 1.0–1.2. So **"≈ diagonal" fails at L ≥ 1**: the prolate vectors become nearly collinear in the Q-metric.
- Fakes on the same vectors: DH Rayleigh/defect is 10 to 1e13 (PD), H₋₂ is 0.8 to 1e12 (PD), and F_0.1 has negative diagonal entries at every L. This is the same separation as S92 (sign for F_θ, size only for DH and H₋ₜ).
- **DEAD END** as a separator, and the mechanism (smallness from the Slepian–Pollak angle between P_λ and P̂_λ) is CC 2106.01715 §3 prior art. 2310.18423 (semilocal prolate/Sonin) has no diagonal statement. An honest test would need high-precision χ_n (mp prolates), which was not done.
