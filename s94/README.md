# Grok S94 (2026-10-22): H93a — anatomy of the S–S^⊥ coupling C = Q_S⊥ Q_⊥⁻¹ Q_⊥S on the Connes–Consani span. RH OPEN.

Status: **HIGH-PRECISION NUMERICAL** (mpmath dps 90 forms, N = 24 unless stated; float64 prolates define S, as in S91–S93).
Positivity of Σ_S = Q_SS − C is **CHECKED** by mp Cholesky at dps 90 (also at half the margin); failures by LDLᵀ negative pivots. Not interval arithmetic; nothing is proved.

## What was run
- `coupling.py` → `coupling.log`, `coupling.json` (~1 m 50 s, 4 processes). Needs `forms_cache.pkl` = symlink to `../s91/forms_cache.pkl` (not committed). Reuses s93 `sperp.py` (orthobasis, complement, Cholesky, LDLᵀ), s92 span, s89 zero swap (MZ 50), s90 prime atoms (P rebuild error ≤ 6e−65).
  - Forms: ζ, DH, F_0.1, H₋₂, H₋₄; ζ with p = 2, 3, 5, 7 removed / sign-flipped (all powers) and Λ × (1 ± 1%). L 0.6, 0.8, 1.0, 1.2, both parities.
  - Off-diagonal block X = B_Sᵀ Q U split **linearly** by source: X = X_a (arch) + X_p (−P_Λ) + X_k (pole) [+ X_d zero swap for H₋ₜ]. With G = Q_⊥⁻¹ (the form's own, not split), C = Σ_ij X_i G X_jᵀ; reported pieces aa, pp, kk, dd (PSD automatically) and the symmetrised crosses ap, ak, pk, ad, … (piece sum reproduces C to 1e−86).
  - Statistics: eigenvalues of C; μ_max = largest generalized eigenvalue of (C, Q_SS) when Q_SS ≻ 0 (Σ_S ⪰ 0 ⇔ μ_max ≤ 1); every piece along v = min eigvec of Σ_S and along v_ζ (ζ's), relative to vᵀQ_SS v; for each form ΔQ_SS and ΔC along v_ζ relative to ζ.
- `dh_extend.py` → `dh_extend.log`, `dh_extend.json`: ζ vs DH at L 1.3, 1.4 (N 24; DH forms built fresh) and ζ at N 32 (L 0.6, 1.0, 1.2).

## Results
| L | par | dim S | ζ μ_max (N24 / N32) | ζ vCv/vQv at v | H₋₂ μ_max; dd/vQv; Σ_S | DH μ_max | F_0.1: Q_SS n₋ / Σ_S n₋ |
|---|---|---|---|---|---|---|---|
| 0.6 | e | 2 | 0.213 / 0.214 | 0.158 | 0.99999; 1.003; PD | 0.421 | 1 / 1 |
| 0.6 | o | 2 | 0.296 / 0.296 | 0.190 | 0.99984; 0.995; PD | 0.023 | 1 / 1 |
| 0.8 | e | 4 | 0.467 | 0.314 | **1.0020**; 1.000; n₋ 2 | 0.429 | 2 / 3 |
| 0.8 | o | 4 | 0.413 | 0.262 | **1.0048**; 0.998; n₋ 1 | 0.139 | 2 / 2 |
| 1.0* | e | 7 | 0.249 / 0.429 | 0.058 | **1.0002**; 0.998; n₋ 2 | 0.502 | 4 / 5 |
| 1.0* | o | 6 | 0.236 / 0.342 | 0.118 | **1.216**; 0.959; n₋ 2 | 0.546 | 4 / 4 |
| 1.2* | e | 10 | 0.357 / **0.9992** | 0.112 | **1.0003**; 0.996; n₋ 3 | 0.979 | 5 / 6 |
| 1.2* | o | 10 | 0.258 / 0.338 | 0.032 | **1.0008**; 0.993; n₋ 2 | 0.894 | 5 / 6 |
| 1.3 | e / o | 13 / 12 | 0.190 / 0.348 (N24) | 0.067 / 0.022 | – | 0.99977 / 0.99605 | – |
| 1.4 | e / o | 16 / 15 | 0.214 / 0.271 (N24) | 0.173 / 0.008 | – | 0.9999996 / 0.999992 | – |

\* N 24 unconverged for q_min at L ≥ 1 (S91 caveat). H₋₄: Q_⊥ not PD (−0.08 to −1.18), split not applicable (S93).

- ζ: Σ_S PD at every L/parity, N 24 and N 32, Cholesky also at half margin (CHECKED). Failing controls: H₋₂ (L ≥ 0.8), F_0.1 and all Λ-perturbations (6 per parity at L ≤ 0.8, 10 at L ≥ 1.0) fail with LDLᵀ n₋ matching the eigen count.
- **ζ's tightness is in Q_SS, not in the coupling.** Along ζ's q_min the coupling is only 0.03–0.31 of vᵀQ_SS v (N 24), and μ_max(ζ) = 0.19–0.47 at N 24 (so Σ_S ⪰ 0.5 Q_SS). Exception: N 32, L 1.2 even gives μ_max 0.9992 with Q_SS min 1.3e−29 — at that level the float64 definition of S (angle ~1e−15) is the floor, so we treat it as a precision artefact (HYPOTHESIS, not verified).
- **H₋₂ fails purely by the zero-swap piece taking back its own S-gain.** Along v the dd piece X_d G X_dᵀ is 0.96–1.00 of vᵀQ_SS v, and along v_ζ ΔQ_SS ≈ ΔC to 3–4 digits (e.g. L 1.0 e: +7.804e−7 vs +7.802e−7). D adds positivity on S that is realisable on S^⊥ and the Schur complement removes it, leaving ζ's tiny q_min minus a margin → μ_max slightly > 1.
- Pieces of ζ's coupling: aa, pp, kk ⪰ 0 automatically (G ≻ 0), crosses indefinite on S: ap ∈ [−0.073, 0.26] (0.6 e), pk ∈ [−0.38, 0.05] (0.6 e), ak ∈ [−0.98, 0.02]; odd ap is almost NSD ([−0.03, −9e−5] at 0.6 o) but has positive eigenvalues at 0.8–1.2 o (up to 4.8e−3). Along v_ζ individual pieces are 1e4–1e24 × vᵀQ_SS v and cancel: no single piece controls the sign.
- Λ-perturbations along v_ζ: ΔC > 0 in all 64 cases, and ±1% scaling gives equal ΔC (second order). This is **automatic**: C ⪰ 0, C_ζ(v_ζ) ≈ 0 because ζ's near-null vector lies in S (S91), and v ↦ vᵀXQ_⊥⁻¹Xᵀv is the jointly convex matrix-fractional function, whose first variation ∝ w* = −Q_⊥⁻¹Q_⊥S v ≈ 0. Removing/flipping p = 3, 5, 7 fails mainly or partly through ΔC (ΔQ_SS comparable in odd parity, 10–10⁴× smaller in even); p = 2 odd and the ±1% scalings fail through ΔQ_SS.

## Answers
1. (H93a, Q1) H₋₂ does fail because its coupling exceeds Q_SS in the q_min direction (μ_max 1.0002–1.22 vs ζ 0.19–0.47), but the excess is the **zero-swap piece itself** (dd ≈ 1), not a Λ-sign effect. ζ's coupling is **not** controlled by a sign-definite Λ ≥ 0 term: the pp piece is PSD for any real Λ (it is quadratic, invariant under Λ → −Λ), and the Λ-linear crosses ap, pk are indefinite. → **DEAD END** (H93a).
2. (Q2) Rules out: μ_max > 1 rules out H₋₂ only, and tautologically (⇔ Σ_S fails when Q_SS ≻ 0). F_θ fails on Q_SS itself (n₋ 1–5) with ΔC ~ 1e−9 — coupling irrelevant; not ruled out by any coupling statistic beyond what Q|S already gives (S92). DH: μ_max → 1 at L 1.2–1.4 (0.98 → 0.9999996) while Σ_S stays PD; this says DH's near-null vector leaves ζ's E-span (E uses Σ_n f(nx), i.e. ζ's Dirichlet series), a "which Dirichlet series" detector, not an RH detector, and not robust (ζ at N 32 L 1.2 e also reaches 0.9992). → **DEAD END** as separator.

## Prior art
Haynsworth (1968) inertia additivity / Schur complement positivity; joint convexity of the matrix-fractional function (Boyd–Vandenberghe §3.1.7); Connes–Consani 2106.01715 §3 for S and the near-radical. Nothing here is new beyond the numbers.

## Next (H94a, HYPOTHESIS)
Since all tightness is in Q_SS = Q(Eφ_n, Eφ_m) on the prolate image, compute Q_SS in the prolate basis itself and test whether it is ≈ diagonal with entries governed by the prolate concentration (1 − |χ_n|) / λ_n of the time–band limiting operator — i.e. in-band positivity = contractivity of Fourier restricted to [−λ, λ] (a unitarity statement), with integers entering only through E and Λ ≥ 0 only via the explicit formula identity. Fakes to test: F_θ (breaks the functional equation, should break the χ_n structure), DH (different Dirichlet series → different E). Prior art to check first: CC 2106.01715 §3–4 and 2310.18423 (prolate wave operator, Sonin space).

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
