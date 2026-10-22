# Grok S93 (2026-10-22): the Weil window form on S^⊥, the complement of the Connes–Consani near-radical span. RH OPEN.

Status: **HIGH-PRECISION NUMERICAL** (mpmath dps 90 forms, N = 24 unless stated; float64 prolates define S). Positivity claims are
**CHECKED** by mp Cholesky at dps 90 (also at half the margin), and failures by LDLᵀ negative pivots. This is not interval arithmetic and nothing is proved.

## What was run
- `sperp.py` → `sperp.log`, `sperp.json`, about 2 min on 4 cores. It needs `forms_cache.pkl`, a symlink to `../s91/forms_cache.pkl` that is not committed.
  - S = span E(φ_n), the S91/S92 CC span (2106.01715 §3, ν from |χ| ≥ 0.9). Π_N(S) is its projection into our Legendre coefficient space, and U is an exact dps-90 orthonormal basis of its complement S^⊥.
  - Forms: ζ, DH, F_0.1, H₋₂ and H₋₄ (zero swap, MZ = 50). For ζ we also ran these prime perturbations:
    - remove or sign-flip each prime with all its powers;
    - scale Λ by 1 ± 1%.
    The rebuilt prime matrix matches the cache to 1e−65.
  - For each form we computed the min eigenvalue and n₋ on the full space, on S and on S^⊥, plus the Haynsworth split of Q|S^⊥ along the arch A|S^⊥. We also compressed to P_⊥(ζ's eigenspace with eig < 1e−6 λ_max).
- `followup.py` → `followup.log`, `followup.json`, about 3 min on 7 cores. It computes three things:
  - the pole-free R = A − P on S^⊥ and ρ_⊥ = λ_max(A_⊥^{−1/2} P_⊥ A_⊥^{−1/2});
  - the complementary exact split Q ⪰ 0 ⇔ Q_⊥ ≻ 0 and Σ_S = Q_SS − Q_S⊥ Q_⊥⁻¹ Q_⊥S ⪰ 0;
  - ζ on S^⊥ at N = 32 (L 0.6, 1.0, 1.2).

## Results (min eigenvalue / λ_max(Q_ζ); n₋ in brackets)
| L | par | dim S / S^⊥ | ζ full | ζ on S^⊥ | F_0.1 full / S / S^⊥ | DH full / S^⊥ | H₋₂ full / S / S^⊥ / Σ_S | H₋₄ full / S / S^⊥ |
|---|---|---|---|---|---|---|---|---|
| 0.6 | e | 2/22 | 3.4e−10 | 0.046 | −2.0e−5 (1) / −1.3e−5 / **+0.046** | 0.029 / 0.47 | 7.8e−11 / 3.0e−7 / 0.017 / 7.8e−11 | −0.20 (4) / 3.9e−7 / **−0.20 (2)** |
| 0.6 | o | 2/22 | 1.1e−7 | 0.140 | −4.9e−6 (1) / −4.7e−6 / +0.140 | 0.23 / 0.46 | 2.5e−8 / 3.5e−5 / 0.082 / 2.5e−8 | −0.17 (4) / 4.3e−5 / **−0.17 (2)** |
| 0.8 | e | 4/20 | 3.9e−18 | 0.142 | −6.2e−4 (3) / −5.3e−4 / +0.141 | 1.5e−3 / 0.49 | −1.9e−8 (2) / 1.0e−10 / 0.058 / **−1.9e−8 (2)** | −0.10 (4) / 1.3e−10 / **−0.10 (2)** |
| 0.8 | o | 4/20 | 3.3e−15 | 0.159 | −6.7e−4 (2) / −6.6e−4 / +0.158 | 0.088 / 0.46 | −5.4e−7 (1) / 1.5e−8 / 0.055 / **−5.7e−7 (1)** | −0.087 (4) / 1.3e−8 / **−0.084 (2)** |
| 1.0* | e | 7/17 | 6.8e−23 | 0.202 | −2.6e−3 (5) / −2.6e−3 / +0.201 | 8.1e−6 / 0.46 | −6.3e−9 (2) / 2.3e−15 / 0.015 / **−6.8e−9 (2)** | −0.75 (5) / 3.0e−16 / **−0.72 (3)** |
| 1.0* | o | 6/18 | 3.7e−20 | 0.135 | −2.2e−3 (4) / −2.2e−3 / +0.134 | 2.5e−3 / 0.38 | −1.9e−4 (2) / 9.9e−13 / 0.023 / **−2.2e−4 (2)** | −0.74 (5) / 2.2e−14 / **−0.74 (3)** |
| 1.2* | e | 10/14 | 1.0e−26 | 0.385 | −4.8e−3 (6) / −4.8e−3 / +0.384 | 2.1e−9 / 0.43 | −1.6e−8 (3) / 4.2e−21 / 0.214 / **−2.0e−8 (3)** | −1.49 (6) / −0.047 / **−0.97 (2)** |
| 1.2* | o | 10/14 | 6.3e−24 | 0.361 | −3.5e−3 (6) / −3.5e−3 / +0.361 | 1.7e−6 / 0.35 | −7.0e−11 (2) / 2.7e−18 / 0.182 / **−7.2e−11 (2)** | −1.42 (5) / −0.27 / **−1.18 (2)** |

\* L ≥ 1: N = 24 is unconverged for the q_min (S91 caveat). On S^⊥ the margin is O(0.1) and stays O(0.1) at N = 32: ζ on S^⊥ is 0.041 / 0.127 at L 0.6, 0.159 / 0.104 at L 1.0 and 0.218 / 0.202 at L 1.2 (e / o). So the S^⊥ margin is stable in sign but not converged in value at L ≥ 1.

Additional facts:
- **The in-band Schur complement restricted to S^⊥ is vacuous.** The arch A compressed to S^⊥ has *no* negative direction at any L (min +0.08 to +0.21 rel), while A has 1–3 per parity on the full space. So S absorbs the whole Γ-band ghost span, and on S^⊥ the Haynsworth split reduces to "out-of-band block PD".
- **The Λ perturbations in (c) never fail on S^⊥.** Every removal or sign flip of p = 2, 3, 5, 7 and both 1% rescalings leave ζ on S^⊥ positive (0.02–0.40), while the full form and Q|S go negative (−4e−3 down to −0.39; 1–7 directions). This includes flipping p = 2 at L 0.6, which gives 0.024 (e) and 0.098 (o). Removing or flipping p = 11 at L = 1.2 changes nothing (window edge, as in S90).
- **Pole-free on S^⊥ (ζ):** R = A − P is PD on S^⊥ (Cholesky) for all odd L and for even L ≥ 0.8, with ρ_⊥ = 0.20–0.65. The one exception is even L = 0.6, where ρ_⊥ = 1.19 and R is −0.023. Compare S89: on the Γ-positive band of A, ρ_out was 1.02–1.85 for ζ even, a dead end. ‖P_⊥‖ / min A_⊥ = 1.06–2.39 > 1, so this is not a sign-blind norm bound.
- **Exact complementary split:** every form except H₋₄ has Q_⊥ ≻ 0 with margin, so Q ⪰ 0 ⇔ Σ_S ⪰ 0 (dim ≤ 10). For ζ, min Σ_S = q_min (to 4 digits). F_θ, H₋₂ and DH are decided entirely by Σ_S. H₋₂ is positive on S *and* on S^⊥ but negative on Σ_S, so it fails **only through the S–S^⊥ coupling**.

## Answers
- (a) ζ: S^⊥ is PD with margin 0.05–0.38 (CHECKED, Cholesky dps 90, also at half the margin), with n₋ = 0. The fakes on S^⊥:
  - F_0.1: PD (n₋ 0).
  - DH: PD (n₋ 0, 0.35–0.49).
  - H₋₂: PD (n₋ 0).
  - H₋₄: negative (n₋ 2–3, LDLᵀ confirms) at every L from 0.6.
  The control that must fail does fail: H₋₄ on S^⊥, and F_θ / H₋₂ on the full form and on Σ_S (Cholesky fails).
- (b) By Cauchy interlacing, S^⊥ can never fail at a smaller L than the full form. F_θ never fails on S^⊥; it fails on S at the same L as the full form. H₋₂ never fails on S^⊥ or on S; it fails only via the coupling. H₋₄ fails on S^⊥ at the same L as the full form (0.6), and before S (1.2). DH does not fail anywhere at L ≤ 1.2.
- (c) ζ's tightness is entirely in S. Positivity on S^⊥ survives every single-prime removal or flip and both 1% rescalings, so it does **not** depend on Λ ≥ 0 at this level. All Λ-rigidity (S90) lives in S / Σ_S.

## Verdict
**DEAD END as an RH detector.** S^⊥ rules out only H₋₄, which is a gross fake already broken at O(1). It does not rule out F_θ, DH or H₋₂. The structural by-product is the decomposition: Weil(L) ⇔ [S^⊥: arch dominates primes, pole-free for L ≥ 0.8, robust, Λ-sign-insensitive] + [Σ_S ⪰ 0 on the ≤ 10-dim CC span, where all tightness, all Λ-rigidity and every failure of F_θ / H₋₂ live]. This is standard Haynsworth on top of CC 2106.01715 §3, so it is not new. It also does not reduce dimension asymptotically, since dim S ≈ 2μ = 2e^{2L} grows.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
