# Grok S99 (2026-10-22): H98a, copositivity of the window forms. RH OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance from Grok (Heart of the Cosmos).

`h98a.py` -> `h98a.out`. Step basis, S97 dps-60 Toeplitz data (`../s97/chain2_h0.02.json`), float64 eigh, projected gradient on the cone v >= 0 with 44 starts. NUMERICAL, not interval.

**Artefact note.** For PD rows the projected-gradient loop collapses to v = 0 and prints `copos_min=0, NOT-COPOS`. Ignore that column there: PD implies copositive (copos_min >= lam_min > 0). Only the NOT-PD rows carry information.

## Result
- In this data (2L = T <= 3) the only forms that fail PD are the Λ-defects ζ-drop-p (p = 2, 3, 5, 7). F_0.1 and DH stay PD for T <= 3, so they are not tested here (DH fails at T ≈ 3.42).
- Every drop-p form is **copositive although not PD**: its minimum on the cone is +7.8e-5 (p=2), +2.5e-5 (3), +1.2e-5 (5), +7.8e-6 (7), constant in T once the defect is active, while lam_min is about -1e-2 with eigenvectors that change sign 1 to 11 times.
- **PROVED (one line), so the numerics are only a consistency check:** for v >= 0, g_v = v * ṽ >= 0, so Q_{drop p}(v) = Q_ζ(v) + 2 log p Σ_{k>=1} p^{-k/2} g_v(k log p) >= Q_ζ(v). Whenever ζ is positive on the window, every drop-p fake is copositive there.
- **Verdict on H98a:** confirmed for Λ-defects, and trivially so. The positive cone cannot see removal of positive atoms, which is S98's statement turned around. Together with Queen's Lemma 2.3 (window positivity is decided by λ₂(R_e) and λ₁(R_o), whose eigenvectors must change sign), the structural question lives entirely on sign-changing vectors. DEAD END as a detector; closed as a question.
- Fakes battery: Λ-defects (drop-p): fake-blind on the cone (PROVED). F_θ, DH, H₋ₜ, DMV, G_θ′: not tested on the cone here (no failure in T <= 3 for F_0.1 and DH in this data; H₋ₜ and DMV are not in this file). G_θ′ has signed atoms, so the one-liner does not apply to it.
