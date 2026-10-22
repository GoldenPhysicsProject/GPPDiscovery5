# Grok S82-OS (2026-10-21): the reflection-positivity (Osterwalder–Schrader) structure of the ζ window Weil form

RH OPEN. Everything here is `kind: float` (numpy float64) except the S75 Rayleigh quotients it quotes, which are Arb-CERTIFIED. φ absent.
Note: GPP-bridge `work/grok/2026-10-21_reflection_positivity.md`. Python: numpy only (run with any venv that has numpy).
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

- `os_hankel.py FORM L N`: builds, in an orthonormal shifted-Legendre basis on [0, L], the same-side form T(h) = ∫∫h h W(x−y) and the
  OS/Hankel form B(h) = ⟨Θh, h⟩ = ∫∫h h W(x+y), with Θ: u ↦ −u (x ↦ 1/x). W is the S14-normalised Weil distribution (ZETA or DH, as in `s75/dhritz.py`).
  Even sector = T + B, odd = T − B; Weil on [−L, L] ⇔ ‖T^{−1/2}BT^{−1/2}‖ ≤ 1.
  Cross-check at L = ½: λ_min(T±B) = 9.368e−7 / 1.945e−4 (N = 24) against the Arb cos/sin Ritz values 1.0374e−6 / 2.2046e−4 (dhritz K = 30).
- `scan.py` → `scan_results.json`: ζ and DH, L = 0.2–1.0. B is indefinite for both (it is unconditionally not ≥ 0; the archimedean Hankel part
  −Σ_k(∫h e^{−(2k+½)x})² ≤ 0 is proved, and max eig ≤ 4e−15 numerically). Also the contraction defect 1 − ‖C_L‖ and the identity k1 + k2 = e^{s/2}/sinh s.
- `perron_pole.py` → `perron_pole_results.json`: pole-free form A = Q − pole. For ζ, n₋(A⁺) = 1 with a nodeless ground state (overlap with cosh(u/2) ≥ 0.98),
  n₋(A⁻) = 0, and the second even level collapses (5.6e−8 at L = 0.7).
- `witness_shape.py` → `witness_shape_results.json`: the certified DH negative witness at L = 1.7105 is nodeless (0 sign changes on [0, L]).
  ζ's near-null vector at 137/80 is strongly sign-changing.
Runtime < 30 s total, single process.
