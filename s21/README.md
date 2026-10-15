# S21 — complex-safety of E_comb (closes the S20c GAP) + ζ first-window minimizer locks onto ζ zeros

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN. GRH OPEN.**

- `s21_zeta_minimizer_zeros.py M h...` (h = decimal, a/b, or `hstar` = ½log3): Ritz minimizer of Q_ζ (S17 builder,
  real g, per parity), zero-side check 2Σ_{γ<200}|ĝ(γ)|²/‖g‖², real roots of ĝ_min on (0,60) vs ζ ordinates, and the
  share of the zero sum carried by γ₁..γ₆. NUMERICAL (mpmath dps 60 on Arb midpoints; ĝ by 600-pt Gauss–Legendre).
- Outputs: `log_M32.txt`, `log_M48.txt`, `out_minimizer_M{32,48}.json`.
- The complex-safety lemma is a proof (see the bridge note); no script is needed for it.
