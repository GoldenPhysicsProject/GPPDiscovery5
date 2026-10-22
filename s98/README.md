# Grok S98: H97a test (sign of the prime atoms on ζ's soft mode). RH OPEN.

`python h97a.py` writes `h97a.out` from the S97 Toeplitz data. It runs in about 15 s.

- **NUMERICAL:** ζ's soft mode on [0, T] is nodeless at every tested T (44/44 cases, h = 0.02 and 0.01, T from log p⁺ to 3).
- **Automatic consequence (exact):** if every v_i > 0, then its autocorrelation g_v is > 0 at every lag |z| < T. So every atom with a_n > 0 contributes −2a_n g_v(log n)/√n < 0.
- **Result:** H97a is FALSE. On ζ's soft mode each Λ(p) atom lowers the form; it does not add a positive local jump. Q_p(v)/λ_min runs from −0.02 just past log p down to about −250 at T = log p + 0.4. The arch and pole part carries all the positivity there.
- **DH:** on its own soft mode (also nodeless), the negative coefficients (n = 3, 4, 9, 12, 13, 14, 19) help and the positive ones hurt. So Λ ≥ 0 is the *worst* sign for the soft mode, not the forcing sign.
- **Drop-p failure (S97) is a different direction.** The drop-p form fails on a vector where Q_p > 0. This is automatic, since ζ is PD there.
- **Fakes:** this is not a separator. It says nothing about F_θ or H₋ₜ, which have the same Λ. For DH it separates by size only, again.
