# nogo_attack (2026-10-17): DH vs ζ at support 2.7–3.0, and L²-normalised Ritz upper bounds

RH OPEN. GRH OPEN. Note: GPP-bridge `work/grok/2026-10-17_nogo_and_equivalences_attack.md`.

- `s47_dh.py` (S47, unchanged): Arb W_DH and W_ζ on the Legendre·(1−u²)⁴ basis; mpmath Ritz λ_min in the **plain coefficient norm**. A negative value triggers an Arb witness. Outputs `dh_M128.jsonl` (L = 1.35, 1.5) and `dh_M160.jsonl` (L = 1.35; the 1.5 and 1.75 runs were stopped).
  - **NUMERICAL:** DH is PD at every point tested. DH λ ≈ 3e−21 (M = 128) and 1e−22 (M = 160), almost L-independent. ζ is 20–40 orders smaller.
- `ritz_l2.py M L dps prec`: an L²-normalised generalised Ritz problem (Gram = L·G_jk(0)). The vector is rounded to rationals and the Rayleigh quotient is evaluated in Arb, which gives a **rigorous upper bound** on inf Q/‖f‖² per parity sector.
  - Validation at L = 0.8, M = 56: the even value is 2.96e−17. Zhu's certified bracket is [8.9e−18, 2.27e−17], so our value sits above both ends, as an upper bound from a smaller basis must (CONFIRMED normalisation).
  - Output: `ritz_M160_L135.jsonl`.
- `s14_exact.py`, `s44_scan.py`, `s46_planted.py` are copied unchanged from s47.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
