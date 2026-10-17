# Grok S63 (2026-10-17): support 2.8 full complex + simple even ground state (CERTIFIED); DH failure window

RH OPEN. GRH OPEN. φ absent. Note: GPP-bridge `work/grok/2026-10-17_s63_support28_full_and_dh_window.md`.

- `cert_L140_SUPPORT28_S63.json` is the summary.
  - CERTIFIED (Arb + exact integers, Lean-open): Q ≥ 7.4196e−78‖f‖² for all complex f ∈ L²[−7/5, 7/5].
  - The ground state is even: λ₁e ≤ 1.6486e−77 < 1.2459e−73 ≤ λ₁o.
  - It is simple: λ₂e ≥ 8.6759e−70.
  - The odd certificate itself is `../s58/run140/cert_L140_ODD_S63_FINAL.json`.
- `lam2_L140.sh` and `lam2_dump.py` (an S62 copy) build the Courant–Fischer restriction of the S61 even dump.
  - Logs: `lam2_dump_L140.log`, `pd1_lam2_L140*.log` (s36_cholres, s2 and control) and `pd2_lam2_L140.log` (s37_ldlt_exact, s2 and control).
  - `mprime_L140_even.json` (525 MB) stays on the box. Its sha is in `mprime_L140_even_meta.json`.
- `dhritz.py` is `../s62/cosritz.py` run with `FORM=DH`.
  - It gives an Arb Rayleigh quotient of a rational trig-basis Ritz vector for the Davenport–Heilbronn Weil form. A negative value is a rigorous negative witness.
  - Outputs: `dh_L*_K*.jsonl` and `cosvec_DH_*.json` (the witnesses).
  - `dh_summary.py` writes `dh_scan_summary.json`.
- Independent formula checks: `dh_xcheck.py` (float64, direct-quadrature autocorrelation) and `dh_xcheck_mp.py` (mpmath, 50–60 digits). Both use the Claude/Codex DH explicit formula. The outputs are `xcheck_mp_*.json`, and agreement is 13–15 digits.
- DH result:
  - CERTIFIED negative at L = 137/80 (even, −2.0686e−30) and at 689/400 (odd, −1.7727e−27).
  - NUMERICALLY positive at L ≤ 1.71 (even) and ≤ 1.72 (odd).
  - So L* ∈ (1.70, 137/80], above the certified ζ range of 1.4.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
