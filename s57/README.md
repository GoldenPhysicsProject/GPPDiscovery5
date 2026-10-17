# Grok S57 (2026-10-17): support 2.76 odd collected; knockout witness certified; knockout margins converge; odd-floor anatomy

See GPP-bridge work/grok/2026-10-17_s57_support276_knockout_witness.md for results and labels. RH OPEN. GRH OPEN.

- `../s56/run138/cert_L138_ODD_S57_FINAL.json`: odd parity of L = 69/50 (support 2.76), T = 1160, K = 1305, collected in S57.
- `s57_witness.py` (Arb ball arithmetic, S14/S53 assembly) and `s57_witness_check.py` (independent: mpmath only; no flint, no moment code):
  certificate that the Weil form of zeta with the prime-5 term removed is NEGATIVE on an explicit odd f supported in [-L, L].
  Outputs in `wit/` (coefficients are exact rationals in the json). `wit/check.jsonl` is the independent evaluation.
- `s57_archcheck.py`: checks the real-space archimedean formula against the digamma/Fourier form (Gaussians, 40 digits) and the full
  explicit-formula normalization (pole - primes + arch = 1.3e-28 at dps 30 for exp(-x^2/2), where the zero side is 2e-43).
- `s57_ko.py`, `ko/`: knockout margins delta_p with the unconstrained edge basis (EDGEK = 0, Legendre), M = 40/56/80/104, and Ritz
  lambda_min(Q_zeta) at L = (1/2) log p.
- `s57_spec.py`, `spec/`: spectral anatomy of the odd/even Weil minimizer (fhat via spherical Bessel functions), plus a zero-side
  check (partial sums over the first 600 zeros, `spec/zeros600.json`).
