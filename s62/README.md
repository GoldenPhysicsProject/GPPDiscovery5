# Grok S62 (2026-10-17): simple even ground state of the ζ Weil form at L = 27/20 (support 2.7), CERTIFIED

RH OPEN. GRH OPEN. φ absent. Note: GPP-bridge `work/grok/2026-10-17_s62_evenness_simplicity_L135.md`.
Summary certificate: `cert_L135_GROUNDSTATE_S62.json`.

- `cosritz.py L K prec dps [even|odd|both]` computes rigorous L²-normalised Ritz upper bounds on the geometric side only, in Arb.
  - Basis: even cos((2m+1)πu/2L), odd sin(mπu/L). Same span as Zhu's sine modes. Gram = L·I.
  - Pole and primes are in closed form. The archimedean term uses acb.integral of analytic (sinc) rewrites of the S14 kernels.
  - The eigenvector is rounded to rationals and the Rayleigh quotient is evaluated in Arb.
  - `FORM=DH` gives the Davenport–Heilbronn form (S47 difference trick). It is NUMERICAL use only.
  - Outputs: `cos_*.jsonl`, `sin_*.jsonl`, `dh_*.jsonl`, and `cosvec_*.json` (the rational vectors).
- `selftest.py` checks the closed-form autocorrelations and pole moments against quadrature (agreement to 4e−31).
- `lam2_dump.py` builds the Courant–Fischer restriction of the S37 even dump. The PD checks then run with `../s36/s36_cholres.py` and `../s37/s37_ldlt_exact.py`.
  - Logs: `pd1_lam2*.log`, `pd2_lam2.log`, `mprime_*_meta.json`.
  - The 142 MB `mprime_L135_even.json` stays on the box. Its sha256 is in the meta file.
- `ritz_w.py`: the S14 Legendre·(1−x²)^w basis (w = 0 vs 4) at small M. It shows that w = 0 is worse, and that both are far from converged (DEAD END for upper bounds).
- `dh_l2_check.py`: L²-normalised DH Ritz on the independent S47 assembly, used to cross-check `FORM=DH`.
- `s14_exact.py`, `s44_scan.py`, `s46_planted.py`, `s47_dh.py` are unchanged copies.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
