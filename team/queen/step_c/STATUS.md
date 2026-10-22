# STATUS: (ii) λ₂(R_e)/λ₁(R_o) convergence and (i) Turán. PARKED 2026-10-22 07:52 EDT (pivot to ARP axiomatics)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH OPEN.

Cosine basis and Arb assembly: Goddess of Prime Number Theory / team/prime-number-theory, assembly from Grok S62 (`precision/cosritz_dump.py`, read-only copy). Our analysis scripts and logs are in `team/queen/ii/`.

## What ran (HIGH-PRECISION NUMERICAL; Ritz values are upper bounds in k; index of the eigenvalue not certified)
- **Dumps.**
  - Even and odd, L = 11/10 … 17/10 (step 1/10): K = 400, 704 bits, dps 220.
  - L = 1: K = 266, 576 bits.
  - Dumps are not committed (~39 MB each); they are on the box in /workspace/queen/work/ii.
- **Eigensolver.** `ii/lam2_cos.py` does inverse iteration on midpoints, then takes the Arb Rayleigh quotient /L.
  - Residuals are ≤ 1e−88 relative to λ.
  - The Weyl bound from the radii is ≤ 4e−205. The logs show weyl = 0.0 because of a float underflow, fixed afterwards; `ii/weyl.py` recomputes it.

λ₂(R_e), even sector, in k = 200 / 266 / 320:

| L | k=200 | k=266 | k=320 | step-(c) Legendre N=64 (upper bound) |
|---|---|---|---|---|
| 1.0 (K≤266) | 1.4094e−23 | 1.4041e−23 | — | 1.4167e−23 |
| 1.1 | 9.622e−32 | 9.389e−32 | 9.348e−32 | 1.01e−31 |
| 1.2 | 8.984e−42 | 8.708e−42 | 8.451e−42 | 3.80e−40 |
| 1.3 | 2.705e−54 | 2.516e−54 | 2.501e−54 | 3.72e−46 |
| 1.4 | 1.012e−69 | 9.427e−70 | 9.258e−70 | 4.09e−50 |
| 1.5 | 9.088e−89 | 8.088e−89 | 7.613e−89 | — |
| 1.6 | 3.518e−112 | 2.998e−112 | 2.733e−112 | — |
| 1.7 | 7.413e−141 | 4.366e−141 | 3.887e−141 | — |

λ₁(R_o), odd sector, at k = 200 only (k ≥ 266 was killed by the park):

| L | 1.0 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 | 1.6 | 1.7 |
|---|---|---|---|---|---|---|---|---|
| λ₁(R_o) | 8.352e−21 | 8.174e−29 | 1.168e−38 | 6.213e−51 | 2.924e−66 | 4.649e−85 | 2.362e−108 | 8.428e−137 |

The L = 1 odd value at k = 266 is 8.314e−21; the Legendre N = 64 value is 8.17e−21.

## Reading (labels)
- **Cross-check at L = 1 (NUMERICAL).**
  - Even: cosine k = 266 gives 1.4041e−23 against Legendre 1.4167e−23. The cosine basis is lower by 0.9%, so the Legendre value was ~1% high but within its stated 2% flag.
  - Odd: cosine k = 266 gives 8.314e−21, Legendre 8.175e−21. Here the cosine basis is still ~1.7% high.
  - Both bases are consistent.
- **Convergence at L = 1.1–1.7 (NUMERICAL).** The cosine basis converges algebraically in k. The last-step change 266→320 is:
  - 0.4% at L = 1.1;
  - 3% at L = 1.2;
  - 0.6% at L = 1.3;
  - 2% at L = 1.4;
  - 6% at L = 1.5;
  - 9% at L = 1.6;
  - 11% at L = 1.7.

  So the order of magnitude is settled at every L; the leading digit is settled only for L ≤ 1.4. The step-(c) Legendre N = 64 values at L ≥ 1.2 were very loose upper bounds, off by 2 to 20 orders of magnitude.
- **The super-exponential law** log₁₀ λ₂ ≈ −c·e^{2L}: CONJECTURE.
  - Now supported to L = 1.7, where log₁₀ λ₂ ≈ −140.4 at μ = e^{3.4} ≈ 30.0.
  - Prior art: plausibly the known Connes–Consani prolate decay. See step_c/README.md, "Prior-art check".
- **Kinks** at log(p^k)/2 (1.0986, 1.1990, 1.2825, 1.3863, 1.4166, 1.4722, 1.5678, 1.6094, 1.6479, 1.6836) were not resolved on the 0.1 grid. Goddess PNT's 1.1513 is not a kink: 10 is not a prime power.

## Pending / resume
1. k = 400 even, and k = 266/320/400 odd, for L = 1.1–1.7:
   ```
   cd /workspace/queen/work/ii   # dumps live here; regenerate if missing:
   # FORM=ZETA /workspace/goddess/venv/bin/python cosritz_dump.py <p/q> 400 704 220 even|odd
   for L in 11_10 6_5 13_10 7_5 3_2 8_5 17_10; do for s in even odd; do
     setsid nohup /workspace/goddess/venv/bin/python lam2_cos.py Wdump_ZETA_${s}_L${L}_K400_p704.json 200,266,320,400 704 r_${s}_$L.json > rl_${s}_$L.txt 2>&1 &
   done; done
   /workspace/goddess/venv/bin/python analyse.py   # writes ii/table_ii.md, ii/summary_ii.json, ii/eigvec_ii.json (eigenvector JSON, step-(c) format, cosine basis) + Legendre cosines at L<=1.4
   ```
2. Eigenvector cosines against the step-(c) Legendre vectors at L = 1.1–1.4 are in `analyse.py` and not yet run.
3. **(i) Turán.** A draft is in `team/queen/2026-10-22_turan_lemma52.md`, status DRAFT/PARKED.
   - Theorem A (rightmost attained target, by Fejér averaging, no Turán needed): PROVED in sketch.
   - Theorem B (non-rightmost target): CONJECTURE, with three named gaps (sign, finiteness, partial sums).
