# Grok S55 (2026-10-17): buffered lift at L = 1.38 (A ≈ 4.954, NUMERICAL); support 2.701 at T = 800 (buffer rule test); knockout margins δ_p ≈ 0.016/log p

RH OPEN. GRH OPEN. φ absent. Code and data: GPPDiscovery5/s55 (commit @D5@).

## H55a: buffered, richer lift at L = 1.38 (support 2.76). Design NUMERICAL, not certified
- New tool `s55_bundle_lp.py` = S36 bundle LP with a buffer: every lift character needs |k·log p| ≥ 2L + ε, ε = 0.0056 (ε·T = 5 at T = 900). At L = 1.38 the closest allowed characters are log 16 = 2.77259 ({2,3,5}, k = (4,0,0)) and 2.82080 ({7,11,13}). log 15, 3log(5/2) and log(169/11) are all excluded.
- Result after 900 s each (warm start from S54; still drifting down about 0.01 per accepted step):
  - {2,3,5}, D9 (461 characters): sampled sup **3.0595** (trust-region LP 3.0566).
  - {7,11,13}, D8 (362 characters): **1.8941** (LP 1.8747).
  - **A ≈ 4.954** (NUMERICAL, NG = 128 grid plus gradient refinement). That is 0.40 better than the S54 re-solve (5.353), but **not ≤ 4.75**.
- LP lower bound (`s55_lp_lower.py`, an unrestricted LP on 25.8k torus points, so it is a valid lower bound for the character set up to LP tolerance). {2,3,5} D8 buffered: **A₁ ≥ 2.883**. Together with the {7,11,13} LP value ≈ 1.87, a D8 separable buffered lift is ≳ 4.75. So D8 cannot reach the target; D ≥ 9 is needed even in the best case. The D10 and D12 runs were killed (memory and CPU, no result).
- T needed with A ≈ 4.954 (+0.004 B&B margin): β = 0.12 (enough for the even sector at L = 1.35) needs **T ≈ 1009**; β = 0.25 (the odd sector needed it at L = 1.35) needs **T ≈ 1149**. That is down from T ≈ 1750 in S54.
- I did not run B&B certification of this lift or assembly at L = 1.38 (time). **Support 2.76 is NOT certified.**

## H55a′: frontier test, L = 2701/2000 = 1.3505 (support 2.701), T = 800, certified S37 lift (A ≤ 4.5928)
- Gap to log 15: ε = 7.05e−3, εT = 5.64. This tests the S54 buffer rule ε·T ≳ 5 at a new (L, T). S54 L = 1.354 had εT = 0.04, and **both** parities failed (odd at pivot 35, found after S54 write-up; λ_est was a misleading +8.2e−66).
- **EVEN: CERTIFIED (Arb, Lean-open, two PD implementations).**
  - cholres p640: λ_min(Ã) ≥ 1.2000e−69. ldlt_exact S640: ≥ 1.2003e−69, ‖E‖_F ≤ 9.5e−115.
  - rE = 5.32e−87. Negative control at 1.2255e−69: @NEG@.
  - Q_ζ ≥ 1.2000e−69 ‖f‖² on even L²[−1.3505, 1.3505]. That is 2.3× the S37 even value at L = 1.35 (T = 700).
- **ODD:** @ODD@
- **Full complex at support 2.701:** @FULL@
- The buffer rule is confirmed in the even sector (HEURISTIC rule, now with 3 data points: εT = 5.6 pass, 5.64 pass, 0.04 fail).
- Using the S37 lift, the rule caps the frontier at 2L ≤ log 15 − 5/T. That is 2.7018 at T = 800, so 2.701 is essentially the end for this lift.

## Fix: lam_estimate
- `s55_lam.py`: lam_estimate_fixed(Ã) returns the smallest-modulus μ (inverse iteration) **and** the Sylvester inertia, i.e. the number of negative pivots of a blocked midpoint LDLᵀ. μ is reported as λ_min only if there are 0 negative pivots, because a PD matrix's smallest-modulus eigenvalue is its smallest eigenvalue. Otherwise it reports INDEFINITE.
- Unit test passed: a 60×60 matrix with eigenvalues {1e−40, −1e−30, 1, …} was flagged INDEFINITE, while the old routine would report +1e−40. The run on the S54 L = 1.354 even dump: @LAMFIX@.
- Old certificates are unaffected (PD was always decided by cholres + ldlt_exact + a failing negative control).

## H55b: knockout margins δ_p = L_p* − ½log p (NUMERICAL, Ritz M = 56; fits CONJECTURE)
- Refined grid (0.005 in L) for p = 5, 11, 13, plus new p = 17, 19. `s55_knockout.py`, `h55b/`, `s55_ko_fit.py` → `h55b/ko_fit.json`.
- A negative Ritz value is a real witness, so each δ interval's right end is an upper bound.

| p | ½log p | odd δ_p ∈ | even δ_p ∈ |
|---|---|---|---|
| 2 | 0.3466 | (0.023, 0.028] | (0.063, 0.068] |
| 3 | 0.5493 | (0.011, 0.016] | (0.026, 0.031] |
| 5 | 0.8047 | (0.008, 0.010] | (0.013, 0.015] |
| 7 | 0.9730 | (0.007, 0.012] | (0.012, 0.017] |
| 11 | 1.1989 | (0.006, 0.011] | (0.006, 0.011] |
| 13 | 1.2825 | (0.0025, 0.0075] | (0.0075, 0.0125] |
| 17 | 1.4166 | (0.0034, 0.0084] | (0.0084, 0.0134] |
| 19 | 1.4722 | (0.0028, 0.0078] | (0.0078, 0.0128] |

- **Fits (log-log on midpoints):**
  - **odd: δ_p ≈ 0.0165·(log p)^−1.02** (rms log-residual 0.14). Equivalently δ_p·log p ∈ [0.013, 0.021] for all 8 primes, so **CONJECTURE: δ_p^odd ≍ c/log p with c ≈ 0.016**.
  - even: δ_p ≈ 0.034·(log p)^−1.29 (rms 0.18).
  - Versus p: odd exponent −0.62 (rms 0.18), worse.
  - Versus w_p = log p/√p: meaningless, because w_p is not monotone (max at p = 7). REFUTED as an organizing variable.
- The odd sector fails first for every p except p = 11 (tie at this resolution). The even margin is larger by about one grid step (0.005) for p ≥ 13.
- HEURISTIC reading: e^{2L_p*} ≈ p(1 + 2δ_p) ≈ p + 0.033·p/log p. The knocked-out form fails once the support passes p by a relative amount ~1/log² p. The δ ≳ 0.002 resolution is a limitation, and M = 56 may overstate δ (M = 72 lowered p = 3 odd). Whether δ_p → 0 is OPEN.
- φ absent (no constant here involves the golden ratio).

## Dead ends
- DEAD END: a D8 separable buffered {2,3,5}+{7,11,13} lift at L = 1.38 reaching A ≤ 4.75 (LP lower bound ≳ 4.75).
- DEAD END (this session): A ≤ 4.75 at L = 1.38 with D9/D8 (A ≈ 4.954 after 900 s, slow drift).
- DEAD END: the S37 lift beyond 2L = log 15 − 5/T (≈ 2.7018 at T = 800).
- REFUTED: w_p = log p/√p as the scaling variable for knockout margins.

## Next (S56)
(a) B&B-certify the buffered A ≈ 4.954 lift (`s54/lift/lift.sh` with the s55 bundles). Then run L = 1.38 even at T ≈ 1010 (β ≈ 0.12) and odd at T ≈ 1150 (β ≈ 0.25), with detached runs or checkpointing; one parity per session if needed.
(b) Or push D10 on {2,3,5} with a memory-light LP (column generation) to bring A toward 4.8.
(c) H55b: M = 72 at p = 13, 17, 19 to test whether δ_p·log p is stable.
