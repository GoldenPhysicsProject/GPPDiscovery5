# Grok S37: A_cert(L=27/20) ≤ 4.5928; ζ Weil positivity on support 2.7 holds for the EVEN sector only (odd sector FAILS); second independent PD implementation (exact-integer LDLᵀ) CONFIRMS S36

RH OPEN. GRH OPEN. Finite window only. φ absent.

## H37a: lift at L = 27/20 (support 2.7). CERTIFIED: two implementations per part
The prime comb is unchanged, since 2L = 2.7 < log 16: primes {2,3,5,7,11,13}, prime powers {4,8,9}.
- **{2,3,5} D8**: bundle LP (`../s36/L135`, NG=160) reached 2.89961.
  - Float B&B `s34_bb_D4`: **sup ≤ 2.9008** (best sample 2.899575, 8.24e6 boxes).
  - **CONFIRMED** by the independent libm-free `s35_interval_bb` (ok, 8.24e6 boxes, min invisible gap 8.05e−3).
  - Files: `bb235/`.
- **{7,11,13} D8**: bundle LP reached 1.687625, converged; I stopped it after 16 iterations.
  - Float B&B: **sup ≤ 1.692** (best sample 1.687895, 4.77e7 boxes, MAXBOXES=4e7).
  - **CONFIRMED** by `s35_interval_bb` (ok, 4.77e7 boxes, min gap 3.20e−2).
  - Files: `bb71113b/`. The 12M-box run at 1.69/1.6912 was abandoned; its log is `bb71113/bb_cap12M_abandoned.log`.
- **Merged: A_cert(L=27/20) ≤ 4.5928**, 701 group-local characters.
  - Checks: `s35_check_merge.py` OK; recomputing the merge from the final B&B part gives identical K/a_q/A.
  - File: `torus_L1.35_D8x2_p235_p71113_cert.json`.
  - The run was launched from a PROVISIONAL {7,11,13} file. It has identical K and a_q to the certified part; only the A_cert value 1.692 was pending at launch, and it was then certified.
- β(T) = log(T/2π) − 1/T − A gives **β(700) = 0.118975 > 0**.

## H37a: Weil cert at L = 27/20, T=700, K=780, ng=92, prec 384 (`cert135/`)
- Assembly: `../s36/s36_onestroke_schurboost.py` with ball Schur disabled (`--schur_bs ''`); wall about 16 min on 5 workers.
- Bounds: M_S ≤ 2334.1 (it was 8.68e6 at L=1.3 with the D4+D20 lift, so the D8 lift is far tamer), M0 ≤ 26.003.
- **rE = 3.65e−87** (quadrature 4.68e−90 per entry). eps_B ≤ 6.0e−127; lamD = 0.11897.
- **EVEN: CERTIFIED (Arb, Lean-open, TWO independent PD implementations).** For all even f ∈ L²[−1.35, 1.35], **Q_ζ(f) ≥ 5.2649e−70 ‖f‖²**.
  - λ_est = 5.2706e−70.
  - impl 1 (`s36_cholres`, 640 bits): λ_min ≥ 5.265e−70, ‖E‖ ≤ 3.6e−115.
  - impl 2 (`s37_ldlt_exact`, S=640): λ_min ≥ 5.265e−70, ‖E‖_F ≤ 8.7e−115. Negative control at l0 = 5.4e−70 FAILS (pivot 429).
- **ODD: FAILED, and the failure is real.** The odd Ã is numerically indefinite: λ_est = −3.28e−65.
  - `s37_ldlt_exact` with l0 = 1e−90 gets a negative pivot at index 114.
  - |λ|/rE ≈ 1e22, so neither quadrature nor rounding explains it.
  - At L=1.3 the odd sector was the *easy* one (1.49e−58 vs even 3.6e−62). The sign flip coincides with the lift change from S34 D4+D20 to D8×2.
  - Diagnosis is CONJECTURE: the D8×2 lift's minorant m is lossy in the odd sector, or K/T are insufficient for odd at L=1.35. Not yet separated.
- So **full complex support 2.7 is NOT certified.** The full-complex window stays at support 2.6 (S36).
- Dumps (75 MB each) stay on the box; their sha256 values are in `cert135/atil_sha256.txt`. Combined result: `cert135/cert_L135_T700_K780_FINAL.json`.

## H37b: second, independent PD implementation (H3), `s37_ldlt_exact.py`
- Shares no code with `s36_cholres.py` and uses no Arb or interval arithmetic: only exact integers (Python int, flint fmpz_mat) and exact Fractions.
- Steps:
  1. sha256 check of the dump.
  2. Exact decimal parse, then N = round(2^S x), with |error| ≤ 2^−S−1.
  3. Fixed-point LDLᵀ of sym(N) − c0 I.
  4. **Exact** integer residual E = 2^{2S}(N − c0 I) − L D Lᵀ, with all D > 0 checked exactly.
  5. λ_min ≥ c0/2^S − ‖E‖_F − ‖Ã − N/2^S‖_F (Weyl). The string error is taken as 1 unit in the last printed digit.
- Self-test on Hilbert(40) (λ_min = 2.6635e−59): certifies 2.66e−59, fails at 2.7e−59.
- **CONFIRMED on the S36 dumps.** sha256 matches `s36/atil_sha256.txt`.
  - even: λ_min ≥ 3.634e−62 (‖E‖_F ≤ 8.5e−115), negative control 3.71e−62 FAILS;
  - odd: λ_min ≥ 1.491e−58 (‖E‖_F ≤ 1.1e−114), negative control 1.52e−58 FAILS.
  - Recomputed certified_lower: even ≥ 3.63398e−62 and odd ≥ 1.49099e−58, the same as S36. Summary: `h37b/h37b_summary.json`.
  - Runtime about 1 min per parity.
- **H3 now has 2 independent implementations** for S36 (support 2.6, both parities) and for S37 (even, support 2.7).

## Dead ends
- DEAD END: B&B with the default MAXBOXES=12M for {7,11,13} D8 at target ≤ 1.6912. Use MAXBOXES ≥ 4e7 or margin ≥ 0.004.
- DEAD END (for now): the D8×2 lift at L=1.35/T=700/K=780 for the ODD sector (indefinite at −3.3e−65).
