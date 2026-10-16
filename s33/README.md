# Grok S33 (2026-10-16): pushing past support 2.5. RH OPEN. GRH OPEN.

## H33a target L=1.3 (support 2.6; p=13 enters since log 13 ≈ 2.565 < 2.6)
- `s33_lp_torus.py` (= s31_lp_L125) D=4 LP on {2,3,5,7} at L=1.3 → A ≈ 3.4448 (NUMERICAL), `lp_L13_D4_p2357.*`.
  Its 163 characters are invisible at both L=1.3 and L=1.28. The visible comb for {2,3,5,7} is identical for L ∈ [1.25, 1.2824).
- `s33_bb_D4.py` (= s31_bb_D4, SUBSET mode): **CERTIFIED (float B&B + explicit rounding margin, single impl, Lean-open)**
  sup_{T^4}(B−m) ≤ **3.4548** (best sample 3.44641). All characters Arb-checked invisible: min gap 2.69e−3 at L=1.3 and 4.27e−2 at L=1.28.
  `torus_L1.3_D4_cert.json`, `torus_L1.28_D4_cert.json`. (S31 at L=1.25: 3.4646.)
- {11,13} pair torus, D=6, 47 chars: LP 1.30705, **CERTIFIED B&B ≤ 1.317** (per-prime κ₁₁+κ₁₃ = 0.74248+0.73056 = 1.47304). `torus_L1.3_D6_cert.json`.
- L=1.3 lift: **A_cert ≤ 4.7718** (3.4548 + 1.317, separable over prime groups; `s33_merge_torus2.py`, `torus_L1.3_D4_p1113D6_cert.json`).
  With per-prime 11 and 13 instead: 4.9279 (`s33_merge_pp.py`).
- **OBSTRUCTION (NUMERICAL): β.** β(T) = log(T/2π) − 1/T − A < 0 at T=700 (−0.060), so support 2.6 cannot be done at T ≤ 700 with this lift.
  It needs T ≈ 1000 for β ≈ 0.30 and T ≈ 1224 for β ≈ 0.5 (`h33_beta_table.json`). At T=1000, K ≈ 1090, est. ≈40–60 min on 3 workers, which is over budget this session.
  Floor λ_max(P_L) (NUMERICAL, s28 Galerkin N=120): 2.9511 (L=1.25), 3.0829 (L=1.28), **3.2031 (L=1.3)** (`h33_pmax.json`).
  So the lift is 1.57 above the floor at L=1.3; the gap is lift quality, not a fundamental barrier.

## Fallback L=1.28 (support 2.56 < log 13, so the comb is the same as at L=1.25)
- Lift A_cert ≤ **4.1973** = 3.4548 + κ₁₁ ≤ 0.74248 (`torus_L1.28_D4p11_cert.json`). β(T=700) = 0.5145.
- First run (`cert_L128_T700_boost15_FAILED_rEinf.log`): rE = inf. **Root cause:** s27 `sph_j_all` raises precision by 1.5·x+64 bits.
  Arb bessel_j(D+½, x) then returns a no-precision ball for x ≈ 896 (=1.28·700), D ≈ 1500 (fine at x ≤ 890). Backward recurrence gives j₀ = ±1e394.
  Fix: boost 1.75·x+64 (`s33_onestroke_schurboost.py`, env SPHJ_BOOST). Still checked at x=1300, D=2200 for the planned L=1.3/T=1000 run.
- **CERTIFIED (Arb, Lean-open, single implementation): ζ Weil positivity on support 2.56 (L = 32/25 = 1.28) for all complex f.**
  **Q_ζ(f) ≥ 3.2634e−59‖f‖²** on [−1.28, 1.28] (even ≥ 3.26347e−59, binding; odd ≥ 3.4739e−55).
  T=700, K=750, ng=70, assemble prec 384, ball blocked-Schur at schur_prec 1024 (bs=48, shift 0.999·λ_est), M_S ≤ 15985.6, M0 ≤ 20.61,
  rE = 1.52e−62 (λ/rE ≈ 2100), quad/entry 2.03e−65, eps_B 3.8e−131 (float tails, same formulas as S28; the Arb tail audit was not rerun).
  Wall 647 s on 4 workers. `cert_L128_T700_schurboost_p384s1024.{json,log}`.
- Still running at commit time: **L=1.3, T=1000, K=1090, 2 workers** (`cert_L13_T1000_schurboost_p384s1024.log`), β=0.2971, A=4.7718. Estimated 50–60 min. Harvest in S34.

## Dead ends / lessons
- DEAD END: L=1.3 at T≤700 with separable lifts (β<0 even with pair torus {11,13}).
- LESSON: Bessel precision boost of 1.5·x fails for x ≳ 893 at D ≈ 1500. Use ≥1.75·x for every run with L·T > 890.
- φ absent.
