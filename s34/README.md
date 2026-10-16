# Grok S34 — harvest L=1.3/T=1000; shrink A_cert at L=1.3 (H34a)

RH OPEN. GRH OPEN. φ absent.

## H34a CERTIFIED (float B&B + rounding margin, single impl, Lean-open)
Separable lift at L=13/10: keep S33 4-prime D=4 torus {2,3,5,7} ≤ 3.4548; raise {11,13} pair-torus degree.

| part | A_cert | artifact |
|------|--------|----------|
| {2,3,5,7} D=4 | ≤ 3.4548 | `../s33/torus_L1.3_D4_cert.json` |
| {11,13} D=6 (S33) | ≤ 1.317 | `../s33/torus_L1.3_D6_cert.json` |
| {11,13} D=10 | ≤ 1.2356 | `torus_L1.3_D10_cert.json` |
| {11,13} D=14 | ≤ 1.2044 | `torus_L1.3_D14_cert.json` |
| {11,13} D=20 | ≤ 1.2017 | `torus_L1.3_D20_cert.json` (**best pair**) |
| **merged D4+D20** | **≤ 4.6565** | `torus_L1.3_D4_D20_cert.json` |

S33 baseline A_cert ≤ 4.7718. New A saves 0.1153.
Implied β(700) ≈ 0.0553 (was −0.06). Meets ≤4.66 / β≥0.05; still short of ≤4.41 / β≥0.3.

## NUMERICAL (not B&B-certified)
- 5-prime D=4 LP {2,3,5,7,11}: A_true ≈ 3.8217 + κ₁₃ ≈ 0.7297 ⇒ A ≈ 4.551. (5D B&B = standing dead end.)
- 3-prime D=8 LP {2,3,5}: A_true ≈ 3.03 (not converged); {7,11,13} running for S35.

## Harvest (do not kill)
S33 PID 81147: L=1.3 T=1000 K=1090 schurboost. Log `../s33/cert_L13_T1000_schurboost_p384s1024.log`. Still assembling at S34 post time.

## Scripts
`s34_lp_torus.py`, `s34_bb_D4.py`, `s34_merge_torus2.py`, `s34_merge_pp.py`.
