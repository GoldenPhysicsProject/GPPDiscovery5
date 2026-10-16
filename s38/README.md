# Grok S38 — Separating test: L=1.3 with D8×2 lift A≤4.4702

RH OPEN. GRH OPEN. φ absent.

## H38a (PRIMARY) — VERDICT: odd stays PD → S37 lift-induced conjecture **WRONG** (for this test)

Re-certified support 2.6 (L=13/10, T=700, K=740, ng=78) using the **D8×2** lift
`../s36/torus_L1.3_D8x2_p235_p71113_cert.json` (A_cert ≤ 4.4702), not the S36 D4+D20 lift (A≤4.6565).

| parity | λ_est | cholres λ_min(Ã) ≥ | ldlt_exact ≥ | rE | certified_lower |
|---|---|---|---|---|---|
| even | 4.264e−62 | 4.260e−62 | 4.260e−62 | 1.53e−72 | **4.260e−62** |
| odd  | 3.717e−58 | 3.713e−58 | 3.713e−58 | 1.53e−72 | **3.713e−58** |

- β(700) = 0.2416 (was 0.0553 with A=4.6565). M_S ≤ 764 (was 8.68e6). Wall ~24 min on 2 workers.
- Negative controls at 1.02·λ_est FAIL (even pivot 261, odd pivot 237).
- **Even improves** vs S36: 4.260e−62 vs 3.634e−62 (~+17%).
- **Odd stays strongly PD** (3.713e−58). So the S37 odd failure at L=1.35 is **not** explained by D8×2 lifts in general; it is support/geometry at L=1.35 (or the specific L=1.35 torus), not a generic D8×2 defect.

## H38b (light) — L=1.4 ladder design only
- 16 = 2⁴ enters the comb when L > log(16)/2 ≈ 1.386294; prime set stays {2,3,5,7,11,13} until L > log(17)/2 ≈ 1.4166.
- Proposed bundles: D8 {2,3,5} + D8 {7,11,13}. No LP started (post H38a first). See `h38b_L14_design.json`.

## Artifacts
- `run/cert_L13_T700_K740_D8x2_FINAL.json` — combined certificate
- `run/atil_sha256.txt` — dump hashes (75 MB dumps not committed)
- `s38_cholres.py`, `s38_ldlt_exact.py`, `s38_onestroke_schurboost.py`

## Dead ends
- None new for H38a. Confirmed dead end remains: treating D8×2 as the cause of odd@L=1.35 without this separating test.
