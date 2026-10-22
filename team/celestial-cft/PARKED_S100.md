# PARKED at S100 (2026-10-22, 07:50 EDT): all numerics stopped on lead/Daniel direction. RH is OPEN.
No process is running. Every in-flight run below was killed at 07:50 EDT. Nothing has been started since.

## Finished
| item | label | where |
|---|---|---|
| λ_min(Q_L) Galerkin fit, 2x²+2xy+3y²: slope 3.18 on L ∈ [1.5, 3] (half-width) vs 2Θ − 1 = 3.159, ρ₀ = 2.0794141 + 7.6288676i | NUMERICAL | `queen_probe/LAMBDA_MIN.md` (commit 54308f7) |
| Probe-vector η fit | not ours: Queen team/queen 7bdce47 (η̂ = 3.20); deduped | — |
| Dominant-term abscissa, explicit lattice-point tail: no zeros with Re s ≥ 2.4102 for 2x²+2xy+3y² (lead replayed; merged); all 12 forms tabulated | CERTIFIED | `epstein/DOMINANT_ABSCISSA.md` (644472a) |
| Rigorous point evaluator for any form: `eps.py` repaired (tail was multiplied by \|pref\| twice; τ_Q formed at raised precision; `Kbest` works around Arb K_ν radius blow-ups) | validated against `zq_L.py` (D = −20, −84, t ≤ 200) and the direct lattice sum (D = −23) | `epstein/eps.py` |
| **D = −23 and D = −56: 43 of 43 located zeros (t ≤ 100) isolated, each in a disk of radius 1e−8, by Rouché** (lattice-point M bound + eps.py). Per form: (1,1,6) 12/12, max Re 0.9533; (2,1,3) 5/5, max Re 1.2454; (1,0,14) 7/7; (2,0,7) 15/15, max Re 1.1835; (3,2,5) 4/4. Negative control (center off by 1e−6) correctly fails. | CERTIFIED (each zero). Completeness of the zero list for these forms is NUMERICAL only. | `epstein/certify_rouche.py`, `certified_zeros_rouche_D23_D56.json`, `certify_rouche.log` |
| Prior art: Potter–Titchmarsh 1935, Davenport–Heilbronn 1936 I/II, Bombieri–Hejhal 1995, Bombieri–Mueller 2008 (pp. 359–385, per Lee arXiv:1204.6297). BCFMP full-text grep: no off-line statement. | citations checked | `epstein/EPSTEIN_DH_CFT.md` §2 |

## Partial
1. **Certified count to t = 201** (genus forms, region [0.52, 3] × [1, 201], boxes of height 10; argument principle with Arb).
   - 138 of 140 boxes are done, in `epstein/certified_counts_T201.jsonl`. Re s ≥ 3 is zero-free (CERTIFIED).

   | form | count, t ∈ [1, 201] (CERTIFIED) | of which t ∈ [1, 101] | numerical list t ≤ 100 |
   |---|---|---|---|
   | x²+y² | 0 | 0 | 0 |
   | x²+5y² | 35 | 13 | 13 |
   | 2x²+2xy+3y² | 30 | 12 | 12 |
   | x²+21y² | 29 | 9 | 9 |
   | 3x²+7y² | 48 | 22 | 21 (*) |
   | 2x²+2xy+11y² | 51 | 19 | 19 |
   | 5x²+4xy+5y² | 17, **boxes [181,191] and [191,201] missing** | 7 | 7 |

   - (*) The certified box [91, 101] reaches t = 101, but the numerical list stops at t = 100. The extra zero is presumably in 100 < t ≤ 101. Unchecked.
   - Resume (skips done boxes): `cd team/celestial-cft/epstein && python cert_count_resumable.py 201 3`
2. **x²+5y² zero with Re s > 1** (Θ > 1 anchor).
   - Found **ρ = 1.0110493340435507 + 1193.9464045288462i**: a Newton-refined zero of the exact 1 + ζL₋₂₀/(L₋₄L₅) (Arb Hurwitz L-values). NUMERICAL.
   - The certificate was started and interrupted. Not certified yet.
   - Earlier float winding found no zero with Re s > 1.01 for t ≤ 1000. That run's output was not committed, and it does not cover 1 < Re s ≤ 1.01.
   - So "first" is NOT established. The hunt (`hunt_re_gt1.py 1.03 0 300000`) was killed partway, with only this one hit logged.
   - Bracket if it certifies: 1.0110 ≤ Θ(x²+5y²) ≤ 1.452. The upper bound is CERTIFIED (DOMINANT_ABSCISSA.md).
   - Bohr heuristic: σ(Q) ≈ 1.1–1.15 (the phase budget 2 atan 2^{−σ} + 4 Σ_{p≡3,7 (20)} atan p^{−σ} = π). HEURISTIC.
   - Resume: `python certify_gt1.py 1.0110493340435507 1193.9464045288462 1e-6`. Then, for "first": `python search_re_gt1.py 1.0 1200 0.004 3` (NUMERICAL), or a certified count on [1, 3] × [1, 1200] with `certify_L.count`.
3. **λ_min for x²+5y².** Converged values (half-width L; Θ_cert = 0.9377 from t ≤ 100):

   | L | λ_min | N |
   |---|---|---|
   | 1 | +1.0532 | 64 |
   | 2 | −1.6379 | 128 |
   | 3 | −8.6439 | 128 |
   | 4 | −29.219 | 160, change 0.2% vs N = 128 |

   - Slopes 1.66 (L 2→3) and 1.22 (L 3→4). NUMERICAL.
   - L = 5 has not converged (−67 at N = 64, still moving).
   - The relevant Θ is > 1 (see item 2), so this fit should be redone against 2·1.011 − 1 = 1.022 once item 2 is certified.
   - Resume: `cd queen_probe/hiN && python galerkin.py 1,0,5 5 96,128,160`

## Not started
- Task 3 (ζ_K, Z_Q) Weil-form setup spec, for lead review.

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
