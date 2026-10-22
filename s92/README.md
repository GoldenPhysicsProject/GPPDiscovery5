# Grok S92 (2026-10-22): backward heat flow H₋ₜ on the Connes–Consani near-radical span. RH OPEN.

Status: **HIGH-PRECISION NUMERICAL** (mpmath dps 90 forms, float64 prolates and angles). Not interval, nothing proved.

## What was run
`ht_espan.py` (needs `../s91/forms_cache.pkl` symlinked as `forms_cache.pkl` in this folder, see LARGE_LOCAL_FILES/s91). For each L ∈ {0.6, 0.8, 1.0, 1.2} and each parity it takes the S91 span S = span E(φ_n) (CC 2106.01715 §3, ν from |χ| ≥ 0.9) and the backward-heat-flow model Q_HT(t) = Q_ζ + Δ_t from s89 (`band_split.delta_ht`: the first MZ = 50 zeros of ζ are swapped for those of H_t, t = −1, −2, −4). Output: `ht_espan.json`, `ht_espan.log`.

## Results (relative to λ_max of Q_ζ)
| L | par | ζ on S | t=−1 full / on S | t=−2 full (n₋) / on S | t=−4 full (n₋) / on S (n₋) | sin(ground of H₋₂, S) |
|---|---|---|---|---|---|---|
| 0.6 | e | 4.1e−10 | 1.8e−10 / 2.2e−7 | 7.8e−11 (0) / 3.0e−7 | −0.20 (4) / 3.9e−7 (0) | 0.007 |
| 0.6 | o | 1.4e−7 | 5.7e−8 / 2.1e−5 | 2.5e−8 (0) / 3.5e−5 | −0.17 (4) / 4.3e−5 (0) | 0.026 |
| 0.8 | e | 5.4e−18 | 1.1e−18 / 5.5e−11 | −1.9e−8 (2) / 1.0e−10 | −0.10 (4) / 1.3e−10 (0) | 0.18 |
| 0.8 | o | 4.4e−15 | 1.1e−15 / 8.5e−9 | −5.4e−7 (1) / 1.5e−8 | −0.087 (4) / 1.3e−8 (0) | 0.24 |
| 1.0 | e | 7.2e−23 | 5.6e−23 / 2.0e−15 | −6.3e−9 (2) / 2.3e−15 | −0.75 (5) / 3.0e−16 (0) | 0.29 |
| 1.0 | o | 4.2e−20 | 4.9e−20 / 8.0e−13 | −1.9e−4 (2) / 9.9e−13 | −0.74 (5) / 2.2e−14 (0) | 0.41 |
| 1.2 | e | 1.1e−26 | 2.6e−27 / 4.3e−21 | −1.6e−8 (3) / 4.2e−21 | −1.49 (6) / **−0.047 (3)** | 0.43 |
| 1.2 | o | 6.5e−24 | 3.3e−24 / 2.8e−18 | −7.0e−11 (2) / 2.7e−18 | −1.42 (5) / **−0.27 (2)** | 0.17 |

## Reading
- H₋₂ and H₋₄ fail on the full window form from L = 0.8 (H₋₄ already at 0.6), but **S does not witness the failure** except H₋₄ at L = 1.2: on S the H₋ₜ forms stay positive, and their failing ground vectors sit at sines 0.2–1 from S.
- On S, H₋ₜ is 10³–10⁵ times *less* null than ζ (e.g. 4e−21 against 1e−26 at L = 1.2 even), so S stops being a near-radical once zeros move. Like DH (S91), H₋ₜ is separated from ζ by S only in size, not by a sign.
- Tally with S91: S witnesses F_θ at every L (cosh weight breaks the Poisson pairing E(f̂)(x) = E(f)(1/x)), but not DH or H₋ₜ, both of which keep the functional equation. **So the CC near-radical span is a functional-equation detector, not an RH detector** (NUMERICAL). Expected from its construction (it is built from Poisson symmetry alone), but now checked against all three standard fakes.
- Caveats: Δ_t swaps only the first 50 zeros, so Q_HT is a truncated model of H₋ₜ; N = 24 basis is unconverged for q_min at L ≥ 1 (S91 caveat); t = −1 does not fail at these L.
