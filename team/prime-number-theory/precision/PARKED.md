# PARKED numerics (2026-10-22 07:55 EDT, on the lead's pivot). RH is OPEN.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Wdump_*.json files are not committed (regenerate with cosritz_dump.py).

Done before parking:
- G_θ' fakes rows (3c7f048).
- Leak export: four runs, `onset_leak.md`.

## Tail bound / Mori strip (MAIN item 2): data ready, write-up not done
- `schur_tail.py` → `runs/tail/schur_even_K400.json` (NUMERICAL), ζ even at L = 1.7105, K400 dump p704:
  - μ(K0) = λ_min(Q_HH) = 0.42–0.69 for K0 = 100–350;
  - η = ‖Q_HL‖ ≈ 1.05–2.1;
  - η²/μ ≈ 1.6–10.7. Against λ_K ≈ 1e-153, the naive Schur lift fails by about 153 orders.
  - The soft vector (K = 266) couples to modes 266..399 with ‖r‖ = 2.9e-77; the second-order shift is 1.7e-154.
  - λ(K) in team units is 3.76e-153 (K200), 1.81e-153 (K266) and 1.41e-153 (K400). An extrapolated fit ~K^-6 gives λ∞ ≈ 1.37e-153 (CONJECTURE-grade).
- `mori_strip.py` → `runs/tail/mori_strip_3421_2000.json` (ESTIMATE). Settings: s = 5.64, from λ_max(C_p) even = 5.629 at K400 (`cp_norm.py`); odd 4.02, still rising.
  - Astra bound Ω > 2916.
  - Ω = 6180 → N_min = 13210; Ω = 7140 → N_min = 15212 (≈ 1.25 ΩL). The team rule 1.671ΩL gives 17.7k–20.4k.
- Cost (ESTIMATE, from s78 scaling at 576 bits):
  - ≈ 1.0–1.3e7 CPU-s, i.e. ≈ 15–19 days on 8 dedicated cores. The 3.5-day plan is undersized by ~5×.
  - Memory: ≈ 13–15 GB per Np² ball matrix at 576 bits, with ≥ 6 live. 16 GB is short by ~5–6× without out-of-core storage or Ozaki slicing.
- Resume: write `precision/tail_bound.md` from the three JSONs above. Rerun with
  `python schur_tail.py runs/tail/Wdump_ZETA_even_L3421_2000_K400_p704.json runs/cos/vec_ZETA_L17105_K266.json "60,100,140,200,266,300,350" runs/tail/schur_even_K400.json`
  and `python mori_strip.py 3421/2000 5.64 "2950,3500,4500,6180,7140" 576`.
- Pending README fix: the λ values in README §1/§3/§4 tables are λ_W = L·(Q/‖f‖²). Divide by L for team units. Thresholds move by less than 1 bit.

## Onset-law split: partial
- DH on-line zeros are done to T = 250 (`runs/dhzeros/z_250.json`, 161 zeros).
- Chunks 250–1450 were killed with no output. Resume:
  `cd runs/dhzeros; python ../../dh_zeros.py 250 450 0.08 z_450.json` (and 450–650, …, 1250–1450; run in parallel).
- Then `python ../onset_split.py <vec.json> dhzeros/z_250.json,dhzeros/z_450.json,... 1450 out.json`.
- ζ f₀ quartet grid is done (`runs/split/zf0_*`).
- The P_X split at 1.6043 was not started. Plan: ζ zero side plus the exact atom term from the difference between the edge and edgeref dumps.

## G1 (Conjecture G1 numerical test): not run
- `g1/g1_test.py alpha X seed` is written but untested. Smoke test: `python g1/g1_test.py 0.6 1e6 1`, then X = 2e8 with several seeds.
