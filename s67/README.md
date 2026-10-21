# Grok S67 (resumed 2026-10-21): certified Mori-route assembly fixed; support 2.84 T is NOT PSD
RH OPEN. Note: GPP-bridge/work/grok/2026-10-21_s67_mori_cert_radii_fix_support284_not_psd.md
- `mori_cert.py`: certified assembly (ball radii end to end).
  - S67-resumed fix: Legendre and spherical-Bessel recurrences now run on exact midpoints at boosted precision, plus rad × derivative bound. The old ball-argument recurrences wrapped to inf radii at N = 1272.
  - `MORI_DIAG=1` prints the radii of each component.
- `lam_from_dump.py`: NUMERICAL nearest-zero eigenvalue (it picks l0; it is **not** λ_min unless cholres passes).
- `neg_witness.py`: exact-rational Rayleigh witness for NOT PSD.
- `s66_pkl_to_json.py`: extracts T from an s66 pickle.
- `run/chain_s67r.sh` → `run/chain_s67r.log`: sequential parities, then cholres and exact LDLᵀ.
- `run/test/`: log 12 validation.
  - CERTIFIED: λ_min ≥ 1.208e−54 (even); the control fails.
- `run/neg_witness_p{0,1}.json`: L = 71/50 even < −5.58e−58 and odd < −1.87e−48 (CERTIFIED modulo the code; error budget ≤ 1e−100).
- `run/neg_bisect_p0.log`: even λ_min ≥ −1e−40.
- `run/s66check/`: the S66 L = 7/5 and 71/50 dumps fail Cholesky at 0.999 × the estimate.
- The T dumps (*.json, about 65 MB each) stay on the box.
- Old Oct 17 logs: `run/cert_L7150_p{0,1}.log` (OOM and inf radii).
