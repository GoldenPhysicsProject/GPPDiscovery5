# s83 — Grok S83 (2026-10-22). RH OPEN.
- `chain_s83.sh` (detached, nohup, first launch 22:37:42 EDT with K 1024 × 4 workers, stopped by the RSS guard at 22:41, see run/*rssguard.log; relaunched ~22:43 EDT with K 256 × 3 workers, chain PID 15311): H80c step 3. Rebuild the even support-2.84 dump
  (S71 reference ../s71/run/T_L7150_m12_O1000_p0.json, sha b609d6f5…) with **batched Arb cells** (now `--chunk 256`, 3 workers: 120 cell jobs instead of 456 of 64; unchanged s71/mori_cert.py, 4 workers, RSS guard kills at MemAvailable < 1.2 GB), then
  `cmp_dumps.py` (entrywise |mid_new − mid_ref| ≤ radT_new + radT_ref, row-sum radius bounds), `lam_from_dump.py`,
  s36 cholres at 0.999 (must pass) / 1.005 (control, must fail), s37 exact LDLᵀ at the same two shifts. Bitwise sha equality is NOT expected
  (different summation order changes midpoints); the within-radius comparison + re-certification is the test. Logs in run/.
- `h82osb_quadruple.py` → `h82osb_results.json` (float64, < 1 min): H82-OS-b probe on a zero-side toy fake (ζ kernel + one off-line quadruple
  1/2 ± δ ± iγ, K(r) = 4cosh(δr)cos(γr)); NOT the real Diamond–Montgomery–Vorhauer prime system.
