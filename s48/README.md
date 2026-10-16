# S48 (Grok): Turbo's signed log-mean gate for d_ψ — Landau reduction + exact drift
- `s48_sieve.py XMAX`: segmented sieve for ψ, ΣΛ/√n, ΣΛ√n, giving G, E₁, d_ψ, log-mean and linear mean at 5·10^k and 10^k checkpoints (run to 1e10: `sieve_1e10.log`, `sieve_10000000000.json`). NUMERICAL (float64).
- `s48_explicit.py`: c_G = 2 − 2ζ'/ζ(½) = 2 − c* (mpmath, 40 digits), the zero identity Σ1/(ρ(ρ−½)) = −(2+γ−log4π) with 100k zeros + Trudgian tail, the RH envelope 0.0462007, and RH-explicit predictions vs sieve (`explicit.log`, `s48_explicit.json`).
- `s48_dh.py N`: DH analogue G_f (q=5; φ via κ only, INPUT-CHOICE), real-axis scan of f, and envelope crossings (`dh_2e7.log`, `s48_dh_2e7.json`).
Note: GPP-bridge work/grok/2026-10-16_s48_dpsi_logmean_landau.md. φ absent (except DH κ input). RH OPEN.
