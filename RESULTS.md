# Results

Tables written up from `results.jsonl`, newest first. Every row says whether it is floating-point or certified.

| date | experiment | result | label | float/certified |
|---|---|---|---|---|
| 2026-10-14 | `weil_small_support/weil_ritz.py` (Ritz, Legendre×(1−u²)⁴, M=14–30) | normalization checked vs 300 ζ zeros (agrees to 1e-12). Even sector is the weak one. ζ: λ_min>0 for h∈[0.2,1.0] (float floor ~1e-14 beyond). F_{1/4}: >0 at h≤0.44 (M=30: +1.5e-5 at 0.44; +3.1e-3 at log2/2), <0 at h≥0.46 (M=18), so h* ∈ (0.44,0.46]. F_{1/20}: >0 at h≤0.44 (M=30), <0 at h≥0.53 (M=18). DH: >0 for h≤1.9 (M=14), unresolved beyond (float floor). All fakes positive at h=log2/2, where no primes enter | NUMERICAL (Ritz = upper bounds on inf; positivity not proved) | float |
| 2026-10-14 | `weil_small_support/weil_ritz.py` (+ `zerocheck.py`: form = 1.59655748142e-4 vs 300-zero sum 1.59655748142e-4) | Ritz min of Weil form on supp ⊂ [−h,h]. All four positive at h ≤ log2/2 = 0.3466 (prime-free window; M = 14/24/30). F_{1/4} first negative (even sector): h = 0.465 (M=14), ≤ 0.46 (M=18), still positive at 0.44 (M=30). F_{1/20} first negative: h = 0.53 (M=18), 0.56 (M=14). DH positive up to h = 1.9 (M=14); beyond that the float noise floor is 1e-14, unresolved. ζ positive to noise floor. Reproduces CC's F_{1/4} −4.1e-5 at h=0.5 | NUMERICAL (Ritz upper bounds; threshold drifts down with M) | float |
| 2026-10-14 | `symmetrized_euler/sym_euler.py N 60` | analytic Λ_N(s)+Λ_N(1−s): off-line zeros (Re>½, t≤60) for ζ: 6/17/32/97/327 at N=2/3/5/13/50; F_{1/4} and DH similar, so no separation | NUMERICAL | float |
| 2026-10-14 | `li_certified/li_arb.py 2001 6000` | λ_1..λ_2001(F_{1/4}) all > 0 (balls exclude 0); L=10600 run for sign change near 9250 in progress | CERTIFIED | certified (Arb) |
| 2026-10-14 | `dh_phi/dh_phi.py` | DH κ = √(φ+2)−φ = 5^{1/4}√φ−φ; ε(χ)=(φ+i)/√(φ+2); φ from choice q=5 | PROVED LEMMA (algebra) + Arb consistency | certified check |
| 2026-10-14 | `ftheta_window` h=19/25 (even) | F_{1/4} −0.021969145467, ζ +0.002030094923, matches expected | NUMERICAL | float |
| 2026-10-14 | `ftheta_window` h=17/25 (odd) | F_{1/4} −0.008863529402, ζ +0.000809109424, matches expected | NUMERICAL | float |
| 2026-10-14 | `li_coefficients/li_ftheta_onset.py 0.25 600` | λ_n(F_{1/4}) is first negative near n≈9250 (250-step grid); −5.5e5 at n=10500 | NUMERICAL (ξ part asymptotic) | float |
| 2026-10-14 | `li_coefficients/li_coeffs.py 100` | λ_1..λ_100 > 0 for ξ, F_{1/4}, F_{1/20}, DH; λ_1(ξ)=0.023096 | NUMERICAL | float |
| 2026-10-14 | `dh_zero/dh_zero.py` | DH zero 0.808517182456637 + 85.699348485377592i, \|f\|≈8e-32, FE residual 2.5e-32 | NUMERICAL | float (dps≈30) |
| 2026-10-14 | `pole_parity/pole_parity.py` | Weil pole term: even +1.71887, odd −0.00528, mixed 1.69775 ≠ 2G(½)²=2.12099 (pole = 2G(½)G(−½)) | NUMERICAL | float |
| 2026-10-14 | `weil_ab_split` M=18 h=0.40..0.52 | F_{1/4} fails h=0.46; on min, dPole>|dA| but dPrimes tips dQ<0; primes/C crosses 1 for F_{1/4} only | NUMERICAL | float |
