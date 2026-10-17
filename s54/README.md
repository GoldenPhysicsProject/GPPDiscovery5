# Grok S54: the D8×2 lift family stops at support ≈ 2.7 (the log 15 wall); knockout test: every prime ≤ 13 is necessary

RH OPEN. GRH OPEN. Finite windows only. φ absent.

## H54a: full-complex certificate at L = 1.38 (support 2.76). NOT achieved
### 1. The S37 lift is not admissible at L = 1.38 (PROVED, elementary; Arb-checked)
- A lift character k is allowed only if its frequency |k·log p| > 2L (invisible). The S37 lift's closest characters are:
  - {2,3,5}: k = (0,1,1), frequency **log 15 = 2.7080502**, coefficient **a = 1.9160** (the largest near-edge one); then (3,0,−3), 3·log(5/2) = 2.74887, a = 0.076.
  - {7,11,13}: k = (0,1,−2), log(169/11) = 2.73200, a = 0.645.
- So the S37 lift is valid only for **L < ½·log 15 = 1.3540251**. At L = 1.38 three characters turn visible, including the two biggest. The comb of prime powers is unchanged (n ≤ 13 below ½log16), but the lift is not.
- The B&B sup does not depend on L (B is the fixed comb). So the S37 bound A ≤ 4.5928 stays CERTIFIED for every L < ½log 15.

### 2. Re-solved lift at L = 1.38 (NUMERICAL, not certified)
- `../s36/s36_bundle_lp.py` at L = 1.38, D = 8, NGRID = 160, warm start from the S36/S37 L = 1.35 bundles with the visible characters removed (336 and 362 characters). Logs and JSONs: `lift/`.
- Start: true max 4.884 ({2,3,5}) and 2.333 ({7,11,13}). Best after about 10 minutes: **3.3084 + 2.0450 = 5.353**. That is +0.76 over 4.5928 at L = 1.35. Not converged. I stopped it (budget overrun, CPU needed); the round-2 restart was also stopped.
- β = log(T/2π) − 1/T − A > 0.25 would need **T ≈ 1750** (≈ 1540 for β = 0.12). That is about 2× the T = 800–850 plan, so assembly would be about 4× the cost. I did not run B&B on a lift that cannot be used at T ≈ 850.
- `lift/lift.sh` (float B&B + libm-free interval B&B + merge check) is ready for S55 if a better lift is found.

### 3. Fallback L = 677/500 = 1.354 (just below ½log 15) with the certified S37 lift. EVEN FAILS
- Run: `run1354/` (T = 800, K = 884, ng = 92, prec 384; same lift A ≤ 4.5928, β = 0.2527). The Arb visibility check passes (log 15 − 2L = 5.0e−5 > 0). Bounds: M_S ≤ 2334.0, M0 ≤ 25.87. rE = 5.39e−87.
- **EVEN: FAILED, and the failure is real.** Both PD implementations stop at pivot 56 (cholres at 640 bits with shifts from 0.999·λ_est down to −1.25e−60; ldlt_exact at S = 640). The leading 60×60 block has λ_min = −7.49e−51, and the leading 100×100 block has −1.66e−31. Its eigenvector has its peak at index 17 (even degree ≈ 34). Control: the S37 even dump (L = 1.35, T = 700) has leading-100 λ_min = +5.13e−68. Files: `run1354/block_eigs_even.txt`, `cholres_even_p640*.log`, `ldlt_even.log`.
- **Tooling note (bug in the λ_est printout).** `lam_estimate` is inverse power iteration, so it returns the eigenvalue of *smallest modulus*, not the smallest one. Here it printed +1.2554e−70 for an indefinite matrix. The certificates were never affected: PD is always decided by the two rigorous implementations, and the negative control must fail. But λ_est > 0 must not be read as "PD". Fix for S55: report also the leading-block λ_min or an LDLᵀ inertia count.
- **Diagnosis (HEURISTIC, with the numbers above).** The character (0,1,1) sits at frequency u = log 15, only ε = u − 2L = 5.0e−5 above the support edge. Its ∫_T^∞ a·cos(ut)|f̂|² tail does not cancel: f̂ of an L²[−L,L] function carries e^{±iLt} boundary waves, so the tail has a slow phase e^{iεt} with εT = 0.04 ≪ 1. The loss is then O(a/T), not O(a/(εT²)). At L = 1.35 the gap was ε = 8.05e−3 (εT = 5.6 at T = 700), and the run was fine. **Design rule for future lifts:** keep every character at gap ε with ε·T ≳ 5, so ε ≳ 0.006 at T = 800.
- ODD at L = 1.354 was still assembling when I wrote this (for diagnosis only; full-complex 1.354 is already lost to even). Log: `run1354/run_odd.log`; PD runs automatically after it (`chain2.sh` → `pd_odd.log`).
- **Conclusion.** With the D8×2 character set, support 2.7 (S53) is essentially the end. Going past log 15 needs a lift without (0,1,1), and then A jumps by about 0.76. That needs T ≈ 1750, or a better lift design.

## Dead ends (H54a)
- DEAD END: the S37 D8×2 lift at any L ≥ ½log 15 = 1.35403 (not admissible).
- DEAD END: L just below ½log 15 (1.354) with that lift at T = 800 (edge resonance, even indefinite at −1.7e−31 on the leading block).
- DEAD END (for T ≤ 850): D8×2 lift re-solved at L = 1.38 (A ≈ 5.35 NUMERICAL, needs T ≈ 1750).
## H54b: knockout test of "the odd floor is prime 2" (NUMERICAL)
`s54_knockout.py` reuses the S53 exact Arb assembly (`../s53/s53_odd_floor.py`, Legendre·(1−u²)⁴ basis). For each L it computes the per-sector λ_min of W_ζ with chosen prime-power terms removed. Variants: full, no 2, no {2,4,8}, no 3, no {3,9}, no 5, no 7, no 11, no 13, no primes. Data: `h54b/ko_*.jsonl`, summary `h54b/summary.json`.
A Ritz value is ≥ the true infimum, so a **negative** Ritz value is a real witness that the knocked-out form is indefinite. The thresholds below are upper bounds for the true thresholds.

| removed | ½log p | odd fails at L in | even fails at L in |
|---|---|---|---|
| 2 (also 2,4,8) | 0.3466 | (0.370, 0.375] (M=56 and M=72) | (0.410, 0.415] (M=56 and M=72) |
| 3 (also 3,9) | 0.5493 | ≤ 0.560 at M=72 (M=56: (0.560, 0.565]) | (0.565, 0.575] at M=72 ((0.575, 0.580] at M=56) |
| 5 | 0.8047 | (0.81, 0.82] (M=56 and M=72) | (0.81, 0.82] at M=72 ((0.81, 0.83] at M=56) |
| 7 | 0.9730 | (0.980, 0.985] | (0.985, 0.990] |
| 11 | 1.1989 | (1.20, 1.24] (coarse) | (1.20, 1.24] (coarse) |
| 13 | 1.2825 | (1.28, 1.30] | (1.28, 1.30] |

- **L_2\* (odd) ≈ 0.372** (between 0.370 and 0.375). It is not log2/2 = 0.347 and not log 2 = 0.693. It sits just above log2/2, where prime 2 first enters. Below log2/2 there are no primes, so the knockout equals the full form there.
- Without prime 2, the form fails at support 2L ≈ 0.75 in the odd sector, 0.05 above log 2. This agrees with S12/S16 (pole+arch alone negative near h ≈ 0.38–0.40). In the window ½log2 < L < ½log3 "no 2" is exactly "no primes".
- **The even sector does NOT stay positive.** Hypothesis REFUTED: even fails too once L > 0.41–0.415 without prime 2. The S53 picture ("primes pull the even minimizer down") was about the full-form minimizer; once 2 is removed the minimizer changes.
- **Prime 2 is necessary, but so is every other prime ≤ 13.** Each knockout breaks positivity once 2L passes log p by a small margin. Margin 2L* − log p (upper bounds): p=2 odd ≤ 0.057, even ≤ 0.137; p=3 odd ≤ 0.022, even ≤ 0.052; p=5 ≤ 0.031/0.031; p=7 ≤ 0.024/0.034; p=13 ≤ 0.035. The odd sector always fails first.
- At L = 1.35 every single-prime knockout is indefinite in both sectors. By size: no 2 (−0.057 even, −0.075 odd) > no 3 (−0.023, −0.058) > no 5 (−2.5e−3, −9.2e−3) > no 7 > no 11 (~1e−7) > no 13 (~1e−9). The full form is +1.1e−39 / +2.3e−37 on this basis. So prime 2 is the biggest support, but it is not the only one that is necessary.
- Powers do not matter: removing 2 alone or 2,4,8 gives the same thresholds; same for 3 vs 3,9.
- HEURISTIC reading: the full Weil form on [−L,L] is "barely positive". Its Ritz λ_min falls from ~1e−16 at L=0.4 to ~1e−39 at L=1.35, the super-exponential decay known from the prolate picture (CLASSICAL in spirit, Connes–Consani–Moscovici). The term w_p⟨f, τ_{log p} f⟩ is indefinite as soon as 2L > log p, and it is O(w_p) in size. So any knockout breaks positivity soon after p enters, and the margin shrinks as the form gets more fragile. This makes "every prime is necessary" generic, not special to 2. Removing all powers of p is the explicit formula of ζ(s)(1−p^{−s}), whose extra zeros at s = 2πik/log p lie on Re s = 0, off the critical line (CLASSICAL). So failure is forced at some support; the data say it comes almost at once.
- φ absent. No value here involves the golden ratio.

## Dead ends (H54b)
- REFUTED: "the even sector stays positive (or improves) without prime 2". It fails at L ≈ 0.41.
- REFUTED as a uniqueness claim: "prime 2 is THE necessary prime". All primes ≤ 13 are necessary near their entry points; prime 2 is only the largest by size.

## Files
- `s54_knockout.py` (H54b), `mk_src.py`, `s14_exact.py` (copy of S44/S53, unchanged).
- `lift/`: bundle LP logs/JSONs at L = 1.38, `lift.sh` (B&B pipeline, not run).
- `run1354/`: run scripts, logs, cert JSON (par0), PD logs, block eigenvalues, `atil_sha256.txt`. The 97 MB dumps stay on the box and are not committed.
- `h54b/`: knockout jsonl files and `summary.json`.
