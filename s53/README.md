# Grok S53: full-complex ζ Weil positivity on support 2.7 (L = 27/20) CERTIFIED; the S37 odd failure was a T/β artifact; the odd-floor mechanism

RH OPEN. GRH OPEN. Finite window only. φ absent.

## H53a: odd sector at L = 27/20 with larger T. CERTIFIED
- Same certified lift as S37: `../s37/torus_L1.35_D8x2_p235_p71113_cert.json`, A_cert ≤ 4.5928 (D8×2, two B&B implementations per part).
- The only change is **T = 800** (β = log(T/2π) − 1/T − A = **0.252685**, up from 0.118975 at T=700, about the S38 value 0.2416). K = 880 is scaled with T (K ≈ 0.815·L·T, as in S36–S38); the degree tail at the first omitted degree is ≈1e−282. ng=92, w=1/2, b=2/5, prec 384.
- Assembly: `s53_onestroke_par.py`, which is `s36_onestroke_schurboost.py` with a `--parities` option. The math is identical; it assembles the odd sector only, at about half the cost. Run with `--parities 1 --schur_bs ''` on 5 workers: bounds 143 s, assembly 507 s.
- Bounds: M_S ≤ 2334.01, M0 ≤ 25.87. **rE = 5.31e−87** (quad/entry 6.0e−90), eps_B ≤ 2.1e−134, eps_CDD ≤ 1.6e−278, lamD = 0.25268.
- **ODD: CERTIFIED (Arb, Lean-open, TWO independent PD implementations).** For all odd f ∈ L²[−1.35, 1.35], **Q_ζ(f) ≥ 1.5659e−65 ‖f‖²**.
  - λ_est = +1.5678e−65. At T=700 it was −3.28e−65.
  - impl 1 (`../s36/s36_cholres.py`, 640 bits): λ_min(Ã) ≥ 1.566e−65, ‖E‖ ≤ 3.8e−115, δ_str ≤ 8.8e−112.
  - impl 2 (`../s37/s37_ldlt_exact.py`, S=640, exact integers): λ_min ≥ 1.56599e−65, ‖E‖_F ≤ 9.9e−115.
  - Negative control: l0 = 1.6e−65 (1.02·λ_est) FAILS at pivot 279.
  - The dump sha256 is in `run/atil_sha256.txt`. The 96 MB dump stays on the box and is not committed.
- Even sector at L=27/20 was already CERTIFIED in S37 (T=700, K=780): ≥ 5.2649e−70, also with two PD implementations.
- **Combined: for all complex f ∈ L²[−1.35, 1.35] (support 2.7), Q_ζ(f) ≥ 5.2649e−70 ‖f‖². CERTIFIED** (Arb, Lean-open). The even sector is binding. The two parities decouple: Q(f) = Q(f_e) + Q(f_o) and Q(f) = Q(Re f) + Q(Im f). Each parity bound is rigorous on its own, so using different T per parity is legitimate.
  - The full-complex frontier moves from support 2.6 (S36/S38) to **2.7**. File: `run/cert_L135_FULLCOMPLEX_S53_FINAL.json`.
- **Diagnosis of the S37 failure (NUMERICAL).** `s53_odd_eigvec_diag.py` ran inverse iteration on the S37 T=700 odd dump. Results are in `run/diag_s37_odd_T700.json`.
  - The negative eigenvector is LOW degree: peak at degree 5, 99% of the coefficient mass at degree ≤ 13, high-degree tail ~1e−65.
  - So the deficit is not a K effect. The matrix is A_T = Q_true − ∫_T^∞ (S_lift − β)|f̂|² dt/π: the invisible lift terms integrate to zero over the full line, so they cost only beyond T.
  - The loss is the near-null odd vector's tiny spectral tail beyond T, weighted by S_lift − β. Moving T to 800 removes it, and the odd λ goes from −3.3e−65 to +1.6e−65.
  - So the S37 failure was a **T/β artifact**, not a lift-structure loss. This agrees with S38 (the D8×2 lift is fine for odd at L=1.3).

## H53b: odd-floor mechanism (NUMERICAL)
`s53_odd_floor.py` is the S14/S44 exact Arb Legendre·(1−u²)⁴ assembly split into Pole, Arch and each −w_n T_n. It takes the per-sector minimizer of W_ζ (mpmath eigsy) and decomposes its Rayleigh quotient. Data: M=56 (dps 130, 1024 bits) and M=72 (dps 160, 1280 bits) at L ∈ {1.0, 1.3, 1.35}, in `h53b/floor_M{56,72}.jsonl`.

| M | L | sector | λ_min | pole | arch | Σ primes | n=2 | n=3 |
|---|---|---|---|---|---|---|---|---|
| 72 | 1.0 | even | 5.53e−30 | +0.5699 | −0.5502 | −0.0197 | −0.0196 | −1.0e−4 |
| 72 | 1.0 | odd | 4.06e−27 | −0.0037 | −0.0319 | +0.0356 | +0.0347 | +8.6e−4 |
| 72 | 1.3 | even | 5.06e−44 | +0.6330 | −0.6024 | −0.0306 | −0.0301 | −4.6e−4 |
| 72 | 1.3 | odd | 5.90e−41 | −0.0058 | −0.0498 | +0.0556 | +0.0524 | +3.2e−3 |
| 72 | 1.35 | even | 2.14e−45 | +0.6505 | −0.6168 | −0.0337 | −0.0330 | −6.1e−4 |
| 72 | 1.35 | odd | 1.02e−42 | −0.0060 | −0.0520 | +0.0581 | +0.0545 | +3.5e−3 |

- **The true ζ odd λ_min at L=1.35 is positive on the Ritz basis (NUMERICAL):** 2.26e−37 (M=56) and 1.02e−42 (M=72). It is about 475× the even value. ζ itself is fine; only the certificate route failed at T=700.
- **Mechanism (NUMERICAL, same at M=56 and 72 and at every L).**
  - In the odd sector the pole is negative (−2P², rank one) and the archimedean part sags (≈ −0.05). The **prime powers push the minimizer up** by +0.056–0.058. Of that, **n=2 gives ≈94%** and n=3 ≈6%. n=4 and 5 are ≤1e−4, and 7–13 are ≤1e−10.
  - In the even sector the roles flip: the pole (+0.65) carries it, and the primes pull down (−0.034, mostly n=2).
  - So the odd floor is a prime-2 effect, consistent with S44's prime-2 coupling window.

## Dead ends
- DEAD END (refuted hypothesis): "the S37 odd failure at L=1.35 is lift/structure-driven." It was T/β-driven.
- Do not re-run L=1.35 odd at T=700.

## Files
- Scripts: `s53_onestroke_par.py`, `s53_odd_eigvec_diag.py`, `s53_odd_floor.py`, and `s14_exact.py` (copied unchanged from s44).
- `run/`: logs, cert JSONs, PD outputs, negative control, sha256, diag, and the final combined JSON.
- `h53b/`: the jsonl outputs.
