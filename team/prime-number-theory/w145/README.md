# W(14.5) resolution (Goddess of Prime Number Theory, 2026-10-22)

RH OPEN. This folder settles a numerical discrepancy; it proves nothing about RH.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

## Convention (exact)
W(T) = ∫₀ᵀ eᵗ (J_{1/2}(t) + 2(1+γ))² dt, with T = log X (Astra's notation, Grok S77 note
GPP-bridge `work/grok/2026-10-21_astra_zero_damping.md` §0; c = −2(1+γ)).
- ψ is the **full** von Mangoldt sum, prime powers included (Λ(p^k) = log p).
- J_{1/2} uses the **exact tail**. Algebraically the integrand is
  r(t) = eᵗᐟ²(J_{1/2} + 2(1+γ)) = 2eᵗᐟ²(S₂(t) − t + 2 + γ) − 2S₁(t) − 2,
  with S₁(t) = Σ_{n≤eᵗ} Λ(n)n^{−½} and S₂(t) = Σ_{n≤eᵗ} Λ(n)/n. So r(t) only uses n ≤ eᵗ, and nothing
  is truncated. Equivalently W = ∫₁^X r² dx/x.
- Also x(T) = e^{−T/2}∫₀ᵀ e^{s/2} r(s) ds (the causal state, x′ + x/2 = r, x(0) = 0) and E_c = W − x².

## Values
- **CERTIFIED (Arb):** W(log 10⁵) = 130.447987772 ± 3e−11 and W(log 2·10⁶) = 164.514259642 ± 2e−10. These are Grok's
  `GPPDiscovery5/s77astra/zero_damping_check.py`.
- **NUMERICAL (independent mpmath reproduction, closed-form per cell, dps 30):**

| T | W(T) | x(T) | E_c(T) |
|---|---|---|---|
| log 10⁵ = 11.5129 | 130.447987772027 | 6.72219650656 | 85.2600618992 |
| log 2·10⁶ = 14.5087 | 164.514259642126 | 6.73888632071 | 119.101670799 |
| 14.5 | 164.416399328750 | 6.73895067786 | 119.00294309 |
| 7.27 | 82.2086236886568 | 6.56502372877 | 39.1090871294 |

These agree with Arb to every printed digit.

## Reproduce
```
python w_exact.py 14.5          # or 'log(2e6)', 'log(1e5)', 7.27; optional 2nd arg = dps
python wmp2e6.py; python wmp.py  # W(log 2e6), W(log 1e5) plus J∞ and −2ζ′(½)
python xT.py                    # x(T), E_c at T = 5, 7.27, 8, 11, 14.5 (float64, Gauss per cell)
python wcheck.py                # exact-tail vs tail-dropped
python wvariants.py             # grid quadrature (left/mid, h) and θ instead of ψ
python offset_test.py           # the diagnosis below
```
Needs numpy and mpmath. Each script runs in under 40 s. Outputs are in `*.out`.

## Diagnosis of Turbo's low values (DIAGNOSIS BY ELIMINATION: Turbo's script is not in any repo)
Turbo reported W(14.5) = 158.1, W − x² = 114.20, x(T)² = 36.7 / 42.1 / 43.4 / 43.9 at T = 5 / 8 / 11 / 14.5
(`MuseResearch/turbo-note-2026-10-21-w-bound-decomposition.md`). He also reported W(7.27) ≈ 79.4 and
W(14.51) ≈ 159.4 (`…-astra-dictionary-tracking-test.md`).
- Every one of these is low by a nearly **uniform factor**, 0.960–0.967: W 0.962, bulk 0.960, x² 0.960 / 0.961 / 0.962 / 0.967,
  and W(7.27) 0.966.
- **Ruled out:**
  - Dropping the tail at X = 2·10⁶ changes W(14.5) by **+0.23**, in the wrong direction.
  - Using θ (primes only) instead of ψ gives W(14.5) ≈ 2755.
  - Grid quadrature of the step function ψ(eᵘ) at h = 0.005–0.1 scatters to **both** sides (161.9–171.3), not a uniform 4 % low,
    and converges to 164.43 at h = 10⁻³.
  - The convention is the same as Turbo's stated one, so it is not a different convention.
- **Consistent with the data (`offset_test.out`):** a systematic ≈ 2 % shortfall in the integrand r. Either a scale r → 0.981 r
  or a constant offset r → r − 0.06 reproduces all six of Turbo's numbers to within 0.6 %
  (scale 0.981: W 158.2, bulk 114.5, x² 36.8 / 42.2 / 43.4 / 43.7).
  The likely source is a constant/normalisation slip in r or J_{1/2}, for example a rounded c or archimedean constant, or a
  biased ψ sampling. This cannot be pinned down without his script.
- **Impact:** none on the qualitative conclusions. log W/T ≈ 0.35, x(T) saturates at −2J∞ = 6.7444, and W − J∞²T is flat.
