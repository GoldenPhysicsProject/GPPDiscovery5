# Grok S91 (2026-10-22): H90a, the Connes–Consani near-radical E(prolates) against our Weil window form, plus an audit of Turbo's Schur scalars. RH OPEN.

Status labels: everything below is **HIGH-PRECISION NUMERICAL**. The forms are built with mpmath at dps 90, the prolates and angles in float64 (floor about 1e-13). None of it is interval arithmetic and nothing here is proved.

## Files
- `forms.py` caches the S87 builder (`s87/hp_forms.py`) forms in `forms_cache.pkl` (dps 90). Sign conventions come from s89/calib.py: Q_e = A − P + ppᵀ and Q_o = A − P − qqᵀ.
- `prolate.py` computes the even prolates PS_{2n,0}(c = 2πμ, ·) from the even-Legendre expansion and gives their Fourier signs χ_n. **Validated:** it reproduces CC's tables χ(5.5, n) and χ(7.5, n) to 10 digits, e.g. χ(5.5, 10) = 0.9506583262 and χ(5.5, 11) = 0.5719706153.
- `h90a.py` runs Task A and writes `h90a.json` and `h90a.log`.
- `schur_audit.py` runs Task B and writes `schur_audit.json` and `schur_audit.log`.
- `tables.md` holds the full tables.

## Conventions
- Our window is x ∈ [e^{−L}, e^{L}], with primes n ≤ e^{2L}.
- In Connes–Consani (CC) terms, λ = e^L, μ = λ² = e^{2L}, and their circle length is L_CC = 2L.
- Our L = 0.6, 0.8, 1.0 and 1.2 correspond to μ = 3.32, 4.95, 7.39 and 11.02.
- Our basis is b_j(s) = √((2j+1)/L) P_j(2s/L − 1), for s = log x ∈ [0, L], extended by parity. The coefficient metric is L²(ds).

## Task A: H90a
**Construction (CC 2106.01715 §3, exactly as written there):**
- Prolates: ψ_n(x) = PS_{2n,0}(2πμ, x/λ). The time limit is [−λ, λ] and the band limit is [−λ, λ], with the self-dual Fourier transform e^{−2πixy}.
- ν is the largest n with |χ_n| ≥ 0.9. This gives ν = 5, 9, 14 and 21, close to CC's ν ≈ 2μ − 1. A strict variant uses 1 − |χ| < 1e−6.
- φ_n = ψ_nψ_0(0) − ψ_0ψ_n(0) for even n, which goes to the even sector. φ_n = ψ_nψ_1(0) − ψ_1ψ_n(0) for odd n, which goes to the odd sector. The range is 2 ≤ n ≤ ν.
- g_n(s) = e^{s/2} Σ_{r≥1} φ_n(r e^s) on [0, L], which is CC eq. 3.4.
- Γ-band: the prolate band [−λ, λ] has nothing to do with our Γ-band |t| < 6.29. The only scale that enters is λ = e^L, through the support of f and of f̂. The E-map turns the time-frequency concentration in (x, ξ) ∈ [−λ, λ]² into multiplicative support [λ⁻¹, λ].

**How the angles are computed:** they are taken in L²[0, L] by exact piecewise quadrature, with breaks at the jumps s = L − log r. The E-vectors are therefore *not* truncated to N = 24. At L = 1.2, up to 19% of the norm of an S direction lies outside the polynomial space.

**Results for ζ** (S = span{E(φ_n)}, with dim S = 2, 4, 7 and 10 even and 2, 4, 6 and 10 odd):

| L | sin(R⁻¹p, S) even | sin(R⁻¹q, S) odd | max principal sine, low eigenspace (eig < 1e−8 rel) e / o | min of Q on S (rel) e / o |
|---|---|---|---|---|
| 0.6 | 1.5e−5 | 2.7e−4 | 1.5e−5 / (none below 1e−8) | 4.1e−10 / 1.4e−7 |
| 0.8 | 2.4e−9 | 7.0e−8 | 1.4e−6 / 2.4e−5 | 5.4e−18 / 4.4e−15 |
| 1.0 | 1.1e−11 | 2.9e−10 | 1.2e−5 / 4.5e−5 | 7.2e−23 / 4.2e−20 |
| 1.2 | 1.6e−13 (near the float floor) | 3.5e−12 | 5.3e−6 / 4.7e−5 | 1.1e−26 / 6.5e−24 |

**Controls** (sine for the ζ near-null vector):
- A Gaussian random subspace of the same dimension: 0.74–0.98 (median).
- E applied to the same number of *non-concentrated* prolates (n > ν): 0.98–1.00.
- E applied to random even polynomial f with f(0) = 0: 0.12–0.51.

So both the E-map *and* the prolate concentration are needed.

**Verdict:** H90a HOLDS at the span level (NUMERICAL), with every angle below 1e−3 at every L. This reproduces CC §3 in our normalisation; it is prior art, not new.

**Dead end (sharper version):** the ordered one-to-one match fails at our N. Gram–Schmidt of the first k CC vectors against the k lowest eigenvectors gives a max sine of 1e−3 to 6e−3 at L = 0.6 and 4e−4 to 3e−2 at 0.8 (k ≤ 3), but 1e−2 to 0.6 at L = 1.0–1.2. The likely cause is basis truncation; see the caveat below.

## Fakes (same CC span S, which depends only on λ)
- **F_θ (θ = 0.1):**
  - Q_{F_θ} restricted to S has NEGATIVE eigenvalues at every L: even −1.3e−5, −5.3e−4, −2.6e−3 and −4.8e−3; odd −4.7e−6, −6.6e−4, −2.2e−3 and −3.5e−3.
  - Its near-null vector is off S (sine 1e−3 to 2e−5), and its most negative eigenvector is off S (sine about 1e−2).
  - So the ζ near-radical span *witnesses* F_θ's failure. This SEPARATES F_θ, consistent with the cosh weight breaking the Poisson pairing E(f̂)(x) = E(f)(1/x).
- **DH:**
  - Positive throughout L ≤ 1.2 (q_min 0.14, 6.6e−3, 3.5e−5 and 7.7e−9 even).
  - Q_DH on S is positive and *not* tiny: relative min 4.0e−2, 1.9e−3, 1.2e−5 and 3.3e−8, against 1e−10 to 1e−26 for ζ.
  - DH's ground vector moves toward S as L grows (sine 0.22, 0.043, 3.9e−3 and 3.5e−4).
  - So S separates DH only *quantitatively*: it is far from a near-radical for DH, but DH does not fail on S. That is no failure witness. We did not reach a DH failure L here.
- **H₋ₜ:** NOT RUN this session (time). The s89 delta_ht engine could be plugged in.

## Caveat: basis truncation at N = 24 (and 32)
- At μ = 3 (L = 0.549) our even q_min is 5.55e−8 at N = 32, consistent with CC's "< 6e−8 at μ = 3".
- At μ = 11 (L = 1.2) CC report 2.389e−48. We get 3.7e−26 at N = 24 and 4.9e−32 at N = 32. So at L ≥ 1.0 our low eigenvalues and eigenvectors are Ritz upper bounds and not converged, while for L ≤ 0.8 they are converged.
- This most likely explains why the ordered match degrades at L ≥ 1.0. It also means the S85–S90 q_min values at L ≥ 1 are basis-limited.

## Task B: audit of Turbo's wave-227 Schur scalars (full table in `tables.md`)
Definitions: s_e = 1 + pᵀR⁻¹p (≤ 0 ⇔ even Weil holds), s_o = 1 − qᵀR⁻¹q (≥ 0 ⇔ odd Weil holds), R = A − P, and rel = |s|/(1 + |t|) with t = pᵀR⁻¹p ≈ −1 or qᵀR⁻¹q ≈ 1, so rel ≈ |s|/2.

| L | μ | s_even | s_odd | rel_e | rel_o | q_min/‖Q‖ e / o (N=24) |
|---|---|---|---|---|---|---|
| 0.549 (μ = 3) | 3.00 | −4.25e−8 | 6.69e−4 | 2.1e−8 | 3.3e−4 | 1.2e−8 / 2.7e−6 |
| 0.6 | 3.32 | −1.22e−9 | 2.50e−5 | 6.1e−10 | 1.2e−5 | 3.4e−10 / 1.1e−7 |
| 0.7 | 4.06 | −3.09e−13 | 9.69e−9 | 1.5e−13 | 4.8e−9 | 9.1e−14 / 5.0e−11 |
| 0.8 | 4.95 | −1.21e−17 | 5.87e−13 | 6.1e−18 | 2.9e−13 | 3.9e−18 / 3.3e−15 |
| 0.9 | 6.05 | −3.59e−21 | 5.47e−17 | 1.8e−21 | 2.7e−17 | 1.2e−21 / 3.4e−19 |
| 1.0* | 7.39 | −1.93e−22 | 5.47e−18 | 9.7e−23 | 2.7e−18 | 6.8e−23 / 3.7e−20 |
| 1.1* | 9.03 | −1.39e−24 | 9.19e−22 | 7.0e−25 | 4.6e−22 | 4.9e−25 / 7.5e−24 |
| 1.2* | 11.0 | −2.15e−26 | 6.07e−22 | 1.1e−26 | 3.0e−22 | 1.0e−26 / 6.3e−24 |
| 1.3* | 13.5 | −9.05e−28 | 4.19e−23 | 4.5e−28 | 2.1e−23 | 4.0e−28 / 4.7e−25 |
| 1.4* | 16.4 | −4.20e−29 | 1.79e−24 | 2.1e−29 | 9.0e−25 | 2.0e−29 / 2.7e−26 |

\* = basis-limited at N = 24. At N = 32: L = 1.0 gives rel_e 2.5e−27 and rel_o 6.8e−23, and L = 1.2 gives 1.5e−32 and 1.3e−27. At L ≤ 0.6, N = 24 and N = 32 agree to 3 digits.

**Both blocks are indefinite at every L.** The prime block P has inertia (12, 12) at N = 24 (strongly indefinite). The boundary block A ± pole has 1–3 negative directions. So "two separately indefinite blocks" is qualitatively right in our normalisation.

**Scaling:**
- log₁₀ rel falls by about 2.4 decades per 0.1 in L at N = 24 (fit over L 0.6–1.4). It falls faster when converged: about 4 decades per 0.1 for L ≤ 0.8, and CC-style roughly exp(−c μ).
- The odd |s| is always 10⁴–10⁵ times larger than the even |s|.

**Questions for Turbo (not a verdict):**
1. Their 4.22e−10 (even) matches our rel at L ≈ 0.60 (μ ≈ 3.35), but their 1.24e−11 (odd) matches ours at L ≈ 0.76 (μ ≈ 4.6). No single L in our normalisation gives both. Their odd < even ordering is the reverse of ours.
2. "L = 3" read as our L = 3 or as L_CC = 3 (our 1.5) would put the true cancellation far below 1e−30 (8e−30 at L = 1.5 already at N = 24, and smaller when converged). A value of 1e−9 there would then be a precision or basis floor.
3. Read as μ = λ² = 3 (our L = 0.549), we get rel_e 2.1e−8 and rel_o 3.3e−4: the even sector is within a factor of 50, the odd far off. Read as λ = 3 (our L = 1.099), we get 7e−25 (N = 24) and 4e−30 (N = 32).
4. Sign: their S_even > 0, whereas our s_even must be ≤ 0 for even positivity. Which sign convention and which normalisation of the pole vector and basis do they use, and at what basis size and precision?

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
