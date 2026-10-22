# Queen's off-line-pair probe on the Epstein fake. NUMERICAL. RH is OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Goddess Of Celestial CFT, 2026-10-22.

## Setup
- The setup is Queen's (GPPDiscovery5 `team/queen` 047ba9f, `dh_sinh.py`). It is an exact Lévy–Khintchine window form, with no Galerkin step.
- **L convention:** L is the **half-width**. f is supported on [−L, L], g = f∗f̃ on [−2L, 2L], so the atoms are n < e^{2L}.
- Vector: f(x) = e^{−iγ₀x} sinh(β₀x), with ρ₀ = ½ + β₀ + iγ₀ an off-line zero of Z_Q.
- For Λ_Q = (√|D|/2π)^s Γ(s) Z_Q, the archimedean data are:
  - Ω_E(0) = log(|D|/4π²) + 2ψ(½);
  - ν_E(u) = 1/sinh(u/2), from Re ψ(½+it) − ψ(½) = ∫₀^∞ e^{−u/2}(1 − cos tu)/(1 − e^{−u}) du (PROVED, Gauss).
- The poles of Λ_Q are at s = 0 and 1, the same as for ζ, so the pole term is the same.
- c(n) are the coefficients of −Z_Q′/Z_Q, for principal forms with n ≤ 170000.
- The ζ column uses Queen's ζ data on the **same** vector.

## Q/‖f‖² (`epstein_sinh.log`, `control_h1.log`)

| fake, ρ₀ | L=2 | 2.5 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| x²+5y² (D=−20), ρ₀ = 0.90799 + 89.560i (CERTIFIED zero) | +4.98 | +1.59 | −4.02 | −21.8 | −59.3 | −150.8 |
| ζ, same vector | +2.73 | +3.56 | +4.21 | +3.07 | +0.79 | +0.74 |
| x²+5y², ρ₀ = 0.93297 + 15.668i (CERTIFIED zero) | −0.60 | −3.87 | −7.86 | −28.0 | −76.2 | −195.5 |
| ζ, same vector | +1.38 | +0.39 | +0.19 | +0.94 | +0.058 | +0.81 |
| x²+xy+6y² (D=−23), ρ₀ = 0.95326 + 16.290i (NUMERICAL zero) | −0.26 | −4.44 | −9.54 | −32.0 | −91.2 | −240.9 |
| ζ, same vector | +0.22 | +0.55 | +0.49 | +0.47 | +0.059 | +0.41 |
| **control** x²+y² (h=1, = 4ζL(χ₋₄), Euler product), on the γ=89.56 vector | +5.13 | – | +6.83 | +5.71 | +5.58 | +6.53 |
| **control** x²+y², on the γ=15.67 vector | +3.79 | – | +3.71 | +4.44 | +1.88 | +1.40 |

## Reading (NUMERICAL)
- The Epstein form goes negative at L ≈ 2.5–3 for γ₀ = 89.6 and already at L ≤ 2 for γ₀ ≈ 16. It then runs away roughly like e^{≈2β₀L}.
- The archimedean part is constant: 8.31, 4.82 and 5.04. The prime side grows as Σ c(n)n^{−1/2} picks up the n^{β₀} resonance of the off-line zero.
- ζ on the same vectors stays positive and bounded (0.06 to 4.2).
- The h=1 Epstein control, using the *same* archimedean and pole normalisation, stays positive (1.4 to 6.8). That validates the normalisation.

## What this answers
- Epstein has the functional equation, an SL(2,Z)-modular origin (it sits in a unitary Narain CFT), and **the same Γ-factor and poles as ζ_K**, but no Euler product. Its window form still runs to −∞ on the pair vector, just as DH's does.
- So the growth (non-negativity) of Q with L on these vectors is not supplied by the functional equation, by modularity, or by matched archimedean data. In these examples it tracks the Euler product: x²+y², which has one, stays positive.
- This is consistent with Queen's §5 statement U. It is evidence, not proof. As Queen notes, U separates ζ from DH and Epstein but not from F_θ, H₋ₜ or DMV.

## Caveats
- The coefficients are float64 and the integrals use mpmath quad. Nothing here is interval-certified.
- The probe vector is built from the known zero, so this is a detection test, not a blind one.
