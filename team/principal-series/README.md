# Principal series (team seat): code and data

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN.** Nothing here proves or implies RH.

Notes live in GPP-bridge `work/grok/team/principal-series.md`, on branch `team/principal-series`.

## twist/: the GL(1) twist calculus for the window Weil form

**PROVED (elementary; consequence of Weil 1952 and Burnol math/9809119).** Let W_ζ be the window Weil kernel in S14 normalisation, and Q_ζ(f) = <W_ζ, f*f~> on [−L, L]. Then:
- G_θ' = Λ(s+iθ')Λ(s−iθ') has kernel 2cos(θ'r) W_ζ, so Q_G(f) = 2Q_ζ(f cos θ'u) + 2Q_ζ(f sin θ'u).
- F_θ = Λ(s−θ)Λ(s+θ) has kernel 2cosh(θr) W_ζ, so Q_F(f) = 2Q_ζ(f cosh θu) − 2Q_ζ(f sinh θu). This is a Krein J-form.

The two consequences tested here:
- **F_θ onset law.** F_θ's even form first goes negative where λ_min^even(Q_ζ)(L) = θ² κ(L), with κ = Q_ζ^odd(x·v) and v the lowest even eigenvector of Q_ζ.
- **G_θ' inheritance.** G_θ' inherits every ζ window certificate (both parities), with constant 2·min(c_even, c_odd).

| file | what | label |
|---|---|---|
| `twist.py`, `twist_results.json` | float64 builder (copy of the s82os conventions) for ζ, F_θ and G_θ'; λ_min and failing-vector overlaps, L = 0.3–0.7 | NUMERICAL |
| `onset.py`, `onset_results.json` | float64 onset: direct bisection vs. predicted | NUMERICAL |
| `hp_onset.py` (+ `hp_onset_lib.py`), `hp_onset_N24.json/.log`, `hp_onset_N32.log` | mpmath dps 80, S87 exact-polynomial builder `s87/hp_forms.py` (imported unchanged), N = 24 and 32 | HIGH-PRECISION NUMERICAL (not interval) |
| `hp_conv.py`, `hp_conv_N40.*` | basis-convergence check at N = 40 | HIGH-PRECISION NUMERICAL |

**Convergence warning.** λ_min^even is converged to 2–3 digits for L ≤ 0.8 (N = 24, 32 and 40 agree). For L ≥ 0.9 it is **not** converged: at L = 1.2, N = 24 gives 3.7e−26, N = 32 gives 4.9e−32 and N = 40 gives 2.5e−36. Do not quote q_min values above L ≈ 0.85 from these runs. Onsets for θ ≤ 1e−8 fall in that range and are therefore not quoted.
