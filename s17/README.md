# S17 — Audit of Codex's uniform all-Dirichlet small-window Weil gap, plus χ3/χ4/χ5 vs ζ/F_{1/4} bounds

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN. Dirichlet GRH OPEN.** Nothing here
tests zero location: every window is prime-free, or holds only the first prime powers.

Form (Codex normalization, real g ∈ C_c^∞((−h,h))): Q_η = 1_{η=1}P + E_ε + (log f − c_ε)‖g‖² − W_η, where
E_ε = ∫₀^∞ ν_ε(t)(F(0)−F(t))dt, ν_ε = 2e^{−(1+2ε)t/2}/(1−e^{−2t}), c_ε = log π − ψ((1+2ε)/4), P = 2M(½)M(−½), and
W_η = 2Σ Λ(n) Re η(n) n^{−1/2} F(log n). F_{1/4} uses the S4 convention (zeros ρ±θ, summed over σ = ±θ), so it counts twice
as many zeros. Basis: Legendre·(1−u²)⁴ (the S14 exact rational autocorrelations G_jk are reused), with Arb moments.

## Files and status labels
| file | what | label |
|---|---|---|
| s17_closed_forms.py → out_closed_forms.txt | T0/T1 closed forms vs quadrature (err ≤ 3e−50); ψ(3/4)−ψ(1/4)=π; c0,c1; B(h) balls; roots | NUMERICAL (dps 50) + Arb balls |
| s17_selftest.py → out_selftest.txt | S17 ζ matrix = S14 independent arch decomposition (diff ≤ 4e−71) | CERTIFIED consistency |
| s17_zero_check.py → out_zero_check.txt | form vs zero side Σ over zeros of L(s,χ) for χ3, χ4, χ5-quadratic and ζ (T<200), h=0.3, 0.5 | NUMERICAL (rel. diff 1e−9 to 1e−13) |
| s17_bounds.py → out_bounds.txt | rigorous prime-free lower bounds: Codex / Astra chains / S17 bathtub / S17 combined | CERTIFIED (Arb), elementary proof in the note |
| s17_ritz.py → ritz_M{16,24,32,48}.json, TABLE.md | Ritz λ_min per g-parity; Arb Rayleigh quotient of a rational g (upper); rigorous lower | NUMERICAL Ritz + CERTIFIED upper and lower |
| s17_window_ends.py → out_window_ends.txt | Arb lower bound at each character's window end | CERTIFIED |
| s17_onset.py, s17_onset_cert.py | where Arch, Arch+Pole and the primes become binding, by sector | NUMERICAL + CERTIFIED negative witnesses |

## Headline
- Codex's theorem (every primitive η, h ≤ 1/25: Q_η > 0.01‖g‖²): **CONFIRMED**. B(1/25) = 0.0707585612 and B(1/32) = 0.3263222785 (Arb). The theorem can be strengthened.
- Normalization (+log f, c_ε, the pole only for η=1, Re η(n) primes): **matches the zero side** of L(s,χ3), L(s,χ4), L(s,χ5) and ζ to between 1e−9 and 1e−13. No mismatch.
- **S17 bathtub + Schur bound (PROVED, elementary):** E_ε ≥ T_ε(h)‖g‖², with the tail taken from h instead of 2h. Combined with Astra's numerical-radius chains this gives a uniform all-η gap: **Q_η ≥ 0.91664‖g‖² at h = 1/25** (Codex's bound there is 0.0708), and **Q_η > 0.01‖g‖² for all h ≤ 0.0962850**.
- Prime-free character windows, certified positive (real g): χ3 on all of h < ½log2; χ5-quadratic on all of h < ½log2; χ4 for h ≤ 0.53463 (its window runs to ½log3); χ5-quartic on all of h < log2 (Re χ(2) = Re χ(3) = 0).
- Two-face mechanism (Ritz, plus Arb witnesses): on even g the pole pays the gamma deficit for ζ, and for nontrivial η the conductor log f takes over that job. On odd g, Arch alone is positive while the window is prime-free. Arch+Pole on odd g goes negative by h = 0.38 (Arb: −0.0038049), after p=2 enters, and from then on the prime current must repay it. That is the H16 odd floor.
- F_{1/4}: certified negative witness at h = 0.45, even g (Arb −2.19197e−5, M=48). ζ, χ3, χ4 and χ5 stay Ritz-positive there.
- φ: does not appear. Only log 5 enters for q = 5, and root numbers drop out of the log-derivative.
