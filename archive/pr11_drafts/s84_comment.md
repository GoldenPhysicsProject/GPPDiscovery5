**Grok S84 (2026-10-22): ζ's near-null direction at support 2.84 is a nodeless pole bump that vanishes at every low zero; DH fails in the same direction (overlap 0.85), through its first off-line zero. RH OPEN.**

**Comms read** (since S83, 6331392543)
- No new PR #11 comments and no new bridge commits from other seats, so nothing is aimed at Grok.
- Still pending: a teammate replay of the 2.84 certificates, and Turbo's W(14.5) convention question.

**Results** (note `work/grok/2026-10-22_s84_nearnull_zero_channels.md`, scripts GPPDiscovery5 `s84/`)

*Method.* Q(f) = Σ_ρ f̂(z_ρ)² for real even f, and f̂ has a closed form in the cosine basis. On the zero side this makes the window form the frame operator of the zero exponentials:
- an on-line zero contributes a projector (a principal-series state);
- an off-line pair contributes a (1,1) indefinite block.

This much is a PROVED restatement of Weil.

1. **CERTIFIED** Arb Ritz upper bounds:
   - ζ at L = 71/50: even λ₁ ≤ 4.44455e−81 (λ₂ 6.6e−73), odd λ₁ ≤ 7.13734e−77.
   - ζ at 137/80: even ≤ 2.62874e−149 (K160, not converged in K).
   - DH at 71/50: ≤ 3.26e−15.
2. **NUMERICAL**, ζ at 71/50:
   - The explicit formula with 8000 zeros reproduces the Arb RQ (ratio 1.0006 even, 1.0004 odd).
   - The near-null vector is nodeless, a low-frequency bump (87 % of the f̂ mass in |t| < 5, pole term 1.52‖f‖²), and stable from K96 to K128.
   - **The lowest zeros carry nothing.** The shares of λ_min are 5e−78 from γ₁, 8e−75 from γ₂, 7e−73 from γ₃, and 4e−32 from all γ ≤ 2πe^{2L} = 107.5. The residual sits at γ ≈ 210–235.
   - So the low zeros are the near-null vector's spectral nodes (|f̂(γ₁)|²/‖f‖² = 1.2e−158).
3. **Correction (NUMERICAL).** The S75/S82-OS ζ vector at 137/80 was dps-60 noise: it has overlap 0.013 with the true one, and its "strongly sign-changing" shape was a float artefact. The true near-null vector is nodeless there too.
4. **NUMERICAL**, same direction as DH:
   - The ζ near-null vector has overlap **0.845** with the certified DH witness at 137/80, and 0.843 with DH's lowest vector at 71/50.
   - Q_ζ(f_DH) from the zeros is 1.500643011e−8, which matches the certified Arb value to 10 digits. γ₁ carries 99.77 % of it.
   - DH zero side: the off-line zero 0.8085 + 85.699i gives −1.234e−28, against +1.19e−28 from the on-line zeros (mostly #42–46, flanking it). With zeros up to T = 520 the total is −1.30e−30 against the certified −2.07e−30, resolved only to ±3e−30: off-line zeros above 177 are missing, plus the tail. The split itself is 40× larger than that uncertainty.
   - So the geometric "pole channel" and the zero-side "first off-line pair" describe the same failure.
5. **Fakes ruled out: none** (F_θ, DH, H_{−t}). H82-OS-a is a DEAD END as a separator, because ζ's own weakest direction is nodeless. This explains why the S83 toy failed nodelessly.
6. **HYPOTHESIS H84a.** RH(L) is a zero-frame inequality in the pole-bump direction: a mirror pair is harmless only if its indefinite block is invisible there. This is a mechanism, not a proof.

**Certificates (secondary line)**
- **CERTIFIED.** The batched-Arb rebuild of the even support-2.84 dump (sha256 dada6eb4…, cells 1323 s against 4323 s in S71) gives:
  - s36 Cholesky: λ_min ≥ 1.071000e−81;
  - s37 exact LDLᵀ: λ_min ≥ 1.0711e−81;
  - control at 1.005·λ_est: FAILS (pivot 1183).

  λ_est 1.07214622023e−81 matches S71 in every printed digit.
- **CHECKED.** The entrywise compare "PASS=false" (max |Δ| 2.78e−144, 3.6× the radius-only tolerance) is the wrong test, not a bug. X and Y are point midpoint approximations, and Δ is ~1e−8 of the inv_err budget of 2.05e−136.
- Next: out-of-core strips at Np ≈ 9.6k for ζ at L ≈ 1.7105 (ESTIMATE ~3.5 days on 8 cores).

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
