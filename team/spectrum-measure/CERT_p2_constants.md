# CERT_p2_constants — odd archimedean H-margin at the p = 2 entry (Arb-certified), Celestial Conj C.13.H certified, (K2) odd direct

Spectrum & Measure seat, GPP RH team. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
**RH is OPEN.** Numerics are evidence, not proof; items marked CERTIFIED are rigorous computer-assisted inequalities
(Arb ball arithmetic, python-flint 0.9.0, scripts and logs in `cert/`). eval. = evaluation of an explicit closed form.

Inputs used, as stated (not re-derived): Prime §N.12–N.13 (team/prime-number-theory 611e39a): N.13.0 notation, S1 = N.12.2
(|2c₂B| ≤ c₂‖ξ‖²_{H(E)}), (Ω2) Ω increasing in |t|, (Ω5); parity split N.12.3. Celestial §C.13.B.5 (team/celestial-cft a7f943a):
Conj C.13.H, Props C.13.H′, C.13.H″; §C.14 Conj C.14.K, Cor C.14.D. Our SECTOR_COERCIVITY_p2 (f1df9ed).
Notation: a = log 2, c₂ = a/√2 = 0.4901290717…, L₀ = a/2, 2L = a + ℓ, ‖u‖²_H = ∫(u′² + ¼u²) = (1/2π)∫(¼+t²)|û|²,
Ω(t) = 2θ′(t) = Re ψ(¼ + it/2) − log π, R(ξ) = (1/2π)∫(¼+t²)Ω|ξ̂|², Q_L = R − 2c₂B (N.13.1; the form is L-independent).

> **Grok audit (S103g, 2026-10-28): SOUND.**
> - *Replay.* `cert/cert_odd_margin.py` was re-run at 80 and 128 bits with N = 2¹⁸ (python-flint 0.9.0). Every printed bound matches the committed log exactly.
> - *Rigour of the chain.* All quantities are arb balls, and inputs are parsed from decimal strings into balls. D(T) is a rigorous upper sum: ball evaluation of dens_o over boxes covering each subinterval, with upper endpoints only. 1 − D enters through its lower endpoint.
> - *No truncation remainder is needed.* ω̃ is nondecreasing on all of [0, ∞), so the step minorant beyond Tmax = 64 is the constant ω̃(64).
> - *Monotonicity is analytic, not read off the grid.* It is (Ω2), PROVED via the series Σ y²/(c_n(c_n² + y²)). The ball increments are used only as enclosures of true increments, which are ≥ 0.
> - *Lemmas A and B are checked line by line.* For Lemma A, the Riesz/Green step in H¹₀ holds with no boundary term, and the closed form dens_o was re-checked against direct quadrature of ∫G s_t·s_t (30 digits; t = 0.3, 4.53, 10, 30). Lemma B is the standard layer-cake bound. T is a grid point, and min(Ω − lo, Ω(T) − hi) is nondecreasing.
> - *Consistency.* A 24-mode Galerkin Ritz value gives inf R/‖·‖²_H ≤ 0.548 (NUMERICAL), which is compatible with the certified ≥ 0.46255.
> - *Inputs.* Thms 1–2 use only Ω (digamma) plus Lemmas A and B. Thm 3 adds Prime S1 (N.12.2) and (Ω5), both PROVED on main and Grok-audited.
> - *"CC-free", made precise.* No CC numerics (ε, q, ε′(1⁺)), no CC §5/Thm 6.11, no κ^odd ≤ 1, and no Mori. The only CC link is the exact form Q_L = R − 2c₂B (main G.4h "One symbol", from CC (102)), which is the classical Weil explicit formula on a window 2L < log 3.
> - *Scope of Thm 3.* It covers only the sub-window 2L ∈ (log 2, log 2 + 0.028] of the p = 2 step (log 2, log 3], and only the odd sector.
> - *Label fixes.* The Cor (R6) values were rounded up. They should be 0.0029 (T = 12.875) and 0.0050 (T = 16); recomputed, ℓ_H = 0.002974 and 0.005061. Separately, C.14.A's pointwise σ_S condition is **not** established here; σ_S stays a GAP. The odd sector no longer needs it.

## 0. Summary (labels)

| # | Statement | Label |
|---|---|---|
| R1 | Lemma A: sharp pointwise uncertainty in H for odd/even ξ ∈ H¹₀(−L, L), closed-form density | PROVED |
| R2 | Lemma B: bathtub bound for increasing symbols against an H-energy CDF majorant | PROVED |
| R3 | **m₀^odd(L₀) := inf_{odd ξ ∈ H¹₀(−L₀,L₀)} R(ξ)/‖ξ‖²_H ≥ 0.46255** | **CERTIFIED** (no CC input) |
| R4 | **Celestial Conj C.13.H holds** (odd, L₀) with (T, μ_hi, μ_lo) = (12.875, 0.6, 0.25) and (16, 0.7, 0.3) | **CERTIFIED** (no CC input; κ^odd ≤ 1 not used) |
| R5 | **Odd (K2) entry, direct:** Q_L(ξ) ≥ 0.0126‖ξ‖²_H for odd ξ ∈ H¹₀(−L, L), all 2L ≤ log 2 + 0.028 | **CERTIFIED** (given Prime S1, (Ω5), proved) |
| R6 | Via C.13.H′ with R4: (K2) odd for ℓ ≤ 0.0050 (eval.); R5 is the stronger range | PROVED given R4 / eval. |
| R7 | Even sector by the same method: bound ≈ 0 at ℓ = 0 (inconclusive) | eval., NUMERICAL |
| R8 | σ_S thresholds t₀(B) for |t| ≥ t₀: σ_S ≥ c₂ | eval. given B; **B GAP (uncertified)** |
| R9 | δ_o⁰, k_∞ for the odd sector: no longer needed (R3/R5 replace SECTOR_COERCIVITY_p2 Thm 2.1-odd's δ_o route) | PROVED (logic) |

GAPs that remain: **GAP B** (2ε′(1⁺) + ‖q‖₁ not certified: needs CC's ε, q, eqs. (102)/(104), not attempted here);
σ_S-route (a),(b) of the C.13.H fold-in therefore not certified, but they are not needed (R4). **Even sector: m₀^e and (P) OPEN** (item 3).

## 1. Lemma A (sharp pointwise uncertainty in H) — PROVED

For L > 0, t ∈ ℝ, odd ξ ∈ H¹₀(−L, L): |ξ̂(t)|² ≤ K_o(t)‖ξ‖²_H, with K_o(t) = ⟨G s_t, s_t⟩, s_t(x) = sin(tx) on (−L, L),
G = (−d² + ¼)^{-1} with Dirichlet conditions at ±L; equality is attained (by u = G s_t).
*Proof.* ξ̂(t) = −i∫ξ s_t for odd real ξ; Riesz in H¹₀ with inner product ∫(u′v′ + ¼uv): ∫ξ s_t = ⟨ξ, G s_t⟩_H, so
|ξ̂|² ≤ ‖ξ‖²_H ‖G s_t‖²_H = ‖ξ‖²_H ⟨G s_t, s_t⟩. ∎ Solving (−d² + ¼)u = sin(t·), u(±L) = 0:
u = sin(tx)/(t²+¼) − A sinh(x/2), A = sin(tL)/((t²+¼) sinh(L/2)), and integrating,
> dens_o(t) := (1/π)(¼+t²)K_o(t) = (1/π)[ L − L·sinc(2tL) − 2 sin(tL)(½coth(L/2) sin(tL) − t cos(tL))/(t²+¼) ],

> dens_e(t) = (1/π)[ L + sin(2tL)/(2t) − 2cos(tL)(½tanh(L/2)cos(tL) + t sin(tL))/(t²+¼) ] (even ξ, cos(t·)).

(Cross-checked numerically against the Dirichlet-mode series Σ_n |⟨s_t, e_n⟩|²/(k_n²+¼), e_n = √(2/2L) sin(nπ(x+L)/2L): agreement to 1e−10.)
Consequence: the H-energy measure dν = (1/π)(¼+t²)|ξ̂|² dt on t ≥ 0 (total mass ‖ξ‖²_H) obeys ν[0, T] ≤ D(T)‖ξ‖²_H, D(T) = ∫₀^T dens.

## 2. Lemma B (bathtub) — PROVED

If ω̃ is nondecreasing on [0, ∞), ν ≥ 0 has mass 1 and ν[0, s] ≤ D(s), then for any grid 0 = e₀ < … < e_N:
∫ω̃ dν ≥ ω̃(e₀) + Σ_i (ω̃(e_{i+1}) − ω̃(e_i))·max(0, 1 − D(e_{i+1})).
*Proof.* The left-endpoint step function ω_step ≤ ω̃ is nondecreasing; ∫ω_step dν = ω̃(e₀) + Σ_i Δ_i·ν[e_{i+1}, ∞) and ν[e_{i+1},∞) ≥ 1 − D(e_{i+1}), Δ_i ≥ 0. ∎
Applied to a band symbol ω (ω = Ω − lo on [0,T), Ω − hi on [T,∞)) through the increasing minorant
ω̃ = min(Ω(t) − lo, Ω(T) − hi) on [0,T), ω̃ = Ω − hi on [T, ∞) (increasing by (Ω2)).

## 3. Certificates (cert/cert_odd_margin.py, log cert/cert_odd_margin.log) — CERTIFIED

Implementation: [0, 64] split into N = 2¹⁸ dyadic subintervals (all T used are grid points); D upper sums from ball
evaluation of dens_o on each subinterval box (sup over the box); Ω by acb digamma; every ball checked finite
(non-finite ⇒ abort: the "bessel lesson"); overlapping balls in the min handled by ball union (never triggered).
Run at 80 and 128 bits (bits ≫ log₂(1/margin) + 13); both identical to the digits shown. The printed ball encloses the computed bound for that grid, i.e. a rigorous lower bound for the infimum.

| run | quantity bounded below (per ‖ξ‖²_H = 1) | certified lower bound |
|---|---|---|
| ℓ = 0 | m₀^odd(L₀) = inf R | **0.4625587612853** |
| ℓ = 0.025, T = 23, η = 1 | inf [R − c₂(1+1/η)(2ℓT/π)ν_{<T} − c₂(1+η)ν_{≥T}] | **0.0461355147** |
| ℓ = 0.028, T = 22, η = 1 | same | **0.0126078043** |
| C.13.H, T = 12.875 | inf [R − 0.6ν_{≥T} − 0.25ν_{<T}] | **0.0491093609** |
| C.13.H, T = 16 | inf [R − 0.7ν_{≥T} − 0.3ν_{<T}] | **0.0758975751** |

Also certified in the log: c₂ = 0.49012907173…, Ω(12.875) = 0.71716…, Ω(16) = 0.93455… (so μ_hi < Ω(T), as must be).
Sanity (eval.): the test vector sin(πx/L₀) has R/‖·‖²_H = 0.609 ≥ 0.4626.

## 4. Theorems

**Thm 1 (R3; CERTIFIED).** For odd ξ ∈ H¹₀(−L₀, L₀): R(ξ) ≥ 0.46255‖ξ‖²_H. Since R is L-independent and the domain grows
with L, m₀^odd(L) is nonincreasing in L; the value at L₀ is the Prime N.13.5 deliverable 1 "explicit m₀^odd > 0" in the H-norm.
Note Thm N.13.A with m₀ = 0.4625 gives only ℓ* ≈ 5×10⁻²⁴; Thm 3 below is the useful form.

**Thm 2 (R4; CERTIFIED). Celestial Conj C.13.H holds** at L₀ in the odd sector: for every odd real ξ ∈ H¹₀(−L₀, L₀),
> R(ξ) ≥ 0.6‖ξ^{>T}‖²_H + 0.25‖ξ^{<T}‖²_H  (T = 12.875), and R(ξ) ≥ 0.7‖ξ^{>T}‖²_H + 0.3‖ξ^{<T}‖²_H  (T = 16),

with μ_hi > c₂ and μ_lo > 0 as required (even with slack 0.049‖ξ‖²_H, resp. 0.076‖ξ‖²_H).
*Proof.* R − μ_hi ν_{≥T} − μ_lo ν_{<T} = ∫ω dν exactly (R is a Fourier multiplier: no cross terms between bands);
Lemmas A, B and the certificate. ∎ **It uses neither CC §5 / Thm 6.11, nor σ_S, nor κ^odd_{L₀} ≤ 1.** So C.13.H is
CERTIFIED unconditionally (relative to the PROVED Ω facts), not merely CERTIFIED-conditional-on-CC.

**Cor (R6; PROVED given Thm 2, eval.).** Celestial Prop C.13.H′ then gives (K2) odd for ℓ ≤ ℓ_H = 0.0029 (T = 12.875) and 0.0050 (T = 16) *(Grok audit S103g: rounded down; exact values 0.002974 and 0.005061)*.

**Thm 3 (R5; CERTIFIED). Odd (K2), direct, polynomial range.** For 0 ≤ ℓ ≤ 0.028 (μ = e^{2L} ≤ 2·e^{0.028} ≈ 2.057) and odd ξ ∈ H¹₀(−L, L):
> Q_L(ξ) ≥ 0.0126‖ξ‖²_H  (≥ 0.0461‖ξ‖²_H for ℓ ≤ 0.025).

*Proof.* S1: |2c₂B| ≤ c₂‖ξ‖²_{H(E)} ≤ c₂[(1+1/η)‖ξ^{<T}‖²_{H(E)} + (1+η)‖ξ^{≥T}‖²_H] and (Ω5) ‖ξ^{<T}‖²_{H(E)} ≤ (2ℓT/π)‖ξ^{<T}‖²_H.
So Q_L ≥ ∫ω dν with lo = c₂(1+1/η)(2ℓT/π), hi = c₂(1+η); Lemma A at L = (a+ℓ)/2 (density at that L), Lemma B, certificate.
It is monotone in ℓ (smaller ℓ: smaller domain and smaller lo), so the ℓ = 0.028 run covers [0, 0.028]. ∎
This range is well inside the CC cap 2L₀ < log 2.44 (NUMERICAL) and 5× the C.13.H′ range; float exploration shows the method stops near ℓ ≈ 0.029.

## 5. C.13.H fold-in, items (a)–(c), plainly

- (a) certified inf σ_S on |t| < T: **not done.** σ_S = Ω + ε̂ needs CC's ε (GAP B). **Not needed**: Thm 2 certifies μ_lo directly.
- (b) analytic/certified ε̂ bound for |t| ≥ T: σ_S ≥ Ω − B/(¼+t²) is PROVED given B; thresholds (eval. given B) for σ_S ≥ c₂:
  B = 46/50/60/70/80/100 ⇒ t₀ = 13.30/13.50/13.96/14.38/14.79/15.53. ε′(1⁺) ≈ 22.9965 (NUMERICAL). B itself: **GAP B**.
- (c) certified m₀^o: given by Thm 1 directly (0.46255), not through C.13.H″'s inf formula.
- **Verdict:** C.13.H is CERTIFIED (odd, L₀), with no dependence on CC; **κ^odd_{L₀} ≤ 1 is not needed** for the odd sector.
  The odd (K2) entry is CERTIFIED for ℓ ≤ 0.028 (Thm 3).
- Cross-check: Ω = c₂ at 10.2615 (eval., our earlier 10.63 Stirling threshold is not sharp), inside Celestial's Arb bracket
  (10.25, 10.30); Ω = 2c₂ at 16.748 inside their (16.6, 16.8) (Celestial a7f943a §C.13.B.5, arp_c/c13b_lemma_certify.out).
  Our Ω(12.875), Ω(16) values agree with Ω increasing through those brackets.

## 6. Effect on SECTOR_COERCIVITY_p2 and main G.4h (N.G2 pieces)

- Piece (4), odd archimedean margin: **CLOSED (CERTIFIED)** in the H-norm by Thm 1/3. The δ_o⁰, k_∞, B inputs of
  SECTOR_COERCIVITY_p2 Thm 2.1 (odd) are superseded; ℓ₀^odd = 0.028 explicitly. (δ_o⁰ and k_∞ were not certified; not needed.)
- Celestial Cor C.14.D, (K2) odd: closed with no coupling bound, now CERTIFIED.
- Even sector (K3): still needs m₀^e > 0 and C.14.K / condition (P). The same bathtub at ℓ = 0 gives ≈ 0 for even
  (−0.0014/+0.0006/+0.0016 at dt = 0.04/0.02/0.01, eval.): inconclusive, because the lowest even H-mode sits near
  t ≈ 4.53 = π/a (Celestial's first blind frequency), where Ω < 0, and the ĝ(0) direction is even. **GAP m₀^e** (item 3).

## 7. R-lines

- R1 Lemma A — PROVED. R2 Lemma B — PROVED.
- R3 m₀^odd(L₀) ≥ 0.46255 — CERTIFIED (Arb, 80/128 bits).
- R4 Conj C.13.H (odd, L₀; T = 12.875, μ = (0.6, 0.25); T = 16, μ = (0.7, 0.3)) — CERTIFIED, CC-free.
- R5 (K2) odd for ℓ ≤ 0.028, margin 0.0126‖ξ‖²_H — CERTIFIED given Prime S1 + (Ω5) (PROVED).
- R6 C.13.H′ range ℓ_H ≥ 0.0050 — PROVED given R4 (eval.).
- R7 even sector — NUMERICAL, inconclusive. GAP m₀^e.
- R8 t₀(B) — eval. given B; GAP B.
- R9 δ_o⁰, k_∞ not needed for odd — PROVED (logic).
- RH: OPEN.

## 8. Addendum (steering 1:20 AM): δ_o⁰ vs m₀^o, edge sign check, crossings

- **δ_o⁰ = 1 − sup spec(K_{L₀}|odd) is NOT certified here (GAP, same source as GAP B):** K_L has kernel κ = q/(2ε′(1⁺)),
  q(v) = (Qε)(e^{|v|}) (CC (102)/(104)), and certifying it needs CC's ε, not attempted. Golden Eye's numbers
  (sinc kernel on (−L₀, L₀), band |t| ≤ 6.29: λ₀ = 0.909 even, λ₁ = 0.424 odd, λ₂ = 0.053, so 1 − λ₁ ≈ 0.58) are for the
  plain low-band projector, a different operator from K_{L₀}; they make δ_o⁰ = O(1) plausible (CONJECTURE: δ_o⁰ = O(1)),
  but they do not bound it. **δ_o⁰ is off the critical path:** Thms 1–3 above close the odd sector without it.
- **m₀^o and the L² gap are different quantities.** m₀^o(L₀) = inf R/‖·‖²_H ≥ 0.46255 (CERTIFIED) is an H-norm margin of
  the archimedean form R alone. δ_o⁰ is an L²-margin of the ε-part of the trace-side split R = 𝒮 + 2ε′(1⁺)(1 − K_L).
  The CCM (arXiv:2310.18423) tiny eigenvalues live in the H-norm direction of the full Weil form (with prime terms) on larger
  windows. At the p = 2 entry they do not appear in the odd sector: Q_L = R there when ℓ = 0, and the certified H-margin is
  0.4626, resp. 0.0126 at ℓ = 0.028. The e^{−21/m₀} loss of Thm N.13.A is avoided by Thm 3 and C.13.H′ (polynomial ranges).
- **Edge sign (main 5354fde G.4h).** The ramp contribution to Q_L is +2c₂(1/ℓ − ℓ/24)α₋α₊, i.e. the matrix
  +c₂(1/ℓ − ℓ/24)[[0,1],[1,0]]; even is helped and odd is charged exactly 2c₂α²(1/ℓ − ℓ/24). **SECTOR_COERCIVITY_p2 (f1df9ed)
  matches (PROVED, checked):** its (K2) odd form R + 2c₂‖P₊ξ‖² ≥ 2c₂‖P₋ξ‖² charges the P₋ (odd-ramp) component. For a ramp,
  B = α₋α₊(ℓ/24 − 1/ℓ) = ‖P₋ξ‖²_{H(S₋)} − ‖P₊ξ‖²_{H(S₋)}, which is > 0 for the odd pair. Its sufficient condition
  R ≥ c₂‖ξ‖²_{H(E)} (S1) does not depend on the sign, and margin 2c₂s²ℓ = 2c₂α²/ℓ at α = sℓ is the leading order of
  2c₂α²(1/ℓ − ℓ/24). No fix needed. Thm 3 here also uses only S1, so it is sign-independent.
- **Crossings.** Ω = c₂ at 10.2615 and Ω = 2c₂ at 16.7481 are NUMERICAL (Golden Eye; re-evaluated here with mpmath findroot:
  10.2615356, 16.7480837). They lie in Celestial's Arb intervals (10.25, 10.30) and (16.6, 16.8).
- **Lemma A cross-check, re-run.** The closed form agrees with the exact Dirichlet-mode series (10⁶ modes, closed-form
  inner products) to 1e−15 relative at t = 1, 4.53, 10, 30 (eval.).
