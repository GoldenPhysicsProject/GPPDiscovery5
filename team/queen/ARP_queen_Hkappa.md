# (H-κ) — Conditioning of E on the prolate space: a polynomial bound, and what Λ^{o(1)} needs (Queen of Mathematics; appendix to ARP_queen_M1.md)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Queen of Mathematics (RH sub-team), 2026-10-28. Against GPPDiscovery5 main `25f8f83` (M1 merged at 17bf04e, fixes 9a83f78).

**RH is OPEN. Nothing here proves or assumes RH** except where a statement is explicitly marked "under RH" (and those are never
used for an RH conclusion). Analytic only: **no numerics were run**.
Labels: PROVED (citation) / PROVED (sketch; constants untracked — audit target) / CONJECTURE / GAP / HEURISTIC / OPINION.
R-lines as in §Q (R1 uniform in L; R2 Euler product + FE; R3 not via weight signs; R4 soft mode/pole; R5 fakes; R6 exact identity).

**Headline.**
- (H.1) **κ_Λ ≤ Λ^B unconditionally**, for an absolute (not explicit) B. PROVED (sketch; constants untracked). Ingredients:
  a two-constants lemma (Balasubramanian–Ramachandra's σ+2 device) and the Logvinenko–Sereda–Kovrijkine theorem. No zero-gap,
  simple-zero or pair-correlation input is used; it is fake-blind (it holds verbatim for Davenport–Heilbronn and Epstein), so the
  polynomial form of H-κ is **confirmed RH-silent**.
- (H.2) **κ_Λ ≤ Λ^{o(1)}** reduces to one named input (Z-avg): uniform window averages of log|ζ(½+is)| over windows of length
  a → 0 are ≥ −o(log T), together with an o(log T) bound on window averages of log⁺|ζ|. The second half is implied by Lindelöf; the
  first is CONJECTURE (expected under RH; not known unconditionally). The obstruction is the *background* size of ζ on a whole window
  (large deviations of log|ζ|), **not** isolated zeros, multiple zeros, or small gaps.
- (H.3) Survivability. With κ ≤ Λ^{o(1)}, M1 needs only ε_Λ/Λ + ρ_Λ ≤ Λ^{−δ}. With the unconditional κ ≤ Λ^B it needs
  ε_Λ/Λ + ρ_Λ ≤ Λ^{−B−δ}. That is available (CONJECTURE, for Spectrum) only for a **trimmed** prolate space; for Connes' untrimmed B_Λ
  (which keeps half the plunge) it is unknown, because B is not explicit.

## H.0 Precise statement

Notation from M1.4: B_{Λ,0} ⊂ L²(ℝ)^ev is the span of the retained even prolates ψ_n (n ∈ I_Λ), truncated to [−Λ, Λ], with
φ(0) = φ̂(0) = 0. Mφ(s) := ∫_0^∞ φ(x) x^{½+is} d^*x (s ∈ ℝ). Then Mφ(s) = ∫ F_φ(u) e^{isu} du with F_φ(u) := e^{u/2}φ(e^u), and
∫|Mφ|² ds = 2π ∫_0^∞|φ|² dx = π‖φ‖²_{L²(ℝ)} (Plancherel; PROVED).

  **m_Λ := min_{0≠φ∈B_{Λ,0}} ∫ |ζ(½+is)|² |Mφ(s)|² ds / ∫ |Mφ(s)|² ds,  κ_Λ = 1/(2 m_Λ)** (M1.4 Mellin–Plancherel, PROVED on S_0^ev).

- **(H-κ, polynomial form):** m_Λ ≥ Λ^{−B}. **(H-κ, sharp form):** m_Λ ≥ Λ^{−o(1)}.
- m_Λ > 0 for each Λ (M1.4; finite dimension). The issue is only the rate.
- Since m_Λ is a minimum over a subspace, trimming I_Λ can only increase m_Λ (PROVED, trivial).

## H.1 Structural facts about the profiles |Mφ|²

Write Φ for the entire extension of φ ∈ B_{Λ,0} (Φ = Σ c_n ψ_n^{ent}, where ψ_n^{ent} is the whole-line prolate. Φ̂ is supported in
[−Λ, Λ] and ‖Φ‖²_{L²(ℝ)} = Σ |c_n|²/λ_n ≤ ‖φ‖²/λ_min, with λ_min := min_{I_Λ} λ_n; double orthogonality, Slepian–Pollak 1961).
λ_min ≥ ½ − o(1) for Connes' retention (n ≲ 2c/π, c = 2πΛ²), and λ_min ≥ 1 − Λ^{−K} for a trimmed retention.

**Fact H.1a (u-support; PROVED).** F_φ(u) = 0 for u > log Λ. For u → −∞, Φ even with Φ(0) = 0 gives
|Φ(x)| ≤ ½x²‖Φ″‖_∞ ≤ ½x²(2πΛ)²(2Λ)^{1/2}‖Φ‖ (Bernstein, plus ‖Φ‖_∞ ≤ (2Λ)^{1/2}‖Φ‖ for band-limited Φ). Hence
∫_{u<−log Λ−A} |F_φ|² du ≤ C e^{−5A} ‖φ‖²/λ_min.
So F_φ is, up to relative L²-error O(e^{−5A/2}λ_min^{−1/2}), supported in an interval of length **ℓ = 2 log Λ + A**. Consequently
Mφ is, to that accuracy, a Paley–Wiener function of type ~ log Λ. It is **not** exactly Paley–Wiener: the left tail decays like
e^{5u/2}, and Mφ is meromorphic with poles at s = i(½ + 2k), k ≥ 1.

**Fact H.1b (resolution; PROVED in Kovrijkine form, HEURISTIC in the phase-space form).** By H.1a and H.2 below, no profile can put
more than a fixed fraction of its mass on a set that misses a fixed fraction of every window of length ≍ 1/log Λ. Phase-space
picture (HEURISTIC): the time–band region {|x| ≤ Λ, |ξ| ≤ Λ} maps to {(u, s) : s ≲ 2πΛe^u, u ≤ log Λ}, of area 4πΛ²; the dimension
is 2Λ², as it should be. So at height s the available u-extent is only ≈ log(2πΛ²/s). Profiles at heights ≍ Λ² therefore have
**O(1)** resolution, which averages over ≍ log Λ zero spacings (spread case). Resolution ≈ 1/log Λ is reached only at heights Λ^{2−α}
with α bounded below (hard case).

**Fact H.1c (height range and escape to large |s|; PROVED).** Integrating by parts in Mφ:
Mφ(s) = [φ(Λ)Λ^{½+is} − M(xφ′·1_{[0,Λ]})(s)]/(½ + is). Using |φ(Λ)|² ≤ 2Λ‖φ‖²/λ_min and ∫_0^Λ x²|φ′|² ≤ ½Λ²(2πΛ)²‖φ‖²/λ_min
(Bernstein on Φ):

  ∫_{|s|>T} |Mφ|² ds / ∫ |Mφ|² ds ≤ C (Λ²/T + Λ⁴/T²)/λ_min  for every T > 0.

So the lead's height ~Λ² is confirmed as the mass scale. Mass **can** escape above Λ², but only through the 1/|s| jump tail at the
endpoint x = Λ, and only a fraction O(Λ²/T) beyond T. Below |s| ≤ 20, |ζ(½+is)| ≥ c₀ > 0 (no zeros below 14.13), so low heights are
harmless.

R-line H.1: R1 Λ-uniform, no L | R2 not engaged | R3 n/a | R4 the pole conditions φ(0) = φ̂(0) = 0 are what make F_φ decay at −∞ |
R5 archimedean only | R6 n/a.

## H.2 Prior art (read via secondary sources this session; cite, don't duplicate)

| Source | What it gives | Does it bound the MIN m_Λ? |
|---|---|---|
| Balasubramanian–Conrey–Heath-Brown 1985 (twisted/mollified 2nd moment ∫\|ζA\|² for Dirichlet polynomials A of length T^θ, θ < ½; extended to θ < 17/33 by Bettin–Chandee–Radziwiłł 2017) | asymptotics for one fixed weight \|A\|² | No. It is a fixed weight, and our profiles are not Dirichlet polynomials. |
| Montgomery–Vaughan mean value theorem ∫_0^T \|Σa_n n^{−it}\|² = Σ(T + O(n))\|a_n\|² | upper and lower bounds for Dirichlet polynomials over long ranges | No. It covers long intervals only, and a minimizer can sit in a resolution cell. |
| Short-interval 2nd moment: ∫_t^{t+H}\|ζ\|² ~ H log t for H ≥ T^{θ}, θ ≈ 1/3 (Ingham–Balasubramanian–Huxley error term; Ivić's book, mean-value chapters) | window averages at long scales | No. The range H ≥ T^{1/3} is far above our scale 1/log T. |
| Balasubramanian–Ramachandra lower bounds: ∫_{T−H}^{T+H}\|ζ(½+it)\|^{2k} ≫ H(log H)^{k²} for C_k log log T ≤ H ≤ T, k ∈ ℕ, unconditional (Ramachandra's monograph; Ivić, Hardy-Ramanujan J. 2009, Thms 1–2; under RH for real k) | short-window lower bounds down to H ≍ log log T | Only for weights flat at scale ≥ C log log T (Prop H.4). The *method* (compare with σ + 2, where ζ ≫ 1) is what Lemma H.2 uses. |
| Ramachandra Ω-results; Balasubramanian–Ramachandra max_{[T,T+H]}\|ζ\| ≥ exp(c(log H/log log H)^{1/2}) for H ≥ 100 log log T | large values | No (wrong direction). |
| Heath-Brown 1981 fractional moments; Selberg CLT (log\|ζ\| ≈ N(0, ½ log log T)); Radziwiłł–Soundararajan (simple proof of Selberg's CLT; one-sided lower-bound method) | distribution and lower tails *on average over t* | No uniform-in-window statement. They motivate (Z-avg) (HEURISTIC). |
| Logvinenko–Sereda 1974; **Kovrijkine 2001** (Proc. AMS 129, 3037–3047, Thm 1): if E is (γ, a)-thick and supp ĝ ⊂ an interval of length b, then ‖g‖_{L^p(E)} ≥ (γ/C)^{C(ab+1)}‖g‖_p | uniform over all band-limited g | **Yes.** This is the tool that turns a pointwise-in-window statement about ζ into a bound on the minimum. |

Why fixed-weight results fail. m_Λ = min over a 2Λ²-dimensional family of weights, and the minimizer adapts to ζ. Any bound of the
form ∫|ζ|²w ≈ (log)∫w for one w (or for an average over w) says nothing about the worst w. A min bound needs a statement that holds
**in every window** at the resolution scale, which is exactly what thickness gives. (This refines the M1 audit note: the min does need
local information on the line, but **thickness of a superlevel set**, not zero spacing.)

## H.3 Theorem H.1 (polynomial conditioning; PROVED (sketch; constants untracked — audit target))

**Lemma H.2 (two-constants; PROVED).** There are absolute c₀, C₀ > 0 with the following property. Let T ≥ 30, |t₀| ∈ [20, T],
0 < a ≤ 1, and let F ⊂ [t₀ − a/2, t₀ + a/2] be measurable with |F| ≥ a/2. Then
sup_{s∈F} |ζ(½+is)| ≥ η_a(T) := T^{−C₀/a}.
Equivalently, E_a(T) := {s : |ζ(½+is)| ≥ η_a(T)} ∪ {|s| ≤ 20} ∪ {|s| > T} is (½, a)-thick.

*Proof.* Let D be the rectangle (½, 5/2) × (t₀ − 2, t₀ + 2) and z₀ = 2 + it₀. log|ζ| is subharmonic on a neighbourhood of D̄ (the pole
at 1 lies outside), so log|ζ(z₀)| ≤ ∫_{∂D} log|ζ| dω_{z₀}.
- On ∂D, |ζ| ≤ C T^{1/4} log T (convexity bound for σ ≥ ½).
- On F ⊂ left edge, if |ζ| ≤ η there, then ω_{z₀}(F) ≥ c₁|F| ≥ c₁a/2, since the harmonic-measure density on the central part of the
  left edge is bounded below by an absolute c₁ (fixed shape).
- |ζ(2+it₀)| ≥ ζ(4)/ζ(2) > 0.65 (Euler product).
- So log 0.65 ≤ (c₁a/2) log η + (1/4 + o(1)) log T, i.e. log(1/η) ≤ C₀ log T/a.
This is Balasubramanian–Ramachandra's σ + 2 device made pointwise on a set. ∎

**Theorem H.1.** m_Λ ≥ Λ^{−B} for all large Λ, with B absolute. (For Connes' retention λ_min ≥ ½ − o(1); for any retention, the
loss is at most a factor λ_min^{−O(1)}.)

*Proof sketch.*
1. Fix small a > 0 and A = A₁ log Λ (A₁ chosen below). By H.1a, G := Mφ = G_PW + r, where G_PW is the transform of F_φ·1_{[−log Λ−A, log Λ]},
   band-limited to an interval of length ℓ = (2 + A₁) log Λ, and ‖r‖ ≤ CΛ^{−5A₁/2}‖G‖ (λ_min ≍ 1).
2. Kovrijkine with E = E_a(T), γ = ½, b = ℓ:
   ‖G_PW‖_{L²(E)} ≥ K^{−1}‖G_PW‖, K := (2C)^{C(aℓ+1)} = Λ^{C a(2+A₁) log(2C)} · O(1).
3. Choose a so that C a(2 + A₁) log(2C) < 5A₁/4. Then ‖r‖ ≪ K^{−1}‖G‖, so ‖G‖_{L²(E)} ≥ (2K)^{−1}‖G‖.
4. Remove {|s| > T} using H.1c with T := C′Λ²K². The removed mass is ≤ (4K)^{−2}‖G‖², which leaves
   ∫_{E∩[−T,T]}|G|² ≥ (4K)^{−2}‖G‖².
5. On E ∩ [−T, T], |ζ|² ≥ min(η_a(T)², c₀²). So m_Λ ≥ η_a(T)²(4K)^{−2}. Here T = Λ^{O(1)} and η_a(T) = T^{−C₀/a} = Λ^{−O(1/a)}.
   ∎ (B ≍ 1/a + a; any fixed admissible a gives an absolute B.)

R-line H.1: R1 Λ-uniform, no L | R2 the Euler product enters only as |ζ(2+it)| ≫ 1, the FE only via the convexity bound;
not jointly in any essential way | R3 n/a | R4 the pole conditions give H.1a | R5 **fake-blind** (H.6) | R6 n/a.

**Prop H.4 (spread case via Balasubramanian–Ramachandra; PROVED (citation)).** If w ≥ 0 is θ-flat at scale H ≥ C log log T on
[20, T] (w(s) ≥ θ w(s′) whenever |s − s′| ≤ H), then ∫_{20}^T |ζ|² w ≥ c θ log H ∫_{20}^T w. *Proof:* tile by windows of length 2H and
apply B–R (k = 1) on each. ∎ This gives the "spread case" at the correct size log H, without loss. Our profiles are, however, **not**
θ-flat in general: they oscillate and vanish, and the top-height profiles are spread only at scale O(1) ≪ log log T (H.1b). So
Prop H.4 does not by itself bound m_Λ; that is why H.1 uses thickness instead. No unconditional window lower bound for ∫|ζ|² below
H ≍ log log T is known to me (GAP in the literature, as far as this survey goes).

## H.4 The hard case and the exact input for Λ^{o(1)}

The loss in H.1 comes from two places: η_a = T^{−C₀/a} (two-constants is lossy) and K = Λ^{O(a)} (needs a → 0). Kovrijkine is
sharp up to C, so **any** improvement must come from ζ.

**Input (Z-avg) (the named input; GAP).** There is a(T) → 0 such that, uniformly over windows I ⊂ [20, T] of length a(T),
- (Z-avg⁻) (1/|I|) ∫_I log|ζ(½+is)| ds ≥ −o(log T), and
- (Z-avg⁺) (1/|I|) ∫_I log⁺|ζ(½+is)| ds ≤ o(log T).

**Prop H.5 (PROVED, given Kovrijkine).** (Z-avg) ⇒ m_Λ ≥ Λ^{−o(1)}.
*Proof:* by Chebyshev, the sublevel set {log|ζ| < −2(V⁺ + V⁻)} has measure ≤ |I|/2 in each window, so E is (½, a(T))-thick at
level η = T^{−o(1)}. Then run Theorem H.1 with a = a(T) → 0, which gives K = Λ^{o(1)}. ∎

Status of the two halves:
- (Z-avg⁺) follows from Lindelöf by Jensen, (1/|I|)∫log⁺|ζ| ≤ ½ log(1 + |I|^{−1}∫_I|ζ|²). RH ⇒ LH. Under RH it also holds pointwise:
  log|ζ| ≤ (log 2/2 + o(1)) log T/log log T (Chandee–Soundararajan; Carneiro–Chandee–Milinovich). Unconditionally only O(log T)
  (Weyl–Hardy–Littlewood t^{1/6}).
- (Z-avg⁻) is **CONJECTURE**. HEURISTIC support: the Selberg CLT lower tail plus a union bound over T/a windows gives a threshold of
  ≍ (log T log log T)^{1/2} = o(log T). Under RH it should follow from Guinand–Weil with Beurling–Selberg minorants of the windowed
  logarithmic kernel (the Carneiro–Chandee–Milinovich method). This is **not verified**: CONJECTURE (under RH). Unconditionally,
  only −O(log T/a) (Lemma H.2).

**What zeros do and do not matter (PROVED, elementary).** Near t₀, log|ζ(½+is)| = Σ_{|γ−t₀|<1} log|s − ρ| + h(s) with h = O(log T)
(Titchmarsh §9.6(A)). The window average of log|s − γ| over a window of length a is ≥ log(a/2e) for *every* γ, real or not. So a
cluster of m zeros in or near the window costs at most m·log(2e/a) in (Z-avg⁻).
- Simple zeros, multiple zeros, small gaps and pair correlation are **not** the obstruction at either level. A double zero or a gap
  of size T^{−100} costs O(log(1/a)).
- Cluster size m ≤ N(t₀+1) − N(t₀−1) ≪ log T (Backlund) gives polynomial loss only, the same as H.1. Under RH,
  S(t) ≪ log T/log log T makes m log(1/a) = o(log T) for a log(1/a) → 0.
- The real obstruction is the background h(s), the Dirichlet-polynomial part, being uniformly ≲ −ε log T on a whole window: a
  uniform large-deviation statement for log|ζ|.
- **Off-line zeros do not create small values on the line**: a zero at distance β − ½ > 0 contributes log|s − ρ| ≥ log(β − ½) > −∞.
  RH enters (Z-avg) only as a *tool* (S(t), Lindelöf), not structurally.
- Labels: (Z-avg) GAP. "Implied by RH": (Z-avg⁺) PROVED; (Z-avg⁻) CONJECTURE (expected). "Implies RH": no (OPINION, but see H.6:
  any route to (Z-avg) that works for DH cannot imply RH).

## H.5 Payoff for M1 (PROVED-conditional)

From M1.C (ℓ_Λ ≤ 2κ_Λ(ε_Λ/Λ + 4ζ(1+σ)²ρ_Λ)) and M1.5 (the diagonal Λ = e^L needs ℓ_Λ = o(1/log²Λ)):
- **Sharp form.** If m_Λ ≥ Λ^{−o(1)} (from (Z-avg)) and Spectrum delivers ε_Λ/Λ + ρ_Λ ≤ Λ^{−δ} for some δ > 0, then
  ℓ_Λ ≤ Λ^{−δ+o(1)} = o(1/log²Λ). So M1 holds, including uniformly on Q.8.E's diagonal. PROVED-conditional.
- **Unconditional form.** Theorem H.1 gives κ_Λ ≤ ½Λ^B. Then M1 (diagonal version included) holds if ε_Λ/Λ + ρ_Λ ≤ Λ^{−B−δ}.
  PROVED-conditional.
- **Is polynomial loss survivable? Honest answer: only with a trimmed prolate space, as far as is known.**
  - With Connes' retention (n up to the middle of the plunge, λ_n ≈ ½ at the edge), the plunge states have endpoint values
    ≍ Λ^{−1/2}. HEURISTIC: ε_Λ/Λ ≍ log Λ/Λ², polynomial with exponent ≈ 2. Since B is not explicit (Kovrijkine's C, two-constants'
    C₀), B < 2 is not available. So: **unknown**.
  - With trimmed retention I_Λ^{(K)} := {n : 1 − λ_n ≤ Λ^{−K}}: κ only improves (H.0). Plateau endpoint values and remainder tails
    should then be O(Λ^{−K′}) with K′ → ∞ as K → ∞ (CONJECTURE, for Spectrum: via ∂_cλ_n = (2/c)λ_nψ_n(1)² and plateau bounds on
    1 − λ_n). That makes any polynomial κ survivable.
  - Cost of trimming. It removes O(log Λ · K) plunge states from Q′_Λ. Their Mellin mass sits at heights ≍ Λ² (H.1b), where ĥ(s) is
    negligible for smooth h. So the change in Tr(Q′V(h)) should be o(1) for fixed f (CONJECTURE; it needs a plunge-state Mellin
    localization bound, which is Spectrum's type of input). It must then be checked against M2's counting term.

## H.6 Falsifiers (R5)

- **Davenport–Heilbronn** f_DH (a 5-periodic Dirichlet series, a₁ = 1, |a_n| ≤ 1; satisfies the FE, no Euler product, has off-line
  zeros). Replace ζ by f_DH in m_Λ (this is the conditioning of E_DH(φ)(x) = x^{1/2}Σ a_n φ(nx)). Lemma H.2 holds verbatim with
  z₀ = 3 + it₀ and the rectangle (½, 7/2) × (t₀ − 3, t₀ + 3):
  - |f_DH(3+it)| ≥ 1 − (ζ(3) − 1) > 0.79;
  - f_DH is entire with polynomial growth on the line (FE + Phragmén–Lindelöf).
  Kovrijkine is pure Fourier analysis. So **Theorem H.1 holds for DH: m_Λ^{DH} ≥ Λ^{−B′}**. PROVED (sketch).
- **Epstein ζ_Q** (class number > 1: FE, no Euler product, off-line zeros): same argument, anchoring at σ₀ large enough that the
  leading term of the Dirichlet series dominates. PROVED (sketch).
- **F_θ / fakes with an Euler product but no FE** (§Q.8.5): H.2 needs only polynomial growth on the line and a nonvanishing anchor,
  so it typically also passes. That is consistent: M1 is not supposed to see RH.
- **Conclusion.** The polynomial form of H-κ is **fake-blind**, and therefore **RH-silent, confirmed** (it cannot carry RH content).
  For the sharp form, the RH-based route to (Z-avg⁻) is not available for DH (DH violates RH, and Lindelöf for DH is open), so its
  fake-blindness is **not established**. OPINION: (Z-avg) is a statement about the background size of the Dirichlet polynomial, and
  should hold for DH as well; if a proof of (Z-avg) for ζ used RH essentially, M1's sharp form would become RH-conditional. That
  would be harmless for M1 as a tool, but it must then **not** be used on any route to RH (circularity flag).

## H.7 Summary table

| Claim | Label |
|---|---|
| H.0 m_Λ definition, κ_Λ = 1/(2m_Λ), monotone under trimming | PROVED |
| H.1a F_φ effectively supported in length 2log Λ + A, error e^{−5A/2} | PROVED |
| H.1b resolution ≈ 1/log Λ at low heights, O(1) at heights ≍ Λ² | Kovrijkine form PROVED; phase-space form HEURISTIC |
| H.1c mass fraction beyond T is O(Λ²/T) | PROVED |
| Lemma H.2 two-constants thickness at level T^{−C₀/a} | PROVED |
| Theorem H.1 m_Λ ≥ Λ^{−B}, B absolute | PROVED (sketch; constants untracked — audit target) |
| Prop H.4 spread weights (flat at scale ≥ C log log T) | PROVED (citation: Balasubramanian–Ramachandra) |
| Prop H.5 (Z-avg) ⇒ m_Λ ≥ Λ^{−o(1)} | PROVED-conditional |
| (Z-avg⁺) | PROVED under Lindelöf (so under RH); GAP unconditionally |
| (Z-avg⁻) | CONJECTURE (expected under RH; HEURISTIC via Selberg CLT); GAP |
| zeros/gaps/multiplicity not the obstruction | PROVED (elementary) |
| H.5 sharp payoff: ℓ_Λ = o(log^{−2}Λ) given (Z-avg) and Spectrum's Λ^{−δ} | PROVED-conditional |
| H.5 polynomial loss survivable with trimmed retention | CONJECTURE (needs Spectrum's trimmed bounds and the trimming-harmless check) |
| H.6 polynomial H-κ holds for DH and Epstein ⇒ RH-silent | PROVED (sketch) |

Not cited: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.
