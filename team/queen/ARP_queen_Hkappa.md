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
- **Circularity flag (Grok audit, plain statement).** The sharp form κ_Λ ≤ Λ^{o(1)} may be RH-conditional: its only known supports are Lindelöf, which follows from RH, for (Z-avg⁺), and an expected RH-based argument for (Z-avg⁻). **It cannot be used on any route to RH.** Only the polynomial form (Theorem H.1) is shown fake-blind (DH, Epstein, H.6), and only it is RH-silent.
- (H.3) Survivability. With κ ≤ Λ^{o(1)}, M1 needs only ε_Λ/Λ + ρ_Λ ≤ Λ^{−δ}. With the unconditional κ ≤ Λ^B it needs
  ε_Λ/Λ + ρ_Λ ≤ Λ^{−B−δ}. That is available (CONJECTURE, for Spectrum) only for a **trimmed** prolate space; for Connes' untrimmed B_Λ
  (which keeps half the plunge) it is unknown, because B is not explicit.
- **Part II (below, H.8–H.11).** Window-average Lemma H.2′ (PROVED (citation)) ⇒ Theorem H.1′ κ ≤ Λ^{B′} (PROVED (sketch)). Λ^{o(1)} is not reachable by this route. Trimmed space: M1^{(K)} holds for K ≥ 2B + 8 (Theorem H.9, PROVED-conditional on the H.1 sketch + M1.1 GAP-minor). The only open interface item is (T3), the trimming-harmless check.

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
endpoint x = Λ, and only a fraction O(Λ²/T) beyond T. ~~Below |s| ≤ 20, |ζ(½+is)| ≥ c₀ > 0~~ *(Grok audit fix: false as written, since γ₁ = 14.1347… < 20.)* Corrected: for |s| ≤ 14,
|ζ(½+is)| ≥ c₀ > 0 (no zeros with |γ| < 14.13). Heights 14 ≤ |s| ≤ 20 are handled by Lemma H.2, which needs only |t₀| ≥ 3, so low
heights are harmless.

R-line H.1: R1 Λ-uniform, no L | R2 not engaged | R3 n/a | R4 the pole conditions φ(0) = φ̂(0) = 0 are what make F_φ decay at −∞ |
R5 archimedean only | R6 n/a.

## H.2 Prior art (read via secondary sources this session; cite, don't duplicate)

| Source | What it gives | Does it bound the MIN m_Λ? |
|---|---|---|
| Balasubramanian–Conrey–Heath-Brown 1985 (twisted/mollified 2nd moment ∫\|ζA\|² for Dirichlet polynomials A of length T^θ, θ < ½; extended to θ < 17/33 by Bettin–Chandee–Radziwiłł 2017) | asymptotics for one fixed weight \|A\|² | No. It is a fixed weight, and our profiles are not Dirichlet polynomials. |
| Montgomery–Vaughan mean value theorem ∫_0^T \|Σa_n n^{−it}\|² = Σ(T + O(n))\|a_n\|² | upper and lower bounds for Dirichlet polynomials over long ranges | No. It covers long intervals only, and a minimizer can sit in a resolution cell. |
| Short-interval 2nd moment: ∫_t^{t+H}\|ζ\|² ~ H log t for H ≥ T^{θ}, θ ≈ 1/3 (Ingham–Balasubramanian–Huxley error term; Ivić's book, mean-value chapters) | window averages at long scales | No. The range H ≥ T^{1/3} is far above our scale 1/log T. |
| Balasubramanian–Ramachandra lower bounds: ∫_{T−H}^{T+H}\|ζ(½+it)\|^{2k} ≫ H(log H)^{k²} for C_k log log T ≤ H ≤ T, k ∈ ℕ, unconditional (Ramachandra's monograph; Ivić, Hardy-Ramanujan J. 32 (2009) 4–23, Thm 2 = Balasubramanian–Ramachandra, k ∈ ℕ, C(ε,k) log log T ≤ H ≤ T — checked this session (Grok audit); under RH for real k) | short-window lower bounds down to H ≍ log log T | Only for weights flat at scale ≥ C log log T (Prop H.4). The *method* (compare with σ + 2, where ζ ≫ 1) is what Lemma H.2 uses. |
| Ramachandra Ω-results; Balasubramanian–Ramachandra max_{[T,T+H]}\|ζ\| ≥ exp(c(log H/log log H)^{1/2}) for H ≥ 100 log log T | large values | No (wrong direction). |
| Heath-Brown 1981 fractional moments; Selberg CLT (log\|ζ\| ≈ N(0, ½ log log T)); Radziwiłł–Soundararajan (simple proof of Selberg's CLT; one-sided lower-bound method) | distribution and lower tails *on average over t* | No uniform-in-window statement. They motivate (Z-avg) (HEURISTIC). |
| Logvinenko–Sereda 1974; **Kovrijkine 2001** (Proc. AMS 129, 3037–3047, Thm 1): if E is (γ, a)-thick and supp ĝ ⊂ an interval of length b, then ‖g‖_{L^p(E)} ≥ (γ/C)^{C(ab+1)}‖g‖_p | uniform over all band-limited g | **Yes.** This is the tool that turns a pointwise-in-window statement about ζ into a bound on the minimum. |

Why fixed-weight results fail. m_Λ = min over a 2Λ²-dimensional family of weights, and the minimizer adapts to ζ. Any bound of the
form ∫|ζ|²w ≈ (log)∫w for one w (or for an average over w) says nothing about the worst w. A min bound needs a statement that holds
**in every window** at the resolution scale, which is exactly what thickness gives. (This refines the M1 audit note: the min does need
local information on the line, but **thickness of a superlevel set**, not zero spacing.)

## H.3 Theorem H.1 (polynomial conditioning; PROVED (sketch; constants untracked — audit target))

**Lemma H.2 (two-constants; PROVED).** There are absolute c₀, C₀ > 0 with the following property. Let T ≥ 30, |t₀| ∈ [3, T] *(Grok audit: was [20, T]; any |t₀| ≥ 3 keeps the rectangle off the pole at s = 1)*,
0 < a ≤ 1, and let F ⊂ [t₀ − a/2, t₀ + a/2] be measurable with |F| ≥ a/2. Then
sup_{s∈F} |ζ(½+is)| ≥ η_a(T) := T^{−C₀/a}.
Equivalently, E_a(T) := {s : |ζ(½+is)| ≥ η_a(T)} ∪ {|s| ≤ 4} ∪ {|s| > T} is (½, a)-thick.
*(Grok audit: the set was {|s| ≤ 20}, which contains the zero at 14.13 and breaks step 5 of H.1. With {|s| ≤ 4}, |ζ| ≥ c₀ there (no zeros), and every window of length a ≤ 1 not inside {|s| ≤ 4} ∪ {|s| > T} has its centre in 3 ≤ |t₀| ≤ T + 1, where the lemma applies (take T + 3 in place of T).*
*Thickness is the contrapositive with F := the sublevel set {|ζ| < η} ∩ window: if it had measure ≥ a/2, the lemma would be contradicted. ζ(½−is) is the conjugate of ζ(½+is), so negative s are covered.)*

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
*(Grok audit: the parenthetical holds only if λ_min ≥ Λ^{−O(1)}. λ_min enters step 1's error as λ_min^{−1/2}, and absorbing it needs a larger A₁, so it changes B, not just a prefactor. For Connes' and trimmed retention λ_min ≍ 1, so B is absolute as stated.)*

*Proof sketch.*
1. Fix small a > 0 and A = A₁ log Λ (A₁ chosen below). By H.1a, G := Mφ = G_PW + r, where G_PW is the transform of F_φ·1_{[−log Λ−A, log Λ]},
   band-limited to an interval of length ℓ = (2 + A₁) log Λ, and ‖r‖ ≤ CΛ^{−5A₁/2}‖G‖ (λ_min ≍ 1).
2. Kovrijkine with E = E_a(T), γ = ½, b = ℓ:
   ‖G_PW‖_{L²(E)} ≥ K^{−1}‖G_PW‖, K := (2C)^{C(aℓ+1)} = Λ^{C a(2+A₁) log(2C)} · O(1).
3. Choose a so that C a(2 + A₁) log(2C) < 5A₁/4. Then ‖r‖ ≪ K^{−1}‖G‖, so ‖G‖_{L²(E)} ≥ (2K)^{−1}‖G‖.
4. Remove {|s| > T} using H.1c with T := C′Λ²K². The removed mass is ≤ (4K)^{−2}‖G‖², which leaves
   ∫_{E∩[−T,T]}|G|² ≥ (4K)^{−2}‖G‖².
5. On E ∩ [−T, T], |ζ|² ≥ min(η_a(T)², c₀²) *(c₀ := min_{|s|≤4}|ζ(½+is)| > 0; Grok audit fix)*. So m_Λ ≥ η_a(T)²(4K)^{−2}. Here T = Λ^{O(1)} and η_a(T) = T^{−C₀/a} = Λ^{−O(1/a)}.
   ∎ (B ≍ 1/a + a; any fixed admissible a gives an absolute B.)

*Grok audit (S103e): steps checked.*
- **Lemma H.2.** Sub-mean value for subharmonic log|ζ| on D̄. The pole at s = 1 is outside D for |t₀| ≥ 3.
  - The harmonic-measure density from z₀ = 2 + it₀ on the central part (|t − t₀| ≤ ½) of the left edge is ≥ c₁. This holds because the shape is fixed and the edge is smooth away from the corners.
  - Convexity gives |ζ| ≪ T^{1/4} log T on ∂D for σ ≥ ½. |ζ(2+it)| ≥ ζ(4)/ζ(2) = 0.658.
  - Hence log(1/η) ≤ (¼ log T + O(1))/(c₁a/2). ✓
- **H.1a.** Φ(0) = Φ′(0) = 0, Bernstein ‖Φ″‖_∞ ≤ (2πΛ)²‖Φ‖_∞, and ‖Φ‖_∞ ≤ (2Λ)^{1/2}‖Φ‖. These give ∫_{u<−logΛ−A}|F|² ≤ Ce^{−5A}‖Φ‖². ✓
- **Kovrijkine.** Kovrijkine, Proc. AMS 129 (2001) no. 10, 3037–3047, Thm 1: ‖g‖_{L^p(E)} ≥ (γ/C)^{C(ab+1)}‖g‖_p for supp ĝ in an interval of length b and E (γ, a)-thick. Bibliographic data and statement checked this session; C is absolute.
  - It is applied to G_PW, which is exactly band-limited (u-support of length ℓ), and not to Mφ itself. ✓
- **Error e^{−5A/2}.** ‖G‖_E ≥ K^{−1}‖G_PW‖ − ‖r‖ ≥ K^{−1}(1 − o(1))‖G‖ − CΛ^{−5A₁/2}‖G‖.
  - Step 3's constraint makes K^{−1} ≥ Λ^{−5A₁/4} ≫ Λ^{−5A₁/2}. An admissible a exists for every A₁ > 0. ✓
- **Step 4.** T = C′Λ²K² is fixed before E_a(T), and K does not depend on T, so there is no circularity. H.1c's integration by parts and Plancherel check. ✓
- **Hole found and fixed: the zero at 14.13** (the low-height set). The sketch is otherwise complete up to constants, so the label is **PROVED (sketch)**.

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
  log|ζ| ≤ (log 2/2 + o(1)) log T/log log T (Chandee–Soundararajan 2009, constant log 2/2; Carneiro–Chandee–Milinovich cited for the method family — *attribution not verified this session (Grok audit)*). Unconditionally only O(log T)
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
    *(Part II, H.9: now settled without Spectrum. The endpoint bound is Lemma H.9.1 (elementary, exponent 2/3), and Theorem H.9 gives ℓ^{(K)} ≤ CΛ^{−1} for K ≥ 2B + 8.)*
  - Cost of trimming. It removes O(K log²Λ) plunge states from Q′_Λ *(Queen correction, Part II (T1): was "O(log Λ · K)"; with ε = Λ^{−K} the KRD count is β_Λ log(1/ε) ≍ K log²Λ)*. Their Mellin mass sits at heights ≍ Λ² (H.1b), where ĥ(s) is
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
- **Circularity flag (Grok audit; not OPINION).** As things stand, the sharp Λ^{o(1)} form **may be RH-conditional** (its only supports are Lindelöf ⇐ RH and an expected RH argument), so it **cannot be used on any route to RH**. Only the polynomial form is shown fake-blind (DH, Epstein).

## H.7 Summary table

| Claim | Label |
|---|---|
| H.0 m_Λ definition, κ_Λ = 1/(2m_Λ), monotone under trimming | PROVED |
| H.1a F_φ effectively supported in length 2log Λ + A, error e^{−5A/2} | PROVED |
| H.1b resolution ≈ 1/log Λ at low heights, O(1) at heights ≍ Λ² | Kovrijkine form PROVED; phase-space form HEURISTIC |
| H.1c mass fraction beyond T is O(Λ²/T) | PROVED |
| Lemma H.2 two-constants thickness at level T^{−C₀/a} | PROVED |
| Theorem H.1 m_Λ ≥ Λ^{−B}, B absolute | PROVED (sketch; constants untracked; Grok audit S103e: steps checked, low-height zero fix applied) |
| Prop H.4 spread weights (flat at scale ≥ C log log T) | PROVED (citation: Balasubramanian–Ramachandra) |
| Prop H.5 (Z-avg) ⇒ m_Λ ≥ Λ^{−o(1)} | PROVED-conditional |
| (Z-avg⁺) | PROVED under Lindelöf (so under RH); GAP unconditionally |
| (Z-avg⁻) | CONJECTURE (expected under RH; HEURISTIC via Selberg CLT); GAP |
| zeros/gaps/multiplicity not the obstruction | PROVED (elementary) |
| H.5 sharp payoff: ℓ_Λ = o(log^{−2}Λ) given (Z-avg) and Spectrum's Λ^{−δ} | PROVED-conditional; may be RH-conditional through (Z-avg), so barred from any route to RH (circularity flag) |
| H.5 polynomial loss survivable with trimmed retention | superseded by H.9: M1 for the trimmed space PROVED-conditional (H.1 sketch + M1.1 GAP-minor); trimming-harmless (T3) open, split (T3a) HEURISTIC/Spectrum + (T3b) GAP/Queen |
| H.6 polynomial H-κ holds for DH and Epstein ⇒ RH-silent | PROVED (sketch) |

Not cited: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.

---

# Part II (2026-10-28, after the H-κ audit): window-average conditioning, and the trimmed-space interface

Against GPPDiscovery5 main `0000f61` (contains 8ecee4f; H-κ merged at 47a8109, main a864c92). Analytic only: **no numerics were run.**
RH is OPEN. Nothing in Part II assumes RH, except the one flagged line in H.10, which is barred from every route.

**Part II headline.**
- (H.8) **Lemma H.2′ (window average; PROVED (citation)).** On every window of length a ≤ 1 at height ≤ T, the average of log|ζ(½+is)| is ≥ −C log T · log(e/a). The source is Titchmarsh's local factorization (Thm 9.6(B)) plus the zero count (Thm 9.2). Off-line zeros only help.
  - Consequence: thickness at level T^{−O(log(e/a))}, instead of T^{−C₀/a}.
  - **Theorem H.1′ (PROVED (sketch)):** κ_Λ ≤ Λ^{B′}, B′ absolute. The η-part of the exponent drops from O(1/a) to O(log(e/a)), but the shape stays polynomial.
  - The other regime, a ≍ 1/log Λ, gives exp(O(log Λ · log log Λ)), which is worse.
  - **Λ^{o(1)} is not reachable by this route.** The O(log t) error in 9.6(B) and the zero count both cost T^{−c}.
  - Fake-blind: DH and Epstein (PROVED (sketch)).
- (H.9) **Trimmed space.** Let I^{(K)} := {n ∈ I_Λ : 1 − λ_n ≤ Λ^{−K}} and K ≥ 2B + 8, where B is any exponent with κ_Λ ≤ CΛ^B (H.1/H.1′). Then **ℓ_Λ^{(K)} ≤ CΛ^{−1}, so M1 holds for the trimmed space, diagonal included, at rate O(Λ^{−1/2} log^{1/2}Λ).**
  - Label: PROVED-conditional on Theorem H.1/H.1′ (PROVED (sketch)) and on the GAP-minor extension of Lemma M1.1. **No Spectrum input is needed for this.**
  - Every archimedean input ((T1), (T2a–c)) is proved here. They come from KRD as verified in Spectrum's PROLATE_BOUNDS §5, from an elementary endpoint lemma, and from Bernstein.
  - What is left is the **trimming-harmless check (T3)**: Tr(S_Λ ΔQ S_Λ V(h)) = o(1), with ΔQ = Q′_Λ − Q′^{(K)}_Λ, of rank O(K log²Λ). It is exactly equivalent to Q.8.C ⇔ Q.8.C^{(K)}, and it is the M2-compatibility condition. Split: (T3a) archimedean low-Mellin mass of plunge states (Spectrum; HEURISTIC) and (T3b) E-transport (Queen; GAP).
  - Correction to H.5: trimming removes O(K log²Λ) plunge states, not O(K log Λ).
- (H.10) (Z-avg⁻) under RH: one flagged line only.

## H.8 Lemma H.2′ and Theorem H.1′ (window-average form)

**Citation (PROVED (citation)).** E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed., revised by D. R. Heath-Brown (Oxford, 1986).
- §9.6, Theorem 9.6(B): log ζ(s) = Σ_{|t−γ|<1} log(s − ρ) + O(log t), **uniformly for −1 ≤ σ ≤ 2** (t ≥ 2), summed over zeros ρ = β + iγ with |t − γ| < 1.
- §9.2, Theorem 9.2: N(T + 1) − N(T) = O(log T).
- *What I checked.* The form and uniformity above are as I know them from §9.6, where (B) is obtained by integrating 9.6(A), ζ′/ζ(s) = Σ_{|t−γ|<1} 1/(s−ρ) + O(log t), along σ from 2. The exact page and the strict-vs-non-strict "<1" were **not re-verified against the book this session** (no copy on the box). Both are immaterial below: only the **real part** is used, which is free of branch conventions, and the zero set is enlarged to a fixed S₀ anyway.
- Real part: log|ζ(σ+it)| = Σ_{|t−γ|<1} log|σ + it − ρ| + O(log t), uniformly for −1 ≤ σ ≤ 2, t ≥ 2, including at zeros (where it reads −∞ = −∞).

**Lemma H.2′ (window average; PROVED (citation: Titchmarsh 9.6(B), 9.2), elementary given these).** There are absolute C₁, T₀ with the following property. Let T ≥ T₀, |t₀| ∈ [3, T], 0 < a ≤ 1, and W := [t₀ − a/2, t₀ + a/2]. Then

  (1/a) ∫_W log|ζ(½ + it)| dt ≥ −C₁ log T · log(e/a).

*Proof.* By conjugation take t₀ > 0.
1. **Fix the zero set.** S₀ := {ρ : |γ − t₀| < 2}. For t ∈ W, every zero with |t − γ| < 1 lies in S₀.
   - For ρ ∈ S₀ with |t − γ| ≥ 1: 1 ≤ |½ + it − ρ| ≤ ½ + 5/2 = 3, so log|½ + it − ρ| ∈ [0, log 3].
   - #S₀ ≤ N(t₀ + 2) − N(t₀ − 2) ≤ C log T (Thm 9.2, four times; for t₀ ≤ T₀ + 2 it is bounded).
   - So, uniformly in t ∈ W: log|ζ(½+it)| = Σ_{ρ∈S₀} log|½ + it − ρ| + θ(t)·C log T with |θ| ≤ 1. The added terms are ≤ #S₀ log 3 = O(log T).
2. **Off-line zeros only help.** |½ + it − ρ| = |(½ − β) + i(t − γ)| ≥ |t − γ|, for any β.
3. **One zero.** For c ∈ ℝ, g(c) := (1/a)∫_{−a/2}^{a/2} log|x − c| dx satisfies g′(c) = (1/a)(log|c + a/2| − log|c − a/2|). This has the sign of c, so g is minimized at c = 0, where g(0) = log(a/2) − 1 = −log(2e/a).
   Hence (1/a)∫_W log|½ + it − ρ| dt ≥ −log(2e/a) for every ρ (a multiple zero counts with multiplicity, and each copy pays the same).
4. Summing: average ≥ −#S₀ · log(2e/a) − C log T ≥ −C₁ log T · log(e/a), using log(2e/a) ≤ 2 log(e/a) since log(e/a) ≥ 1. ∎

R-line H.2′: R1 window-uniform, no L | R2 Euler product only through Thm 9.6's anchor |ζ(2+it)| ≫ 1; FE through growth/zero count | R3 n/a | R4 n/a | R5 fake-blind (below) | R6 n/a.

**Corollary H.2′′ (thickness; PROVED).** Let U := sup_{|t|≤T+1} log|ζ(½+it)| ≤ (¼ + o(1)) log T (convexity) and V := C₁ log T · log(e/a). Put

  η′_a(T) := exp(−U − 2V) = T^{−¼ − 2C₁ log(e/a) − o(1)}.

Then E′_a(T) := {s : |ζ(½+is)| ≥ η′_a(T)} ∪ {|s| ≤ 4} ∪ {|s| > T} is (½, a)-thick.
*Proof.* Take a window W of centre t₀ with 3 ≤ |t₀| ≤ T + ½ (other windows lie in {|s| ≤ 4} ∪ {|s| > T}, as in the audit note to Lemma H.2), and apply H.2′ with T + 1. Let F := {t ∈ W : log|ζ| < −Y} and θ := |F|/a. Then −V ≤ avg_W log|ζ| ≤ (1 − θ)U − θY, so θ ≤ (U + V)/(U + Y). With Y := U + 2V this gives θ ≤ ½. ∎
Compare Lemma H.2: there the level was T^{−C₀/a}. Here it is T^{−O(log(e/a))}, so the dependence on a goes from 1/a to log(1/a).

**Theorem H.1′ (PROVED (sketch); constants untracked, like H.1).** m_Λ ≥ cΛ^{−B′} for all large Λ, with B′ absolute. Same proof as Theorem H.1, with E′_a(T), η′_a(T) in place of E_a(T), η_a(T). Bookkeeping of the exponent (no explicit Kovrijkine constants):
- Step 2 (Kovrijkine) gives K = Λ^{e_K}·O(1), where e_K := C a(2 + A₁) log(2C).
- Step 3 (left-tail/remainder) requires e_K < 5A₁/4. This is the only constraint, and it forces a < a* := 5/(4C log(2C)) (absolute) for every A₁.
- Steps 4–5: T = C′Λ²K² = Λ^{2+2e_K+o(1)}, and

  m_Λ ≥ η′_a(T)²(4K)^{−2} = Λ^{−B′}, B′ = 2e_K + (2 + 2e_K)(½ + 4C₁ log(e/a)) + o(1).

**Shape (the question asked).**
- **a fixed (0 < a < a*, absolute; A₁ fixed with e_K < 5A₁/4):** κ_Λ ≤ Λ^{B′}, with B′ = O(1 + log(e/a)) at fixed e_K. This is the best shape from this route.
  - Versus H.1: the η-part of B goes from O(1/a) to O(log(e/a)). The Kovrijkine part e_K is unchanged.
  - Taking a small now costs only log(1/a), so one may take a ≪ a* and A₁ = O(1) with e_K ≍ a·A₁ small. Then B′ = 1 + O(log(e/a)) + O(e_K).
  - "B′ absolute" is meant as before: absolute, not explicit.
- **a ≍ 1/log Λ (with A = O(1), so aℓ = O(1) and K = O(1)):** the remainder constraint becomes e^{−5A/2} ≪ 1 (fine), but η′ = T^{−O(log log Λ)} with T ≍ Λ². So m_Λ ≥ exp(−O(log Λ · log log Λ)). This is **worse** than polynomial, so a ~ 1/log Λ is not the right choice.
- **Λ^{o(1)} is not reachable by the factorization route (PROVED (elementary), as a statement about the method).**
  - Its inputs allow a background loss of size O(log T) in 9.6(B).
  - The count #S₀ ≍ log T allows c log T zeros (with multiplicity) to sit in one window of length a. That configuration is not excluded by 9.2 and 9.6 alone.
  - Either effect makes the lower bound on avg log|ζ| as small as −c log T, i.e. a loss of T^{−c} = Λ^{−2c}.
  - Reaching Λ^{o(1)} needs average information beyond counting: (Z-avg⁻) of H.4 (GAP), or an equivalent zero-clustering bound.

**Fake-blind check (R5; PROVED (sketch)).** The proof of Thm 9.6 uses three inputs:
- (i) an anchor where the function is ≫ 1;
- (ii) polynomial growth in a fixed disc around the anchor;
- (iii) the Borel–Carathéodory / Jensen lemma (Titchmarsh §3.9, Lemma α type) on that disc. It gives ζ′/ζ = Σ_{|ρ−s₀|≤r/2} 1/(s−ρ) + O(M/r) on |s − s₀| ≤ r/4, with log of the growth ratio ≤ M = O(log t), and also the count of zeros in the disc ≪ M = O(log t).

These carry over to the fakes as follows.
- **Davenport–Heilbronn.**
  - Anchor s₀ = 3 + it₀: |f_DH(3+it)| ≥ 1 − (ζ(3) − 1) > 0.79, as in H.6.
  - f_DH is entire with |f_DH| ≪ |t|^{A} on −2 ≤ σ ≤ 4 (FE + Phragmén–Lindelöf).
  - Take r = 16, so the disc |s − s₀| ≤ r/4 = 4 covers −1 ≤ σ ≤ 2. Integrating from the anchor, then trimming to |t − γ| < 1 exactly as in Thm 9.6, gives log|f_DH(s)| = Σ_{|t−γ|<1} log|s − ρ| + O(log t). The zero count O(log t) comes from Jensen on the same disc.
  - Steps 1–4 of H.2′ then go through verbatim. Step 2 is where DH's **off-line zeros** enter, and they only increase |s − ρ|.
  - So Lemma H.2′, Corollary H.2′′ and Theorem H.1′ hold for DH. Their constants depend on DH's growth exponent, not on its zeros' real parts.
- **Epstein ζ_Q (x² + 5y²).**
  - Anchor at σ₀ large: Z_Q(s) = Σ r_Q(n) n^{−s} with r_Q(1) = 2, and |Z_Q(σ₀ + it)/2 − 1| ≤ ½ for σ₀ ≥ σ_Q.
  - Polynomial growth comes from the FE and Phragmén–Lindelöf. The pole at s = 1 is avoided by |t₀| ≥ 3, as for ζ.
  - Same conclusion: PROVED (sketch).
- Hence the window-average form is **fake-blind**, so it carries no RH content. That is consistent with H.6.

## H.9 Trimmed prolate space: the interface with Spectrum (and what is already proved)

**Conventions.**
- (ψ_n) are the even prolates of M1.4: L²[−Λ, Λ]-normalized truncations of whole-line prolates g_n.
- g_n is band-limited to |ξ| ≤ Λ, with ‖g_n‖²_{L²(ℝ)} = 1/λ_n by double orthogonality.
- Spectrum's normalization: c = 2πΛ², M = 4Λ², β_Λ := (2/π²) log(200Λ² + 25) = (4/π² + o(1)) log Λ.
- I_Λ is Connes' retention (λ_n ≥ ½ − o(1)).
- Define I^{(K)} := {n ∈ I_Λ : 1 − λ_n ≤ Λ^{−K}}, B^{(K)}_{Λ,0} := span{ψ_n : n ∈ I^{(K)}} ∩ {φ(0) = φ̂(0) = 0} ⊂ B_{Λ,0}, and Q′^{(K)}_Λ := projection onto E(B^{(K)}_{Λ,0}), with ΔQ := Q′_Λ − Q′^{(K)}_Λ ≥ 0.
- Write κ^{(K)}, ℓ^{(K)}, ε^{(K)} for the trimmed quantities.

**Requested inputs, with owners and labels.**

| Input | Statement | Label | Owner |
|---|---|---|---|
| (T0) | κ^{(K)} ≤ κ_Λ ≤ CΛ^B | PROVED (H.0 monotonicity) + Thm H.1/H.1′ PROVED (sketch) | Queen (done) |
| (T1) count | r_K := #(I_Λ \ I^{(K)}) ≤ β_Λ(K log Λ + log 10) + 7 = (4/π² + o(1)) K log²Λ; I^{(K)} ⊇ {even k ≤ ⌊M⌋ − 7 − β_Λ(K log Λ + log 10)} | PROVED (citation: KRD Thm 3 / Cor 3, as verified in Spectrum PROLATE_BOUNDS §5) | Spectrum (done; quoted) |
| (T2a) endpoints | |ψ_n(Λ)|² ≤ (64π)^{2/3} Λ (1 − λ_n)^{2/3} for λ_n ≥ ½; so ε^{(K)} := Σ_{I^{(K)}} |ψ_n(Λ)|² ≤ CΛ^{3−2K/3} | PROVED (elementary; Lemma H.9.1 below) | Queen (done). *Optional* for Spectrum: exponent 1 in place of 2/3 via ∂_cλ_n = (2/c)λ_nψ_n(1)² (would relax only the non-binding endpoint constraints; not needed) |
| (T2b) tails | ∫_{|ξ|>X} |Fψ_n|² ≤ ∫_{|ξ|>Λ} |Fψ_n|² = 1 − λ_n ≤ Λ^{−K} for X ≥ Λ | PROVED (M1.1 bullet; definition of I^{(K)}) | done |
| (T2c) smoothness | Σ_{n∈I^{(K)}} q₂(ψ_n) ≤ CΛ^7, with q₂(φ) := |φ′(Λ)|² + Λ‖φ″‖²_{L²[−Λ,Λ]} | PROVED (Bernstein; below) | done |
| (T3a) archimedean low-Mellin mass of the removed states | ω^{arch}_Λ(S) := sup over unit φ ∈ span{ψ_n : n ∈ I_Λ \ I^{(K)}} of ∫_{|s|≤S} |Mφ|² / ∫|Mφ|² ≤ Λ^{−B−2−δ} for S = log Λ (any polynomial rate faster than Λ^{−B−2} suffices) | HEURISTIC (plunge states live in phase space near |x| ≈ Λ, |ξ| ≈ Λ, i.e. at Mellin height s ≍ 2πΛ², H.1b) | **Spectrum (requested; RH-silent)** |
| ~~(T3a)~~ → (T3a*) *(Part III)* | **Withdrawn**: the rate Λ^{−B−2−δ} and the edge-localization premise are expected FALSE (Spectrum 65bba74: plunge states oscillate at |ξ| ≈ Λ at all x; low-Mellin fraction ≍ (log Λ)/Λ²). Replaced by the trace-form (T3a*), 𝒟_rem(s) ≤ C r_K(1+|s|)^AΛ^{−2}, H.13.3 | CONJECTURE (Spectrum: HEURISTIC + NUMERICAL) | Spectrum (after Arb work) |
| (T3b) E-transport | the E-orthogonal complement of B^{(K)}_{Λ,0} in B_{Λ,0}, and its S_Λ-cut, inherit (T3a) up to a factor Λ^{O(1)} | **GAP** | **Queen** |
| (T3b) *(Part III)* | reduced to (PD) (profile domination; GAP) + (L) (HEURISTIC); worst case closes iff B < 2 | see H.13 | Queen |

**Lemma H.9.1 (endpoint value from out-of-band mass; PROVED (elementary)).** If λ_n ≥ ½, then |ψ_n(Λ)|² ≤ (64π)^{2/3} Λ (1 − λ_n)^{2/3}.
*Proof.* Put g := g_n and v := |g(Λ)| = |ψ_n(Λ)|.
- Bernstein and Cauchy–Schwarz: ‖g‖_∞ ≤ ‖ĝ‖_1 ≤ (2Λ)^{1/2}‖g‖_2 ≤ (2Λ)^{1/2}λ_n^{−1/2} ≤ 2Λ^{1/2}, and ‖g′‖_∞ ≤ 2πΛ‖g‖_∞ ≤ 4πΛ^{3/2} =: D.
- So |g| ≥ v/2 on [Λ, Λ + v/(2D)], and ∫_{x>Λ}|g|² ≥ v³/(8D).
- But ∫_{|x|>Λ}|g|² = 1/λ_n − 1 ≤ 2(1 − λ_n). Hence v³ ≤ 16D(1 − λ_n) = 64πΛ^{3/2}(1 − λ_n). ∎
- *Grok audit (S103g), checked line by line.*
  - Under ĝ(ξ) = ∫g e^{−2πixξ} with band |ξ| ≤ Λ (so c = 2πΛ²), ‖g‖_∞ ≤ ‖ĝ‖₁ ≤ (2Λ)^{1/2}‖g‖₂, and ‖g‖²₂ = 1/λ_n ≤ 2.
  - Bernstein gives ‖g′‖_∞ ≤ 2πΛ‖g‖_∞. Since |g| is D-Lipschitz, the bound |g| ≥ v/2 on a length v/(2D) holds.
  - ∫_{|x|>Λ}|g|² = (1 − λ_n)/λ_n.
  - The exponent: (16D)^{2/3} = (64π)^{2/3}Λ. PROVED (elementary) stands.
- With |I^{(K)}| ≤ |I_Λ| ≤ 4Λ² + O(log Λ) (Spectrum §4, index location; the even states are about half of these), this gives ε^{(K)} ≤ CΛ^{3−2K/3}.

**(T2c) proof.** ‖g′‖²_∞ ≤ 16π²Λ³ (above), and ‖g″‖²_{L²(ℝ)} ≤ (2πΛ)⁴‖g‖²_2 ≤ 32π⁴Λ⁴. So q₂(ψ_n) ≤ CΛ⁵ for each n, and summing over ≤ 5Λ² indices gives CΛ⁷. ∎

**Theorem H.9 (M1 on the trimmed space; PROVED-conditional on Thm H.1/H.1′ (PROVED (sketch)) and on the Lemma M1.1 extension (GAP-minor)).** Let κ_Λ ≤ CΛ^B and K ≥ 2B + 8. Then ℓ^{(K)}_Λ := ‖(1 − S_Λ)Q′^{(K)}_Λ‖²_HS ≤ CΛ^{−1}. Hence, by Theorem M1.A (which holds for any finite-rank projection Q′), for every f ∈ L¹ ∩ L²:

  |Lk^{(K)}_Λ(f ∗ f^*)| ≤ C‖f‖₁² Λ^{−1} + C‖f‖₁‖f‖₂ Λ^{−1/2}(log Λ)^{1/2} → 0, including on the diagonal Λ = e^L.

*Proof.* Run the proof of Theorem M1.C on B^{(K)}_{Λ,0} with σ = ½, but bound each remainder trace by the **smaller of two positive forms**. For X ≥ Λ, by Lemma M1.2:
- Q1_X(φ) := 2∫_{|ξ|>X}|Fφ|² + (4/(π²X))|φ(Λ)|² ≥ ∫_{|ξ|>X}|R_φ|², because R = Fφ − jump and ∫_{|ξ|>X}|jump|² ≤ 2|φ(Λ)|²/(π²X).
- Q2_X(φ) := q₂(φ)/(3π⁴X³) ≥ ∫_{|ξ|>X}|R_φ|², because (|φ′(Λ)|/2 + ‖φ″‖₁/4)² ≤ q₂(φ)/2 and ∫_{|ξ|>X} ξ^{−4} = 2/(3X³).

Both are positive forms. So, as in M1.C, Σ_j ∫_{|ξ|>X}|R_{φ_j}|² ≤ κ^{(K)} min(Tr Q1_X, Tr Q2_X). By (T2a–c):
- Tr Q1_X ≤ C(Λ^{2−K} + Λ^{3−2K/3}/X);
- Tr Q2_X ≤ CΛ^7/X³.

With X = mΛ, split at m = Λ^p, p := 2(B + 5)/3:

  Σ_m m^{1/2} min(·, ·) ≤ C(Λ^{2−K+3p/2} + Λ^{2−2K/3+p/2} + Λ^{4−3p/2}).

The jump part is 2κ^{(K)}ε^{(K)}/Λ ≤ CΛ^{B+2−2K/3}. Multiplying by κ^{(K)} ≤ CΛ^B, the four exponents are:
- 2B + 7 − K;
- (4B + 11 − 2K)/3 *(Grok audit S103g: was (4B + 8 − 2K)/3. Since p/2 = (B + 5)/3, the exponent is B + 2 − 2K/3 + p/2 = (4B + 11 − 2K)/3. This is still ≤ −5/3 at K = 2B + 8, so the conclusion is unchanged.)*;
- −1;
- B + 2 − 2K/3.

All are ≤ −1 once K ≥ 2B + 8. ∎

- B is absolute but not explicit, so K(B) = 2B + 8 is absolute, not explicit. Spectrum's optional exponent-1 endpoint bound would **not** change K(B): the binding constraint, 2B + 7 − K ≤ −1, comes from the tail term (T2b), not from the endpoints.
- R-line H.9: R1 Λ-uniform, no L | R2 E's ℤ-sum only via the sawtooth (as in M1) | R3 n/a | R4 φ(0) = φ̂(0) = 0 kept | R5 fake-blind: the archimedean inputs are identical, and κ for DH/Epstein comes from H.1′ | R6 n/a.

**(T3): trimming-harmless check, and M2-compatibility.**
- *Exact reduction (PROVED, algebra).* Since Q′ = Q′^{(K)} + ΔQ, we have Tr(S Q′ S V(h)) − Tr(S Q′^{(K)} S V(h)) = Tr(S ΔQ S V(h)) = ‖V(f)^* S ΔQ‖²_HS ≥ 0. Define
  **Conj Q.8.C^{(K)}** (CONJECTURE): Tr(S_Λ Q′^{(K)}_Λ S_Λ V(h)) = 2‖f‖² log Λ − 𝒮(f) + o(1). Then
  **Q.8.C ⇔ Q.8.C^{(K)} ⇔ (T3): Tr(S_Λ ΔQ S_Λ V(h)) = o(1) for each fixed f ∈ C_c^∞.** (Each equivalence assumes the other side's statement; precisely, Q.8.C and Q.8.C^{(K)} differ by exactly this term.) *Grok audit (S103g): read the chain as "given Q.8.C, Q.8.C^{(K)} ⇔ (T3); given Q.8.C^{(K)}, Q.8.C ⇔ (T3)". The three statements are not mutually equivalent: (T3) alone implies neither conjecture. The algebra is exact, Tr(SΔQSV(h)) = ‖ΔQ S V(f)‖²_HS ≥ 0, because Q′^{(K)} ≤ Q′ are nested projections.*
- *M2-compatibility (PROVED, same algebra).* S ⊖ S E(B^{(K)}) is S ⊖ S E(B) plus a piece of rank ≤ r_K, and its V(h)-trace exceeds M2's by exactly Tr(S ΔQ S V(h)). So M2 for the trimmed space has the **same** counting term 2‖f‖² log Λ − PV and the same zero-vector Gram/Poisson sum **iff (T3)**. M3 is unaffected, since the zero vectors live at fixed heights.
- *Sufficient condition (PROVED, elementary).* For x ∈ L²(du), ‖V(f)^*x‖² = (1/2π)∫|f̂(s)|²|x̂(s)|² ds, and |f̂(s)| ≤ C_N(1+|s|)^{−N}. So, for S = log Λ,
  Tr(S ΔQ S V(h)) ≤ r_K (‖f‖₁² ω_Λ(S) + C_N S^{−2N}), where ω_Λ(S) := sup over unit x ∈ ran(S_Λ ΔQ) of (1/2π)∫_{|s|≤S}|x̂|².
  With r_K = O(K log²Λ), **(T3) ⇐ ω_Λ(log Λ) = o(1/log²Λ).**
- *Reduction of ω_Λ to (T3a) + (T3b).*
  - If x were E(φ) with φ ∈ span{removed ψ_n}, then x̂ ∝ ζ(½+is)Mφ(s). By the definition of m_Λ, ∫|ζ|²|Mφ|² ≥ m_Λ∫|Mφ|² ≥ cΛ^{−B}∫|Mφ|².
  - With max_{|s|≤log Λ}|ζ(½+is)|² ≤ C log Λ, this gives ω ≤ CΛ^{B} log Λ · ω^{arch}_Λ(log Λ), which is o(1/log²Λ) under (T3a).
  - The two real issues are:
    - **(T3b)**: ran ΔQ is the E-orthogonal complement of E(B^{(K)}_0) inside E(B_0). Its φ's are not in the removed span; they are mixed with kept states through the E-Gram matrix.
    - The S_Λ-cut, which smears s at scale 1/log Λ (harmless at height log Λ, but must be written).
  - Owner Queen; GAP.
- **Recommendation (OPINION).** Make the trimmed space the primary object of the program: Q.8.C^{(K)}, M2^{(K)}, M3^{(K)}. M1^{(K)} is now proved (modulo sketch/GAP-minor). (T3) is then needed only to identify the limit with Connes' own Q′_Λ.
  - Reason: for Connes' untrimmed retention, M1 is **not** expected to follow from M1.A. Spectrum §5′ shows that the archimedean HS² leak of the full prolate span, with N pushed to the plunge, is ≍ (2/π²) log Λ, not small.
  - M1 itself may fail at O(1) there. This is OPINION; there is no lower bound for ℓ_Λ (untrimmed) here.

## H.10 (Z-avg⁻) under RH (flagged; one line only)

**Under RH (flagged: barred from every route to RH by the circularity flag, and untestable on DH/Epstein, which violate RH):** (Z-avg⁻) keeps the status of H.4. It is CONJECTURE, expected via the RH-conditional Selberg/Soundararajan-type explicit formulae (HEURISTIC). Nothing in H.8–H.9 uses it.

## H.11 Part II summary table

| Claim | Label |
|---|---|
| Titchmarsh 9.6(B) real-part form, uniform on −1 ≤ σ ≤ 2; 9.2 zero count | PROVED (citation; page/strictness not re-verified this session, immaterial) |
| Lemma H.2′ window average ≥ −C₁ log T log(e/a); off-line zeros only help | PROVED (citation + elementary) |
| Cor H.2′′ (½, a)-thickness at level T^{−¼−2C₁ log(e/a)} | PROVED |
| Thm H.1′ κ_Λ ≤ Λ^{B′}, B′ absolute; η-part O(log(e/a)) instead of O(1/a) | PROVED (sketch; constants untracked) |
| a ≍ 1/log Λ gives only exp(O(log Λ log log Λ)) | PROVED (sketch; as a statement about the method) |
| Λ^{o(1)} unreachable via 9.6(B) + 9.2 alone | PROVED (elementary; method limitation) |
| H.2′/H.1′ for DH and Epstein (fake-blind) | PROVED (sketch) |
| (T1) r_K ≤ (4/π² + o(1)) K log²Λ (corrects H.5's O(K log Λ)) | PROVED (citation: KRD via Spectrum §5) |
| Lemma H.9.1 / (T2a) endpoint bound, ε^{(K)} ≤ CΛ^{3−2K/3} | PROVED (elementary) |
| (T2b), (T2c) | PROVED |
| Thm H.9: K ≥ 2B + 8 ⇒ ℓ^{(K)} ≤ CΛ^{−1}; M1 for the trimmed space, diagonal included | PROVED-conditional on Thm H.1/H.1′ (sketch) and the Lemma M1.1 extension (GAP-minor) |
| Q.8.C ⇔ Q.8.C^{(K)} ⇔ (T3); M2-compatibility ⇔ (T3) | PROVED (algebra) |
| (T3) ⇐ ω_Λ(log Λ) = o(1/log²Λ) | PROVED (elementary) |
| (T3a) plunge-state low-Mellin mass ≤ Λ^{−B−2−δ} | **withdrawn (Part III)**: expected FALSE per Spectrum 65bba74; replaced by (T3a*) at rate r_K/Λ² |
| (T3b) E-transport of (T3a) to ran ΔQ and its S_Λ-cut | GAP (Queen) |
| (T3b) *(Part III)* | reduced to (PD) (profile domination; GAP) + (L) (HEURISTIC); worst case closes iff B < 2 | see H.13 | Queen |
| Conj Q.8.C^{(K)} | CONJECTURE (unconditional; not RH-equivalent) |
| M1 for Connes' untrimmed space may fail at O(1) | OPINION |
| (Z-avg⁻) under RH | flagged one-liner; barred; CONJECTURE/HEURISTIC |

Not cited: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.

---

# Part III (2026-10-28): (T3b) E-transport, profile-weighted; density extension of Lemma M1.1

Against GPPDiscovery5 main `9cd3187` (Part II audited SOUND, merged at f847b98; the lead's fixes are kept). Analytic only: **no numerics were run by Queen.** RH is OPEN. Nothing here assumes RH.
Spectrum input: team/spectrum-measure `65bba74`, T3a_plunge_mellin.md (credit: God of the Spectrum & Measure Theory). (T3a) itself is Spectrum's job, after their Arb work. Below I only **state the interface I need** from it.

**Part III headline.**
- **(B) closed.** Lemma M1.1 holds on all of B_{Λ,0} (truncated prolates, with φ(0) = φ̂(0) = 0). Mellin–Plancherel and κ are well defined there. PROVED (H.12; Dirichlet–Jordan + Mellin mollification).
  - So Theorem H.9 (M1 for the trimmed space) is now **PROVED-conditional only on Theorem H.1/H.1′ (PROVED (sketch))**. The GAP-minor is gone.
- **(T3a) corrected.** My Part II request (T3a): ω^{arch} ≤ Λ^{−B−2−δ}, together with its HEURISTIC premise (plunge states localized near |x| ≈ Λ), is **withdrawn**.
  - Spectrum 65bba74 shows that plunge states oscillate at |ξ| ≈ Λ at *all* x. Their model |M cos(c·)|² = π/(2c) is PROVED. Their low-Mellin fraction is ≍ (log Λ)/Λ² (CONJECTURE: HEURISTIC + NUMERICAL).
  - New interface **(T3a*)**, a basis-free **Mellin density of states** of the removed span at the corrected rate r_K/Λ².
- **(T3b) profile-weighted transport (H.13).** **Theorem H.13 (PROVED-conditional):** (T3a*) + (PD) + (L) ⇒ Tr(S_Λ ΔQ S_Λ V(h)) = O(K (log Λ)^{A+2} Λ^{−2}) + o(1) → 0, **with no Λ^B loss and no B < 2 condition**.
  - (L) is the low-pass of the leak. Its jump (sawtooth) part is PROVED; this is exactly where the integers enter. The rest of (L) is HEURISTIC.
  - **(PD), profile domination of the E-orthonormalized complement, is the one remaining GAP (Queen).** This is where the ζ-weight enters. It is supported by fixed-weight mean values (Hardy–Littlewood / Balasubramanian–Ramachandra; PROVED pieces H.13.3) rather than by the worst-case m_Λ.
  - The worst-case route (PROVED) closes only if B < 2, matching Spectrum's T4.
- Fake-blind (DH, Epstein): every ingredient is a second-moment, Poisson/FE, or archimedean statement. So this stays **RH-silent**.

## H.12 (B) Density extension of Lemma M1.1 to B_{Λ,0}: PROVED

Setting: φ ∈ B_{Λ,0}, i.e. even, supported in [−Λ, Λ], smooth on [−Λ, Λ] (a finite combination of truncated prolates), with a jump of size φ(Λ) at ±Λ, φ(0) = 0 and φ̂(0) := ∫φ = 0. Write E(φ)(u) := e^{u/2} Σ_{n≠0} φ(n e^u). For each u this is a finite sum, defined at all u with Λe^{−u} ∉ ℤ, using midpoint values at jumps.

**Step 1 (the pointwise Poisson identity, without density; PROVED).** For x > 0 with Λ/x ∉ ℤ:

  Σ_{n≠0} φ(nx) = x^{−1} lim_{N→∞} Σ_{0<|k|≤N} Fφ(k/x).

*Proof.* g(t) := Σ_{n∈ℤ} φ(x(t + n)) is 1-periodic and of bounded variation, since locally it is a finite sum of BV functions. Its Fourier coefficients are ∫_ℝ φ(xy)e^{−2πiky}dy = x^{−1}Fφ(k/x).
- By the Dirichlet–Jordan test (Zygmund, *Trigonometric Series*, Vol. I, Ch. II §8; theorem number not re-verified this session), the symmetric partial sums at t = 0 converge to ½(g(0+) + g(0−)).
- This equals g(0) unless some nx = ±Λ, i.e. unless Λ/x ∈ ℤ.
- The **pole terms** are the n = 0 and k = 0 terms, φ(0) and x^{−1}Fφ(0) = x^{−1}φ̂(0). Both vanish by the constraints. ∎

With x = e^u this is E(φ)(u) = E(Fφ)(−u) (symmetric summation) for a.e. u. That is Lemma Q.8.2 on B_{Λ,0}, with **no density argument needed** for the identity itself.

**Step 2 (the leak formula; PROVED).** For u < −log Λ, substitute y = e^{−u}:

  ‖(1 − S_Λ)E(φ)‖² = ∫_{−∞}^{−log Λ} |e^{−u/2} Σ_{n≠0} Fφ(n e^{−u})|² du = ∫_Λ^∞ |Σ_{n≠0} Fφ(ny)|² dy,

with the symmetric sum. Lemma M1.2 then splits it into:
- the sawtooth part: boundedly convergent, Σ_{n≥1} sin(2πnΛy)/n = π(½ − {Λy}) for Λy ∉ ℤ;
- an absolutely convergent remainder, |R(ξ)| ≪ ξ^{−2}.

So the right side is finite, and Lemma M1.1 holds as stated on B_{Λ,0}.

**Step 3 (consistency with the closure of E, and Mellin–Plancherel on B_{Λ,0}; PROVED).**
- **Decay of Mφ.** F_φ(u) = e^{u/2}φ(e^u) has F_φ′ ∈ L¹ apart from one jump of size Λ^{1/2}φ(Λ) at u = log Λ. So |Mφ(s)| ≤ (‖F_φ′‖₁ + Λ^{1/2}|φ(Λ)|)/|s|.
- With ∫_{|s|≤T}|ζ(½+is)|² ≪ T log T (Hardy–Littlewood; Titchmarsh §7.3), this gives ‖φ‖²_ζ := (1/2π)∫|ζ(½+is)|²|Mφ(s)|² ds < ∞.
- **Approximation.** Let ρ_k ∈ C_c^∞(ℝ_+^*) be a multiplicative approximate identity (∫ρ_k d^*y = 1, supp ρ_k → {1}), and put φ_k(x) := ∫ φ(x/y)ρ_k(y) d^*y.
  - φ_k is even and C_c^∞; the jump is smoothed. Since φ is smooth near 0, φ_k is smooth at 0.
  - φ_k(0) = φ(0) = 0 and ∫φ_k = φ̂(0)·∫yρ_k(y)d^*y = 0. So φ_k ∈ S_0^ev.
  - Mφ_k = Mφ · Mρ_k with |Mρ_k| ≤ 1 and Mρ_k → 1 pointwise. By dominated convergence, ‖φ_k − φ‖_ζ → 0.
- **Identification of the limit.** On S_0^ev, ‖E(·)‖² = 2‖·‖²_ζ-type (M1.4; Connes 1999 §III), so E(φ_k) is Cauchy in L²(du).
  - Its limit equals the pointwise E(φ) a.e.: E(φ_k)(u) → E(φ)(u) at every u with Λe^{−u} ∉ ℤ, since these are finite sums and φ_k → φ at continuity points. Pass to an a.e.-convergent subsequence.
  - Hence the closure of E coincides with the pointwise E on B_{Λ,0}. Mellin–Plancherel, ‖E(φ)‖² = (2/π)∫|ζ(½+is)|²|Mφ(s)|² ds (M1.4 normalization, using ∫|Mφ|² = π‖φ‖²), holds on B_{Λ,0}. So does the definition of κ_Λ.
- **Continuity of both sides of M1.1 in ‖·‖_ζ (the alternative route).**
  - LHS: ‖(1 − S)E(·)‖ is ‖·‖_ζ-continuous.
  - RHS: by Step 1, the RHS equals the LHS identically, so it is continuous as well. No separate continuity estimate is needed.
- **Status.** The GAP-minor in M1 (Lemma M1.1 extension) is **closed: PROVED**.

**Effect.** Theorem H.9 (M1^{(K)} for K ≥ 2B + 8) becomes PROVED-conditional on Theorem H.1/H.1′ only. Since those are PROVED (sketch), M1^{(K)} is **PROVED (sketch)**.
R-line H.12: R1 Λ-uniform | R2 Poisson (FE half) only, through the ℤ-sum | R3 n/a | R4 the pole terms φ(0) and φ̂(0) are exactly what vanish | R5 DH: periodic a_n, so the Poisson dual is the FE-dual periodic sequence (Gauss sums), and the same BV argument applies; Epstein: 2-dimensional Poisson on the lattice, BV in each fibre | R6 exact identity.

## H.13 (T3b) E-transport, profile-weighted

### H.13.0 Prior art (checked in local copies: /workspace/queen/papers)
- **Connes–Consani, arXiv 2106.01715 ("Spectral triples and ζ-cycles"), §3 "Eigenfunctions and the prolate projection Π(λ,k)".**
  - E on S_0^ev, (3.1), and the Poisson relation E(f̂)(x) = E(f)(x^{−1}), (3.2), quoted: "E(f̂)(x) = E(f)(x^{−1}), ∀f ∈ S_0^ev".
  - Slepian–Pollak prolates as the "approximate intersection" P_λ ∩₀ P̂_λ, used to build the near-radical of QW_λ.
  - They give **no** conditioning, leak, or trimming estimates for E on the prolate span. (T3b) is not there.
- **Connes–Consani, arXiv 2006.13771 ("Weil positivity and trace formula, the archimedean place"; Selecta 2021), Lemma 1.4, Prop. 1.5 and the justification after (33).**
  - The **commutator trace-class technique**: "it is enough … to show that the commutator [P, ϑ_m(f)] is of trace class". Here ϑ_m(f) is multiplicative convolution, i.e. a Mellin multiplier f̂(s), and P is a cutoff. The kernel of [P, T] is (P(λ) − P(μ))k(λ, μ).
  - I use this in H.13.2 for [S_Λ, V(f)]. It is trace class with a Λ-independent norm, because the kernel is (1_S(u) − 1_S(v))f(u − v) at two boundary points. It is used **only as bookkeeping**: an O(1) bound, not o(1).
- **Connes–Consani–Moscovici, arXiv 2310.18423 ("Zeta zeros and prolate wave operators").**
  - Semilocal prolate operators, Sonin spaces, and the Jacobi-matrix form of W_λ in the Mellin picture. Abstract: the negative spectrum reproduces "the ultraviolet behavior" of the zeros, and the Sonin space corresponds to it.
  - This is consistent with Spectrum's finding that transition-band states spread across all Mellin heights. There is **no** E-conditioning estimate.
- **Conclusion.** No prior art bounds the E-transport of low-Mellin mass. (T3b) is new as stated.

### H.13.1 Exact structure (PROVED)
- Let R₀ := B_{Λ,0} ⊖_{L²} B^{(K)}_{Λ,0}, the L²-orthogonal complement inside B_{Λ,0}. Then dim R₀ = dim B_{Λ,0} − dim B^{(K)}_{Λ,0} ≤ r_K + 2.
- R₀ ⊂ span{removed ψ_n} ⊕ span{k₀, k₁}, where k₀ := Σ_{kept} ψ_n(0)ψ_n and k₁ := Σ_{kept}(∫ψ_n)ψ_n are the constraint-correction vectors in the kept span.
- Since E is injective, E(B_{Λ,0}) = E(B^{(K)}_{Λ,0}) ∔ E(R₀), so

  ran ΔQ = (1 − Q′^{(K)}) E(R₀) = E(R̃), R̃ := {φ ∈ B_{Λ,0} : Eφ ⊥ E(B^{(K)}_{Λ,0})}.

- Let (φ_i) be a basis of R̃ that is orthonormal for the **E-inner product**, so (Eφ_i) is an ONB of ran ΔQ. Then, exactly,

  Tr(S ΔQ S V(h)) = Σ_i ‖V(f)^* S E φ_i‖²  (H.13.a)

### H.13.2 Where each ingredient enters

Split S Eφ_i = Eφ_i − (1 − S)Eφ_i. Then (H.13.a) ≤ 2Σ_i‖V^* Eφ_i‖² + 2Σ_i‖V^*(1 − S)Eφ_i‖².

**(i) The ζ-weight enters here, exactly (PROVED, Mellin–Plancherel via H.12):**

  Σ_i ‖V(f)^* Eφ_i‖² = c₀ ∫ |f̂(s)|² |ζ(½+is)|² 𝒟(s) ds (c₀ absolute, M1.4 normalization), where 𝒟(s) := Σ_i |Mφ_i(s)|².

𝒟 is basis-free for the E-inner product on R̃: it is the diagonal of the kernel of the E-orthogonal projection, in Mellin form. Since f is fixed and |f̂| decays faster than any power, only |s| ≤ Λ^ε matters. The tail is ≤ C_N Λ^{−2Nε} Σ_i‖Eφ_i‖² = C_N Λ^{−2Nε}(r_K + 2).

**(ii) The integers (sawtooth, Lemma M1.2) enter here: the low-pass of the leak.** By H.12, for u < −log Λ, (1 − S)Eφ(u) = e^{−u/2} Σ_{n≠0}Fφ(ne^{−u}), the jump part plus the remainder.
- **Jump part (PROVED, elementary).** It equals 2φ(Λ)e^{u/2}(½ − {Λe^{−u}}). Its Fourier transform in u at frequency s is 2φ(Λ)Λ^{½−is} ∫_{Λ²}^∞ w^{−3/2+is}(½ − {w}) dw (substitute y = e^{−u}, then w = Λy).
  - ½ − {w} has mean 0 and a bounded primitive (≤ 1/8). One integration by parts gives modulus ≤ C|φ(Λ)|(1 + |s|)Λ^{−5/2}.
  - Hence ‖V^*·jump‖² ≤ C_f |φ(Λ)|² Λ^{−5}.
  - **This is the only use of ℤ, and it helps:** the sawtooth's mean zero kills the low-s content.
  - Summed over an L²-ONB this is ≤ C_f Λ^{−5} Σ|ψ_n(Λ)|² ≤ C_f Λ^{−2} (Lemma H.9.1 and |I_Λ| ≤ 5Λ²). Over the E-ONB (φ_i) there is an extra factor sup‖φ_i‖²_{L²}; see (PD).
- **Remainder part: input (L), HEURISTIC.** Σ_i ‖V^*·(remainder leak of Eφ_i)‖² = O(r_K Λ^{−2}(log Λ)^{O(1)}).
  - Reason: by Poisson the leak near u = −log Λ is y^{1/2}Fφ(y) with y ≳ Λ. The out-of-band tail of a prolate oscillates at frequency ≈ Λ with a slowly varying amplitude. In u this is Mellin frequency ≍ Λ², so the low-pass gains Λ^{−2} per integration by parts.
  - A proof needs derivative bounds on the non-oscillatory amplitude of the tail. That is an archimedean, RH-silent input (Spectrum type). Crude Bernstein bounds do **not** suffice, because they lose Λ per derivative.
- **[S, V(f)] (CC 2006.13771 technique; PROVED, O(1) only).** Moving S through V costs ‖[S, V(f)]‖₁ ≤ C_f, uniformly in Λ. It is not used for the o(1) claim. It is recorded because the naive "Tr(SΔQSV) ≈ Tr(ΔQV)" is false to O(1) without (ii).

**(iii) Where κ would enter (worst case; PROVED).** E-orthonormality gives ‖φ_i‖²_{L²} ≤ κ_Λ. Hence 𝒟(s) ≤ κ_Λ · sup_{unit φ∈R̃}|Mφ(s)|²·(r_K + 2). Also the φ_i are not in the removed span: they carry kept components Q′^{(K)}-feedback. Combined with Spectrum's corrected ω^{arch}, the worst case gives (T3) only if Λ^{B−2}(log Λ)^{O(1)} = o(1), i.e. **B < 2**. This is PROVED given (T3a*), and agrees with Spectrum 65bba74 T4. B is not explicit, so the worst case does not close.

### H.13.3 The profile-weighted route (avoids Λ^B)

**Interface requested from Spectrum: (T3a*), at the corrected rate.** CONJECTURE (Spectrum: HEURISTIC + NUMERICAL, 65bba74 §§1–2). First-principles proof: GAP on their side.

  𝒟_rem(s) := Σ_{n ∈ I_Λ \ I^{(K)}} |Mψ_n(s)|² ≤ C_A r_K (1 + |s|)^{A} Λ^{−2} for |s| ≤ Λ^{ε}, and ∫𝒟_rem = π r_K.

- This is the **trace** (basis-free) form. It is weaker than a sup over the span, and it is exactly what (H.13.a) consumes.
- It matches the model |M cos(c·)|² = π/(2c), c = 2πΛ², and the measured per-state low-Mellin fraction (0.1–0.2)(log Λ)/Λ² at S = log Λ.
- The same bound is asked for k₀ and k₁ (the two constraint-correction vectors), with r_K replaced by O(1). Spectrum's ψ_n(0)² = O(1) is relevant here (NUMERICAL).

**Input (PD), profile domination (GAP; owner Queen; this is where the ζ-weight must be used against the actual profile).** For |s| ≤ Λ^{ε},

  𝒟(s) = Σ_i |Mφ_i(s)|² ≤ (log Λ)^{A} (𝒟_rem(s) + 𝒟_{k₀,k₁}(s)).

In words, the E-orthogonal complement R̃, normalized in the E-norm, has low-s Mellin density no larger than the removed prolates themselves, up to a polylog.

**Theorem H.13 (PROVED-conditional on (T3a*), (PD), (L)).** For fixed f ∈ C_c^∞:

  Tr(S_Λ ΔQ S_Λ V(h)) ≤ C_f (log Λ)^{A} r_K Λ^{−2} ∫|f̂|²|ζ(½+is)|²(1 + |s|)^{A} ds + O(r_K Λ^{−2}(log Λ)^{O(1)}) + o(1) = O(K(log Λ)^{A+2}Λ^{−2}) → 0.

*Proof.*
- (H.13.a) with (i): use (PD), then (T3a*). ∫|f̂|²|ζ|²(1+|s|)^A < ∞ because f̂ is Schwartz and |ζ(½+is)| ≪ (1+|s|)^{1/6}.
- The leak terms: (ii), using the jump part (PROVED) and (L). The factor ‖φ_i‖²_{L²} in the jump sum is absorbed by (PD) as well, via Plancherel ‖φ_i‖² = (1/π)∫|Mφ_i|².
- r_K = O(K log²Λ) by (T1). ∎

So (T3) closes **without B < 2**, given (T3a*), (PD) and (L). Then Q.8.C ⇔ Q.8.C^{(K)}, read as in the lead's audit: given either one, the other ⇔ (T3). M2-compatibility follows.

**Evidence for (PD): what is PROVED, and what is HEURISTIC.**
- *(PD-1) Lower E-norm for flat profiles; PROVED (citation: Balasubramanian–Ramachandra via Prop H.4).* If |Mρ|² is θ-flat at scale H ≥ C log log Λ on its support, then ‖Eρ‖² ≥ cθ log H·‖ρ‖².
  - So E-normalization *shrinks* flat profiles, by 1/log H. It does not inflate them by Λ^B.
  - This replaces m_Λ by a fixed-weight mean, which is exactly Spectrum's point (65bba74 T5): the plunge profile is spread, so B–R and short-interval second moments apply to it, unlike the minimum.
- *(PD-2) Upper mean; PROVED (citation: Hardy–Littlewood mean value, Titchmarsh §7.3, plus Lemma H.2′ style windowing).* For a profile with |Mρ|² ≤ C‖ρ‖²Λ^{−2} on |s| ≤ CΛ² and negligible beyond, ∫|ζ|²|Mρ|² ≤ C log Λ·‖ρ‖².
- *(PD-3) The obstruction: kept feedback; GAP.*
  - φ_i = ρ_i − χ_i, where χ_i ∈ B^{(K)} is the ζ-weighted best approximation of ρ_i by kept states.
  - If |ζ|² were constant, χ_i = 0, since kept ⊥ removed in L². So χ_i is driven only by the **fluctuation** |ζ|² − (mean). In the approximate functional equation, the off-diagonal terms (m/n)^{is}, m ≠ n, are the ℤ-structure of E again, now as dilation-correlations ⟨ρ_i(m·/n), χ⟩.
  - (PD) asks that this feedback not move Mellin mass from the flat, height-≍Λ² bulk into |s| ≤ Λ^ε beyond a polylog factor.
  - HEURISTIC support: |ζ|² is diagonal in s, the removed profiles are already flat in s, and the slowly varying mean part has trace-class commutators with prolate projections (CC 2006.13771 technique). So the feedback should be s-local.
  - What is missing is a **Mellin-local inverse bound for the ζ-weighted Gram on B_{Λ,0}**: a Demko–Moss–Smith-type decay statement for its inverse in the s-variable, with condition-number loss only polylog on the low band. I do not have it. Label: **GAP**.
- **Honest summary.** (T3b) is reduced to one named statement, (PD). The κ-free route requires it. The worst-case route is PROVED but closes only if B < 2.

### H.13.4 Fake-blindness (R5): RH-silent
- **Davenport–Heilbronn.** Every step has a DH analogue with DH's own data.
  - (H.13.a) is exact for any E_a(φ) = x^{1/2}Σa_nφ(nx).
  - (i) uses Mellin–Plancherel with |f_DH|² in place of |ζ|².
  - (ii) uses Poisson with the FE-dual periodic sequence. The sawtooth becomes a mean-zero periodic function (finite Fourier series of period 5 times sawtooths). Same bound.
  - (PD-1): B–R's σ + 2 method needs only a Dirichlet series anchor and convexity; sketch.
  - (PD-2): Montgomery–Vaughan mean value with |a_n| ≤ 1 gives ∫_0^T|f_DH|² ≪ T log T.
  - Off-line zeros of DH play no role anywhere.
- **Epstein ζ_Q.** Same, with 2-dimensional lattice Poisson. The leak's periodic part becomes a lattice-point discrepancy, which is mean-zero but larger. HEURISTIC for the analogue of (ii); the rest is the same.
- Hence Theorem H.13 is fake-blind as a reduction. (T3), and with it M1/M2 compatibility, carries **no RH content**. That is consistent with Q.8.C being "unconditional, not RH-equivalent".

R-line H.13: R1 Λ-uniform, no L | R2 Poisson through the sawtooth (ii); the Euler product only through ζ's mean values (Hardy–Littlewood, B–R) | R3 n/a | R4 pole terms vanish (H.12), and the constraint vectors k₀, k₁ are tracked | R5 fake-blind | R6 the reduction (H.13.a) is exact.

## H.14 Part III summary table

| Claim | Label |
|---|---|
| H.12 Step 1: pointwise Poisson for BV compactly supported φ, symmetric sums, pole terms vanish | PROVED (citation: Dirichlet–Jordan, Zygmund Vol. I Ch. II §8; thm no. not re-verified) |
| H.12: Lemma M1.1 and Mellin–Plancherel/κ on B_{Λ,0}; closure of E = pointwise E | PROVED |
| ⇒ the GAP-minor in M1 is closed; Thm H.9 (M1^{(K)}) is PROVED-conditional on H.1/H.1′ only, i.e. PROVED (sketch) | PROVED (sketch) |
| (T3a) at Λ^{−B−2−δ} and its "edge-localized" premise | **withdrawn**; expected FALSE (Spectrum 65bba74: CONJECTURE, HEURISTIC + NUMERICAL) |
| (T3a*) 𝒟_rem(s) ≤ C r_K(1+|s|)^A Λ^{−2} on |s| ≤ Λ^ε (trace form, corrected rate) | CONJECTURE; requested from Spectrum (RH-silent) |
| H.13.1 ran ΔQ = (1 − Q′^{(K)})E(R₀) = E(R̃); trace formula (H.13.a) | PROVED |
| H.13.2(i) ζ-weight enters as ∫|f̂|²|ζ|²𝒟 | PROVED (exact) |
| H.13.2(ii) jump/sawtooth leak low-pass ≤ C_f|φ(Λ)|²Λ^{−5} (mean-zero sawtooth) | PROVED (elementary) |
| (L) remainder leak low-pass O(r_KΛ^{−2}polylog) | HEURISTIC (needs tail-amplitude derivative bounds; archimedean) |
| ‖[S_Λ, V(f)]‖₁ ≤ C_f (CC 2006.13771 technique) | PROVED (O(1) bookkeeping only) |
| Worst-case transport: (T3) ⇐ (T3a*) if B < 2 | PROVED-conditional (agrees with Spectrum T4) |
| (PD) profile domination of the E-orthonormal complement | **GAP (Queen)**: needs a Mellin-local inverse bound for the ζ-weighted Gram |
| (PD-1) flat profiles: ‖Eρ‖² ≥ cθ log H‖ρ‖² | PROVED (citation: B–R via Prop H.4) |
| (PD-2) flat profiles: ‖Eρ‖² ≤ C log Λ‖ρ‖² | PROVED (citation: Hardy–Littlewood, Titchmarsh §7.3) |
| Thm H.13: (T3a*) + (PD) + (L) ⇒ (T3), rate O(K(log Λ)^{A+2}Λ^{−2}), no B < 2 | PROVED-conditional |
| H.13 fake-blind (DH, Epstein) ⇒ RH-silent | PROVED (sketch) as a reduction |
| Prior art (CC 2106.01715 §3, CC 2006.13771 §1, CCM 2310.18423) contains no (T3b) | checked in local copies |

Not cited: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.
