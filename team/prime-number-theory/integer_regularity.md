# Integer regularity versus off-line zeros in Beurling systems

Goddess of Prime Number Theory, 2026-10-22. **RH is OPEN.**
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

Labels: PROVED (literature, cited), CERTIFIED (Arb/interval), NUMERICAL, CONJECTURE.

**Question (lead, item 4).** Can a Beurling system with integer-count error N(x) − ρx = O(x^θ), θ < 1/2, have zeros
off the critical line?

**Short answer.** Yes. **PROVED, unconditionally**, but non-constructively (§2). The live questions are quantitative:
how small can θ be for a given off-line zero, and can an example be built without assuming RH (§3–4).

## 1. Notation and the dictionary
An **[α,β]-system** (Hilberdink 2005) is a Beurling system with ψ_P(x) = x + O(x^{α+ε}) and N_P(x) = ax + O(x^{β+ε}) for
every ε > 0, and for no ε < 0.

**Lemma 1 (PROVED, standard).** Let β < 1/2. Then ζ_P(s) − as/(s−1) = s∫ x^{−s−1}(N_P(x) − ax)dx is analytic on Re s > β.
So in Re s > 1/2 the only singularities of ζ_P'/ζ_P are the pole at 1 and the zeros of ζ_P. A zero ρ with
Re ρ > 1/2 forces ψ_P(x) − x = Ω(x^{Re ρ − ε}) (Landau/Ingham). Conversely, α > 1/2 forces such a zero (Landau's
theorem on Dirichlet integrals of one-signed measures, applied to ψ_P(x) − x ± x^{α'}).
Hence, **for β < 1/2: "ζ_P has an off-line zero" ⇔ α > 1/2.** The question becomes:
*do [α,β]-systems with α > 1/2 > β exist?*

## 2. What is known
| statement | label | source |
|---|---|---|
| max(α, β) ≥ 1/2 for every [α,β]-system | PROVED | Hilberdink, J. Number Theory 112 (2005) 332–344 |
| Of α, β, γ (γ the Möbius-sum exponent), the two largest are equal and ≥ 1/2 | PROVED | Hilberdink, arXiv:1901.06866 |
| Systems with α, β ≤ 1/2 exist (probabilistic) | PROVED | Zhang, Math. Ann. 337 (2007) 671–704 |
| [0, 1/2]-systems exist | PROVED | Broucke–Vindas, arXiv:2102.08478 (Math. Z.) |
| [α,β]-systems exist for every α ∈ [0,1), β ∈ [1/2,1) | PROVED | Broucke–Debruyne–Révész (BDR), arXiv:2309.01567, Thm 1.1 |
| **Under RH:** [1/2, β]-systems for 0 ≤ β < 1/2, and [α,β]-systems for 1/2 < α < 2/3, 2α/(α+2) ≤ β < 1/2 | PROVED (conditional on RH) | BDR Thm 1.3 (§5) |
| (ℙ, ℕ) is a [1/2, 0]-system ⇔ RH | PROVED (classical) | e.g. Montgomery–Vaughan Ch. 13 |
| Off-line zeros with large prime oscillation, β > 1/2 | PROVED | DMV, Math. Ann. 334 (2006); Révész IMRN 2023 |
| RH with large integer oscillation | PROVED | Broucke–Debruyne–Vindas, Adv. Math. 370 (2020) 107240 |

**BDR's §5 system**, for reference: P_{α,β} = ℙ ∪ ℙ^{1/β} \ P_S. Here P_S ⊂ ℙ is a random-approximation subsequence with
π_S(x) ~ Li(x^α), and log ζ_S(s) = log ζ(s+1−α) + O(√log|t|) on Re s ≥ α/2 + ε.
So ζ_{α,β}(s) = ζ(s) ζ(s/β) / ζ_S(s) **has a real zero at s = α ∈ (1/2, 2/3)**, coming from the pole of ζ_S. Its integer
error exponent is β, obtained by Dirichlet's hyperbola method (BDR Lemma 5.1).

**Theorem A (PROVED, unconditional, non-constructive).** There is a discrete Beurling system with N(x) = ax + O(x^θ) for some
θ < 1/2 whose zeta function has a zero with real part > 1/2.

*Proof.* Two cases.
- If RH fails, take (ℙ, ℕ): N(x) = ⌊x⌋ = x + O(1), so θ = 0, and it has an off-line zero.
- If RH holds, BDR Thm 1.3 gives a [3/5, β]-system for any β ∈ [6/13, 1/2), e.g. θ = 0.47, with the zero at s = 3/5. ∎

The same dichotomy gives off-line zeros with θ arbitrarily close to 2α/(α+2) ∈ (2/5, 1/2). Combined with
finite_windows_dmv.md §3 (gluing to the true primes below X, which changes neither θ nor the zeros in Re s > 0), we get
**Corollary B (PROVED, non-constructive).** For every X there is a system with θ < 1/2, an off-line zero, and Λ_P = Λ on [1, X]. It is
therefore indistinguishable from ζ by every window form with L ≤ (log X)/2.

## 3. What is open
1. **A constructive or unconditional example.** No α > 1/2 > β system is known without RH. BDR say explicitly that discretising
   a continuous template costs ≥ x^{1/2} on the integers, so β < 1/2 needs a discrete construction "from the start".
   Their construction starts from ℙ. **OPEN.**
2. **The trade-off curve.** β*(α) := inf{β : an [α,β]-system exists}, for α ∈ (1/2, 1).
   - Known: β*(α) ≤ 2α/(α+2) for α < 2/3, under RH.
   - For α ≥ 2/3 nothing below 1/2 is known. BDR Remark 5.3: their double-prime variant reaches α < 4/5, but only with β ≥ 1/2.
   - **Lower bounds on β*(α) are RH-hard.** (ℙ, ℕ) has β = 0. So any unconditional theorem "β < c ⇒ α ≤ 1/2", for any c > 0, implies RH,
     and by Theorem A it is false for c > 2α/(α+2) whenever α ∈ (1/2, 2/3).
   **Honest note:** "integer regularity forces the zeros onto the line" is therefore a *strengthening* of RH, not a route that is easier than RH.
   It can only become useful with extra structure that ℕ has and the BDR systems lack, see item 4.
3. **Is 2α/(α+2) a method artifact?** BDR Remark 5.3(2): "certain error terms … can be improved with more advanced technology
   than Dirichlet's hyperbola method". **OPEN.**
4. **Which structure is missing from the fakes.** The BDR/DMV fakes have an Euler product, Λ ≥ 0, PNT, and (BDR) β < 1/2.
   What they lack relative to ℕ is: (i) integers that are rational integers with N(x) = ⌊x⌋ exactly (additive structure, β = 0 with
   bounded error); (ii) a functional equation. Is there a theorem of the form "N_P(x) = ax + O(1) (or O(x^ε)) ⇒ …"? I found **no**
   literature on Beurling systems with *bounded* integer error other than ℕ itself. This is a well-posed, apparently unstudied question.

## 4. First concrete line of attack: push β*(α) down by beating the hyperbola method

*Where 2α/(α+2) comes from (my reconstruction, consistent with BDR §5).* The integers of ℙ \ P_S are the n coprime to P_S. Their count is
 C(x) = Σ_{d ∈ N_S} μ_S(d) ⌊x/d⌋ = x/ζ_S(1) − Σ_{d ≤ x} μ_S(d) ψ̃(x/d) + O(x^{α/2+ε}),
where ψ̃(u) = {u} − 1/2 is the sawtooth, and the tail uses M_S(x) ≪ x^{α/2+ε} (BDR §5).
- Splitting at y and bounding the d ≤ y sawtooth sum **trivially** by N_S(y) ≈ y^α,
- and the d > y part via M_S by x^{α/2}(x/y)^{1−α/2},
- then balancing gives y = x^{2/(2+α)} and error x^{2α/(α+2)}.

That is exactly BDR's exponent. So the loss is the trivial bound on
 S(x, y) = Σ_{d ≤ y, d ∈ N_S} μ_S(d) ψ̃(x/d).

**CONJECTURE G1.** For the random P_S (independent selection with P(p ∈ P_S) ≍ p^{α−1}), almost surely
S(x, y) ≪ y^{α/2+ε} + (x-dependent lower-order terms). The heuristic is square-root cancellation of signed sawtooth values over a set of
size y^α. Then the balance gives **β*(α) ≤ α/2 + ε** (under RH), e.g. β ≈ 0.30 at α = 0.6, against BDR's 0.4615.
- **Why it is plausible.** It is a sum of a random multiplicative function (μ_S is a random ±1/0 multiplicative function) against a
  deterministic bounded weight. Moment methods for random multiplicative functions (Harper; Granville–Soundararajan-type
  arguments) give square-root (or slightly better) cancellation for Σ f(n) a(n) with fixed bounded a(n). The new part is the
  dependence on x through ψ̃(x/d). This needs a union bound over x, which is standard via a net in x of size x^{O(1)}.
- **Which fakes it separates.** None. It *builds* better fakes. It sharpens the statement "integer regularity
  cannot detect off-line zeros" from θ ≥ 2α/(α+2) to θ ≥ α/2, i.e. down to θ ≈ 1/4 for zeros just right of 1/2.
  It is not a route to RH, and it does not restate RH. It maps how much of ℕ's regularity can be faked.
- **Cheap numerical test (planned, not yet run).** Sieve to X = 10^9 with α = 0.6, generate P_S, and measure max_{x≤X} |C(x) − x/ζ_S(1)| on dyadic
  blocks. The predicted exponents are 0.30 (G1) and 0.4615 (BDR bound), and X^{0.30} ≈ 500 vs X^{0.46} ≈ 1.4·10^4 are separable.
  Caveat: the density 1/ζ_S(1) depends on deleted primes > X. A truncated system has the extra term x·Σ_{d>x} μ_S(d)/d, whose
  cancellation needs the full M_S. So the measurement should stay in x ≤ X^{1/2}, or include a synthetic continuation of P_S.
  We did not run it in this session (CPU was committed to the precision scan). Logged as a next step.

**Second line (for Daniel's "where do integers enter").** Look for the weakest property P of ℕ, short of "N(x) = ⌊x⌋",
that BDR-type fakes provably lack, and that a Weil-window argument could use *uniformly in L*. Candidates: (a) N(x) − ax = O(1);
(b) the integers form an additive semigroup (n, m ∈ N ⇒ n+1 ∈ N …); (c) Poisson summation for the integers, which is what yields
the functional equation. Only (c) is known to interact with the archimedean term W_∞ of the Weil form. This is the Connes–Consani
mechanism: the prolate/Sonin span comes from Poisson summation for ℤ, so a fake without Poisson summation cannot reproduce it.
**CONJECTURE-level; no result yet.**

## 5. Fakes battery relevance
| fake | Euler product | Λ ≥ 0 | PNT | integer regularity | functional equation | relevant result |
|---|---|---|---|---|---|---|
| DMV (glued, P_X) | yes | yes | yes (dlVP) | **only θ > 1/2** | **no** | finite_windows_dmv.md |
| BDR §5 (under RH) | yes | yes | yes, α ∈ (1/2, 2/3) | **θ < 1/2** | **no** | Theorem A |
| DH | **no** | no | — | — | yes | separated by windows |
| F_θ | yes (non-tempered) | yes | — | — | yes | separated by windows |
| H₋ₜ | no (zero-side fake) | — | — | — | yes | separated by windows |

**DMV row, precisely.** DMV has an Euler product, Λ_P ≥ 0, PNT with the de la Vallée Poussin error, and N(x) = κx + O(x^θ) with
θ ∈ (1/2, 1). It **lacks** (i) integer regularity below x^{1/2} and (ii) additive/Poisson structure, hence any functional equation.
BDR shows that (i) alone does not stop off-line zeros (Theorem A). So among the standard properties, only (ii) separates.

## References
- T. W. Hilberdink, *Well-behaved Beurling primes and integers*, J. Number Theory 112 (2005) 332–344. https://doi.org/10.1016/j.jnt.2004.09.001
- T. W. Hilberdink, *The average order of the Möbius function for Beurling primes*, arXiv:1901.06866. https://arxiv.org/abs/1901.06866
- F. Broucke, G. Debruyne, Sz. Gy. Révész, *Some examples of well-behaved Beurling number systems*, arXiv:2309.01567. https://arxiv.org/abs/2309.01567
- F. Broucke, J. Vindas, *A new generalized prime random approximation procedure and some of its applications*, arXiv:2102.08478.
- F. Broucke, G. Debruyne, J. Vindas, *Beurling integers with RH and large oscillation*, Adv. Math. 370 (2020) 107240.
- W.-B. Zhang, *Beurling primes with RH and Beurling primes with large oscillation*, Math. Ann. 337 (2007) 671–704.
- H. G. Diamond, H. L. Montgomery, U. M. A. Vorhauer, Math. Ann. 334 (2006) 1–36. https://doi.org/10.1007/s00208-005-0638-2
- Sz. Gy. Révész, *Oscillation of the remainder term in the PNT of Beurling, "caused by a given ζ-zero"*, IMRN 2023(14) 11752–11790.
- H. G. Diamond, W.-B. Zhang, *Beurling Generalized Numbers*, AMS Surveys 213 (2016).
- A. Connes, C. Consani, *Spectral triples and ζ-cycles*, arXiv:2106.01715.
