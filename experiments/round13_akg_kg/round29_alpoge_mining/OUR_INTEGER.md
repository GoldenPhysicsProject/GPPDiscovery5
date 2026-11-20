# Round 29 (C): "Where is our integer?"

**Brief:** Daniel t525u ("All three"), item C. Daniel's standing rule: Weil positivity is an **output**, never the target.
**Question.** Alpöge's sign (NOTES.md §3) is: *a nonzero integer has |N| ≥ 1*. Primes are forced to divide N with known
multiplicity (ASK), N's size is bounded (MAIN), and therefore ASK − MAIN ≤ error. With a multiplicity amplifier, the error is
√x-sized, and an off-line zero would make ASK too big. **What integer could play N for our D(L)?**
Labels: PROVED (classical, cited), CHECKED (computed here), HEURISTIC, OPEN, DEAD END.

---

## 1. The picture in Daniel's language (HEURISTIC)

- **Positive integers counted once.** Every n ≥ 1 is built exactly once from primes (unique factorisation). In log form:
  log n = Σ_{d|n} Λ(d). Each prime power p^k is a single slot that n either fills or doesn't. Primes are the fermions: a slot
  holds at most what its occupation allows, so you can't put "a bit more p" into n.
- **The two faces x and 1/x are the product formula.** For a positive integer n, the archimedean face is |n|_∞ = n. The prime
  faces multiply to Π_p |n|_p = 1/n. So **the 1/x face of an integer is its prime content**, and n · (1/n) = 1 is the
  conservation law: Σ_{all places v} log|n|_v = 0. Every prime face is one-signed, |n|_p ≤ 1 ("no negatives"). So the
  archimedean face has to carry the balance. That is the whole of "|N| ≥ 1".
- **The balanced ½ weight.** The n^{−1/2} in our S(L) = Σ Λ(n) n^{−1/2}(L − log n) is (Π_p|n|_p)^{1/2}: half of the prime face.
  At s = ½ the two faces are weighted equally. That is where the zeros sit (s = Δ = ½ + iγ: massless boundary modes, ++ / −−
  balanced).
- **What an off-line zero would mean.** A zero at ½ + δ + iγ adds an oscillation of size x^{½+δ} to the prime sums. In
  integer language it asks some window of primes to fill **more slots than the archimedean volume allows** (ASK > MAIN by
  x^{½+δ}), or fewer. Only the "more" side is forbidden by |N| ≥ 1. That is why both Alpöge's bound and F13 are **one-sided**.
- **Zitter.** The Weyl halves of a massive fermion talk to each other. For a prime, the "talk" between its two faces is
  |p|_∞ · |p|_p = 1, and it always happens. The zeros are the modes for which the two faces are exactly balanced (½, ½): no
  Higgs, no mass, the boundary. **This is a restatement, not a mechanism.** The mechanism must be an integer.

## 2. What the integer has to do (the spec)
To give a √x-quality one-sided bound on our prime sums, an integral object N = N(x, m) must have:
1. **N ≠ 0**, with |N| ≥ 1 or an integral-lattice analogue (covolume, index, positive-definite integer Gram determinant).
2. **Forced divisibility:** v_p(N) ≥ (weight of p in the D(L) window) × m for every prime p ≤ x, with a known multiplicity m.
3. **Size:** log|N|_∞ ≤ m·(main term of the window) + (trivial loss of order x, *not multiplied by m*) + (cost of order m²).
4. Then m·(prime-window sum − main) ≤ x + m², which gives an error of x/m + m, about √x at m ≈ √x (TOY_FQ.md, NOTES.md §4).
The hard requirements are (2) with multiplicity m and the "trivial loss not multiplied by m" in (3).

## 3. Candidate survey (blunt)

| # | Object N | Forced divisibility | Size bound | What it gives | Verdict |
|---|---|---|---|---|---|
| 1 | **lcm(1..N) = e^{ψ(N)}**: ψ is *exactly* the Λ-sum in D(L) | v_p = ⌊log N / log p⌋, exact | none independent; it *is* e^{ψ} | alone, a tautology | **DEAD END alone** |
| 1′ | Nair/Gelfond: lcm(1..2n+1)·∫₀¹ P(x)dx ∈ ℤ_{>0} for P ∈ ℤ[x] of degree 2n, e.g. x^n(1−x)^n | the lcm clears all denominators | ∫ P small (integer Chebyshev problem) | ψ(x) ≥ c·x (Chebyshev-type) | **DEAD END for PNT, let alone RH.** PROVED (classical): the method gives ψ(x) ≥ x·log(1/t_ℤ), where t_ℤ([0,1]) is the integer transfinite diameter. Known lower bounds t_ℤ > 0.42 > 1/e make c = 1 unreachable (Gorshkov, Montgomery's *Ten Lectures* ch. 10, Borwein–Pinner–Pritsker; constants as recalled, re-verify before citing). |
| 2 | **Factorials x! and Chebyshev ratios Π (x/n_i)!^{a_i}** (Landau's criterion: integer iff Σ a_i⌊t/n_i⌋ ≥ 0) | Legendre: v_p(x!) = Σ_k ⌊x/p^k⌋, the *same shape as Alpöge's ASK* | Stirling | log x! = Σ Λ(n)⌊x/n⌋ (Chebyshev identity); ratios give θ(2n) − θ(n) ≤ 2n log 2 | **The right shape at multiplicity 1, and stuck there.** CHECKED (`chebyshev_slack.py`): for N = C(2n,n) the slack MAIN − ASK ≈ 0.39·n for n = 10³…10⁶, never √n. Chebyshev-type inequalities can push constants towards 1 only with PNT as input (Diamond–Erdős). No multiplicity amplifier. |
| 3 | Hankel/Toeplitz determinants of Λ-sums, Li/Keiper coefficients, our box-Gram (R23 Thm 7) | none: these are not integers | positivity of the matrix | RH ⇔ PSD | **DEAD END: this is Weil positivity in disguise**, the target, not the input. Excluded by the standing rule. |
| 4 | **Redheffer determinant** det R_n = M(n) (Mertens) | none forced | Hadamard (trivial) | RH ⇔ M(n) ≪ n^{½+ε} | **Right integer, wrong direction.** Integrality only gives lower bounds on \|N\|, and M(n) can be 0. RH needs an *upper* bound on \|M\|. |
| 5 | Gram determinant / discriminant of a theta-series lattice; Epstein zeta | disc ∈ ℤ_{>0} | covolume | lattice point counts | **DEAD END alone:** Epstein zetas have off-line zeros (our EPS fakes), so lattice integrality without an Euler product (Hecke structure) cannot force the line. With CM and Hecke structure it becomes Alpöge's route (#8). |
| 6 | Bernoulli numerators, K-group orders, class numbers, Ш | von Staudt–Clausen, Kummer congruences; Lichtenbaum-type formulas for \|K_{4k−2}(ℤ)\|/\|K_{4k−1}(ℤ)\| | exact functional-equation size | values of ζ at single points | **OPEN, secondary.** Alpöge uses Ш-length bounds (`VWin`, the "Vandiver weakening") as an auxiliary input to kill twisted classes, not as the main N. Point values carry no window information. |
| 7 | **Auxiliary integer polynomial / resultant with simultaneous high-order p-divisibility for all p ≤ x** (Thue–Siegel–Bombieri style, "integer Stepanov lemma") | multiplicity m at every prime | product formula | would give √x | **THE target. OPEN.** This is R16's named gap ("GAP = integer Stepanov lemma"), and it is what Alpöge's CM-theta derivatives supply. |
| 8 | Gram determinant under the **Fock/Bombieri pairing** ⟨f,g⟩ = Σ f_μ g_μ μ! of integer polynomials (Alpöge's `Formalism.A.GramZ` template) | det divisible by Π_a α_a!, i.e. multiplicity as factorials | Hadamard: det ≤ Π diagonal | positive-definite integer Gram matrix ⇒ det ≥ 1 | **Most promising elementary toy. OPEN.** It gives *multiplicity as factorial divisibility* without derivatives on ℤ. It is the one place "x² ≥ 0" (PosDef) and "counted once" (integer det ≥ 1) appear together. |
| 9 | CM theta derivatives D^a θ(τ_K) / Ω^a (Alpöge's actual N) | D^p ≡ 0 mod an inert p (supersingular), so p^{a/p + a/p² + …} | Chowla–Selberg plus Stirling | his ASK − Main ≤ error, then m lifts | Works (if his tree is right) but **is not physics**. It is the arithmetic of CM elliptic curves. |

### Where ψ meets this exactly (CHECKED identities, HEURISTIC reading)
- e^{ψ(x)} = lcm(1..⌊x⌋) and Σ_{n≤x} Λ(n)⌊x/n⌋ = log⌊x⌋! (Chebyshev). Alpöge's Gz uses the Riesz weight (x/n − 1), which is
  ⌊x/n⌋ minus a sawtooth. So **his Gz is "log x! minus Stirling" with the floors smoothed out.** Our S(L) carries n^{−1/2}, and
  **that half-shift destroys integrality**. Π n^{Λ(n) n^{−1/2}} is not an integer. So D(L) itself is not the log of an integer.
  The integer lives at the unshifted level (ψ, x!, lcm), and the ½ appears only after the √x amplifier. That fits Daniel's
  picture: the ½ is the *balance point*, not part of the counting.

## 4. Verdict (blunt)
- **No integer attached to D(L) in an obvious way gives more than multiplicity one.** lcm, factorials and Chebyshev ratios
  are exactly the right objects (ψ is in D(L), and Legendre's formula is Alpöge's ASK). They stall at error of order x
  (CHECKED: slack ≈ 0.39n), and the integer-Chebyshev obstruction (PROVED, classical) shows that refining the integer at
  multiplicity one cannot even reach PNT.
- Positivity-type determinants (#3) are the forbidden target. Lattice integrality without Euler products (#5) is refuted by
  our own Epstein fakes. Redheffer (#4) points the wrong way.
- **The missing piece is a multiplicity amplifier over ℤ:** an integral object depending on m, divisible by p^{m·(…)} for all
  small p, with a trivial size loss that is *shared* across the m copies. Alpöge found one in CM theta arithmetic.
- **Physics reading (HEURISTIC):** we need a *quantised degeneracy*. That means a count of states (an integer) in which each
  prime-fermion mode is forced to appear m times by a symmetry (multiplicity), inside a phase-space volume (size) that does not
  grow with m. Pauli exclusion, "counted once", is the |N| ≥ 1. The "m copies with shared vacuum" is the amplifier.

## 5. Proposed next steps (none aimed at W ≥ 0)
- **R30-a (toy, cheap):** build #8. Take Bombieri-pairing Gram matrices of x^α·E_a(x) families. Measure how much prime-window
  divisibility (primes ≤ x) is forced against the Hadamard size, as a function of the multiplicity |α|, and see whether any
  family gets slack below the 0.39n of C(2n,n). **OPEN.**
- **R30-b (physics, Daniel):** name a quantum system whose ground-state degeneracy is an integer that the primes ≤ x must
  divide with multiplicity m. Candidates: a free-fermion system with one mode per prime power (the primon gas / Julia gas,
  where the partition function is ζ), and its m-fold tensor power with a shared vacuum. Then ask what bounds the degeneracy's
  size. **OPEN, HEURISTIC.**
- Report to the Lean worker: none of this changes the F13-window target.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
