# Round 29: mining Alpöge's RH/GRH Lean release for ideas (step 1 of 3: steal)

**Date:** 2026-11-20. **Brief:** Daniel, t521u–t522u: "I want to steal insight and understanding from his work but make our own proof."
**Standing rule (Daniel):** Weil positivity is an OUTPUT, never the target. Nothing below proposes a round aimed at W ≥ 0.
**Credit:** The mathematics mined here is **Levent Alpöge's** (with Claude, per his note; Lean release of 2026-11-20, `remark.pdf`).
Turbo filed the release in GPP-bridge `research/external/alpoge-rh/`. **Claude Code (D1)** wrote the first understanding notes
(GPP-bridge `CLAUDE_CODE_RESEARCH_NOTES.md`, section "Alpöge RH/GRH release: D1's understanding notes") and made the first
kernel receipt (`Formalism.landauLemma_of`, standard three axioms). This round starts from those notes and goes deeper. Where it
disagrees with them, it says so (§1.3).

**We do not claim RH. We do not claim Alpöge's proof is correct.** Nobody here has built his 175k-file tree. Labels:
* **PROVED [A]**: a `theorem` with a full proof term in Alpöge's source (no `sorry`, per Turbo's scan). We **read** it but did not
  kernel-check it. Its trust is the trust of his unbuilt tree, nothing more.
* **READ**: a definition or statement read in source. **INFERRED**: my reading of what the code means.
* **HEURISTIC**: translation into Daniel's language. **OPEN**: not done.

Archive: the 14 parts reassemble to sha256 `9631b4f0…c7fc` (matches PROVENANCE.md). I extracted them to a scratch dir and indexed
1,050,822 declarations. Import reachability from `Main.lean` covers 175,149 of 175,187 files, so being imported says nothing about
being *used*. The trace below therefore follows actual proof terms, by hand, from `everyCharacter` down.

---

## 0. One-paragraph answer

The engine is **not** a Weil-form positivity argument. It is a **Liouville / product-formula inequality**: a nonzero algebraic
quantity cannot be divisible by more prime powers than its archimedean size allows. The prime divisibility comes from
**derivatives of theta values at CM points**. At a prime p inert in the CM field, the derivation D satisfies D^p ≡ 0 (mod p),
which gives Legendre-type divisibility p^{a/p + a/p² + …}. This is the integer substitute for Stepanov's "vanishing to high
multiplicity". A **chain of m theta lifts** repeats the same divisibility m times. The archimedean "trivial" error does **not**
grow with m. Each lift costs about one unit. So the error per unit of signal is x/m + m, which is smallest at **m ≈ √x**. That is
the "√X chain length". Averaging over a character family (finite Fourier analysis / Parseval on a finite abelian group) turns
"inert primes can't be over-represented" into an **upper bound on the Riesz-weighted von Mangoldt sum**,
G_z(x) ≤ x^{1/2+Cε}, on long windows. A windowed Landau–Pringsheim lemma then gives the zero-free half-plane. Positivity of
anything Weil-like is an output.

---

## 1. The spine, top to bottom (all file paths relative to `lean/src/`)

### 1.1 Chain of proof terms (READ; each line is a PROVED [A] theorem unless marked)

1. `UpperTheory.FinalTheorems.riemannHypothesis := riemannHypothesis_of_everyCharacter everyCharacter`
   (`Development/AutomorphicArithmetic/CharacterBoxes/MainConclusion.lean`).
2. `everyCharacter := everyCharacter_of_row galoisLocalGlobalStrong_recParameter` (same file).
   `everyCharacter_of_row h := LiteratureTheorems.every_character_of_literature_of_vWin ⟨h, recTExists, h⟩ (One.vwin_of_row h)`.
   The structure `OpenLiteratureStatements` (in `…/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean:215`) has three fields:
   `caraianiLocalGlobalDual`, `localLanglandsExists` and `galoisRepOfPolarizedPairLit`. All three are **filled by his own theorems**,
   not assumed. The field names say which literature results he re-proved.
3. `every_character_of_literature_of_vWin` → `Lit.every_character_of_theta_integrals_local_matching`
   (`…/ThirdPlaceSeries/SiegelEisensteinSecondTermHypotheses.lean:69`) → `Top.every_character_of_green_v_win` →
   `Top.zero_free_all_of_green_v_win` (`…/ThirdPlaceSeries/TopPlug.lean:12`), which is
   `zeroFreeAll_of_topG (la := 3/2) … (fun ψ => supplyChi ψ (hV ψ))`.
4. `UpperTheory.Top.zeroFreeAll_of_topG` (`…/CharacterBoxes/ZeroFreeAllOfTopG.lean`) applies
   `Formalism.HC.zeroFreeOn_of_tiedWin_cyc` with three inputs:
   * `UpperTheory.Comparison₂.level_gm_le_on_family_scales`, the **analytic inequality from family data** (§3);
   * `UpperTheory.Comparison₂.formB`, the **windowed Landau lemma** (§2);
   * `TopG`, which supplies, for every ε, a family `Fm` satisfying `ChainTied` (the arithmetic/automorphic bulk: theta lifts,
     slopes, divisibility). TopG comes from `shapeWin` and `Green`. **I did not open TopG's proof. This is where the 140k-file bulk lives.**
5. `VWin` (`…/CharacterBoxes/ExistsPageN.lean:49`, READ) is a bound on **lengths of Tate–Shafarevich-type groups**
   `shaLen (H¹ classes, Q_ℓ/Z_ℓ(−e))` over cyclotomic fields:
   Σ_{ℓ∈S} log ℓ · Σ_{1≤e<t} len Ш ≤ √exp(x^ε). It is used to "kill classes by twisting" (remark steps 14 and 17, the "weakening
   of Vandiver"). It feeds `supplyChi : VWinAt ψ 401 → HostExW (3/2) (GChi ψ).req`.

### 1.2 Where the hard content sits (INFERRED from the names and from what is consumed)
* Arithmetic geometry and automorphic forms: `ChainTied` and `TopG`.
* Analytic number theory: `ch_top`, `ch_core`, `core_gm`, `formB` (§2–§3). That part is short and readable.

### 1.3 Correction to D1's notes (INFERRED, checked by tracing)
D1 guessed the outer layer is `Formalism.LandauStep` / `RHofG` (`Formalism/Statement.lean`), with a **lower** bound
−C x^{σ₀} ≤ G_χ(x). **Those Props exist but are not on the main chain.** No proof term on the chain mentions them; the only file
that mentions them is `Statement.lean` itself. The chain uses an **upper** bound on a weighted family sum `Λ.Gm`, via
`ChebyshevInequality*` and `FormB` (§2). D1's §1.4 (the Landau architecture matches our F12/F13) is right in spirit. The sign
convention and the exact objects differ, as set out below.

---

## 2. The outer layer: exact one-sided bound and the Landau step (task 3)

**Objects (READ).**
* `Formalism.Gz x := Σ_{n≤x} Λ(n)(x/n − 1) − x log x + (1+γ) x` (`Formalism/StatementOC.lean:13`).
  This is the Riesz mean of ψ minus its main term. Its Mellin transform is −ζ′/ζ(s)/(s(s−1)) plus explicit poles.
* `An.gChi χ y := Σ_{n≤y} Λ(n) χ(n)(y/n − 1)_+ + y · L′/L(1,χ)` for χ ≠ 1
  (`Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/Analytic.lean:23`).
* `Level.Gm y := Gz y + Σ_{ψ ∈ Nt} wt ψ · Re gChi(ψ, y)`
  (`…/TotientGaloisClasses/ExplicitFormulaErrorOfHostConditions.lean`, inside `namespace Level`). Here
  `Level.cyc ψ = Level.ofSub (zpowers ψ.χ)`. KEPT is the cyclic group generated by ψ, Nt is KEPT minus the trivial character,
  `psK = none`, and wt = 1/(3o+1) with o = |⟨ψ⟩|.
  **For RH itself (ψ = 1): Nt = ∅, so Gm = Gz exactly** (READ: `Level.ofSub`, `Level.Nt`).

**The input to Landau (READ, statement of a PROVED [A] chain).**
* `Formalism.HC.Level.FormB`: for every la > 1 and θ ∈ (0,1), if `HasWindows {x | Λ.Gm x ≤ x^θ} la`, then every L(ψ,·) with
  ψ ∈ KEPT is nonzero on Re s > θ (`…/TotientGaloisClasses/Window.lean:30`).
  `HasWindows S la` means: for every y there is Y ≥ y such that **all** x ∈ [Y, Y^la] lie in S.
* It is proved as `UpperTheory.Comparison₂.formB` (`Development/AutomorphicArithmetic/PrimitiveCosets/Imports.lean:8`) via
  `Deduction.zeroFree_of_rpow_window` (`…/DedekindLFunctions/LandauTheorems/ZeroFreeOfRpowWindow.lean:11`):
  if ∀Y₀ ∃Y ≥ Y₀ with Gm(y) ≤ C y^θ on [Y, Y^{1+δ}], then ZeroFree above θ. The proof (`lemmaL` → `landau_window`, with files
  `LandauKernelMax`, `LandauNearest`, `IteratedDerivLaplace`) is a **windowed** Landau argument, not the plain
  abscissa-of-convergence one. I read the statements only.
* The bound actually delivered, `ChebyshevInequalityOnWindows` (`…/TotientGaloisClasses/TopWin.lean`), is
  **Λ.Gm(x) ≤ x^{1/2 + (eD+eL+7)ε} for x ≥ x₁ in the family's window set**. It holds for every small ε > 0, with la = 3/2.
* A non-windowed sibling, `Formalism.HC.Landau.landau_family` (`…/TotientGaloisClasses/Landau.lean:380`, PROVED [A]), reads:
  R·Gz(y) + Σ_{ξ∈S} Gchi_ξ(y) ≤ C₀ y^θ for y ≥ y₀ (S closed under inverse) ⇒ ζ and every L(ξ) are nonzero on Re s > θ.
  Its engine `landau_entire` uses H = (non-negative) with Mellin transform (P′/P − (a+bs))/(s(s−1)) + C/(s−θ).

**Why an *upper* bound is the right direction (INFERRED).** H = C x^θ − Gm ≥ 0. Its Mellin transform carries **+P′/P**. At a
real zero β < 1 that term tends to −∞ as σ ↓ β, which contradicts H ≥ 0. So real zeros (including Siegel-type ones) are excluded
by the sign alone, with no separate hypothesis. Complex zeros are excluded by the usual Landau comparison
|F(σ+iγ)| ≤ F(σ).

**Dictionary to our F13 (INFERRED).** R23: D(L) = A(L) − S(L) with S(L) = Σ_{n≤e^L} Λ(n) n^{−1/2}(L − log n). So
"D(L) ≥ −C e^{εL}" is "S(L) ≤ A(L) + C e^{εL}". That is **the same direction** as Gz(x) ≤ x^{1/2+ε} with x = e^L: an upper bound
on a smoothed prime sum. The kernels differ: his is (x/n − 1)_+ with Mellin 1/(s(s−1)), ours is n^{−1/2}(L − log n)_+ with
Laplace 1/σ² at the ½-shift. One-sided bounds do not transfer between kernels pointwise for free. **Steal #1:** his `FormB`
shows the bound is only needed **on infinitely many long windows [Y, Y^{3/2}]**, not for all L. That is a strictly weaker
hypothesis than PR #649's pointwise F13.

---

## 3. The sign: what actually forces it (task 1, steps 4 and 6)

### 3.1 The inequality the bulk delivers (READ, PROVED [A])
`Deduction.ch_top` (`…/LandauTheorems/ChTop.lean:37`): for each node u,

  **gapA(u) = ASK(u) − Main(u) ≤ x^{CErr·ε} · errorBracketSizeWindows(x, ε, n, n_T, n_T2, m, b, k, a_E).**

It is derived from `SizeData₂` (the hypotheses `HidaSetting, DIV, EComp, HComp, Count, Slopes, QuotientSlopeBoundsWith,
LossBounds, SlopeIneq, HHBound, FParams`, defined in `Formalism/RankGapBounds/MeanBoundSecondForm.lean`), together with
`ErrWindows` and `BNBounds`.

* `ASK` (`Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean:84`, READ) is
  Σ_{p ≤ N inert} 2 log p · Σ_{i ≤ N, i·Frob_p ∉ Φ} max(a/p^i − w(partner), 0), plus a ramified term at d.
  Every summand is ≥ 0 (a `max(·,0)`). Its shape, Σ_i (a/p^i − w)_+, is **Legendre's formula** v_p(a!) = Σ_i ⌊a/p^i⌋, shifted
  by w. So ASK is "how much p-adic divisibility the inert primes force".
* `Main` = Σ_{j<m} Σ_t a_t (log a_t − 1 − γ + L′/L(1,χ_C) − mean_ρ log w): a Stirling/size term.
  In `DIV` the same size appears as m·(ph(a_t) − a_t log Ω) with ph(y) = y log y − y and
  log Ω = γ + log(π/v) − L′/L(1,χ_C). INFERRED: Ω is the **Chowla–Selberg period** of the CM curve, so m·log(a!/Ω^a) is the
  archimedean size of the m-fold a-th derivative of θ at the CM point (remark step 2: "D^a Ω·θ(0) divisible by p^{a/p+a/p²+⋯}").
* `DIV`, `Slopes` and `SlopeIneq` are **Arakelov slope / height inequalities** (μ_max, μ_min, heights h, Harder–Narasimhan
  machinery), i.e. the Bost slope inequality (remark step 13). INFERRED: the logic is "a nonzero evaluation map has total
  divisibility ≤ archimedean size". The arithmetic content is the **product formula**: for a nonzero algebraic integer N,
  Σ_p v_p(N) log p ≤ log |N|_∞.

**So the sign is: (forced divisibility) ≤ (size). It is one-sided because a nonzero integer has |N| ≥ 1.** (INFERRED, but every
piece is read.)

### 3.2 Where the "determinant" is (READ / INFERRED / OPEN)
* Remark step 4, "W(1_S) ≥ 0 via a determinant, interpret as Stepanov with multiplicity one", and step 7, "determinants for W(·)'s
  containing the desired periods, GGP periods…". I did **not** locate the specific step-4 lemma on the main chain (**OPEN**).
  The nonvanishing that the slope inequality needs is a determinant of periods (`DivisorInequality/StatementObj.lean`
  `BoundsMu`: `IsUnit (…).det`, `gram (XQ t) v0 (eb t y)).det`, HN slopes). That is READ but not traced to its proof.
* The cleanest template in the tree is **`Formalism.A.GramZ`** (`Formalism/SmithDeterminants/StatementGram.lean`, READ; it is a
  **stated Prop**, with `AssembleZ` reducing it to five sub-Props in `DefsGram.lean`). I found **no theorem proving it and no use
  on the chain**. It says: with the **Fock/Bombieri pairing** ⟨f,g⟩ = Σ_μ f_μ g_μ μ! on integer polynomials, the Gram matrix of
  linearly independent f_a = E_a·x^{α_a} is **positive definite**, so det > 0 and det is an **integer**, hence det ≥ 1. It is
  **divisible by ∏_a α_a!** (the multiplicity factorials) and ≤ ∏ diagonal (Hadamard). Hence
  v_p(det) ≤ Σ_a [v_p(α_a!) + e_a log(k²s)/log p]. This is "Stepanov with multiplicity, written as a Gram determinant": the
  Fock pairing ⟨x^α, f⟩ = α!·coeff = ∂^α f(0) **is differentiation**, so Gram positivity plus integrality turns "vanishing to
  order α" into "divisible by α!", capped by size.
* Step 6, "W(1_(0,1] + (~X^η)·1_(0,1/X]) ≥ 0 gives > 1 − η". INFERRED match: `NodeSizes` in
  `GrowthBoundGivesRiemannHypothesis.lean` uses **two scales of Riesz weights**: "thick" nodes with w ∈ [b, b+η+1], i.e. scale
  x = a/b, and "thin" nodes with w ≤ 2b′, where b′X′ ≤ xb, i.e. the long scale X′. The thin nodes' character sums are
  controlled by `ThinScales` (‖GPsi(y)‖ ≤ yθ for y ≥ X′/3, i.e. PNT-quality only). Only the thick part produces the Gz term.
  This is the "indicator at scale 1 plus X^η × indicator at scale 1/X" test function, carried out in the Lean as a family average.

---

## 4. The multiplicity amplifier and the √X (task 2)

### 4.1 Where m lives (READ)
* `ChainTied` (`Formalism/SharpRestriction/GlobalSchurBounds/TieGlobal.lean:97`) and `CoreAtH`
  (`…/PrimitiveCosets/CoreAtH.lean:29`) carry `m f : ℕ`, `Dd : Data U0 W Top m f` and weights `w : H → Fin m → E0×H → ℝ`, one
  per lift j < m. `ASK` and `Main` are summed over `j : Fin m` (`Ident`/`IdentU`).
* The theta-lift tower: `Development/AbelianVarietyTheory/BlockTheory/TwistedDiscreteSeries/IteratedLiftChains.lean`
  (`liftSpace k := span {F | IsIterLift θ ψ th ν 𝓕 σ k F}` on the unitary group `Yc D k`, with `DividedPowerChartLift`
  coefficients). INFERRED: the j-th lift is U(X_j) → U(X_{j+1}), π ↦ π ⊞ μ_j (remark step 8). Divided powers are where the
  factorials, i.e. the multiplicity, are kept integral.

### 4.2 The size constraints that force m ≈ √x (READ, PROVED [A])
* `RestSize` (`Formalism/SizeWindowBrackets/SChain.lean:53`) requires **x^{1/2−ε} ≤ m**, k ≤ n²·x^{1/2−ε}, a_E ≤ 2xb,
  n_T ≥ x^{−ε/8} n/2 and n_T2/n ≥ x^{−ε}.
* `Deduction.ch_core` (`…/LandauTheorems/ChCore.lean:343`) requires **m ≤ x^{1/2+ε/2}** and x·m ≤ k. Its conclusion is
  GmA ≤ x^{1/2+(C+3)ε} + SigA.
* `errorBracketSizeWindows` (`Formalism/SizeWindowBrackets/ErrorFormsWeakWindowsFalse.lean:133`) is
  n_T·a_E·(1 + m²/k) + m·(2n·m·b·x^ε + n_T2·k·b·x^ε/(n²x^{−ε})).
* `GrowthSizeBound`, proved as `proposition_size_bound` (`…/ChainSlopeWindows/LipschitzWeightAndMinimum.lean:507`), says:
  n_T·m·b·(δ₂/2)·(Gx − 4x^{1/2+ε}) ≤ Err + 2n_T·m·b/x with the sizes above ⇒ **Gx ≤ x^{1/2+(C+3)ε}**.

### 4.3 Why √x: the arithmetic (INFERRED from the formulas above)
The signal is Gz(x) × n_T·m·b·δ₂. Divide the error bracket by it:
* **trivial term:** n_T·a_E / (n_T m b δ₂) ≈ 2x/(m δ₂). It is the size-x "no cancellation" loss, and **it does not grow with m**.
* **per-lift cost:** 2n·m²·b·x^ε / (n_T m b δ₂) ≈ m·x^{O(ε)}. **Each lift costs about one unit.**
* **degree term:** ≈ k·x^{O(ε)}/(n·n_T), which is ≤ x^{1/2+O(ε)} by k ≤ n²x^{1/2−ε}.

So Gz(x) ≲ x^{O(ε)}·(x/m + m + √x), minimised at **m = √x** (AM–GM: x/m + m ≥ 2√x). That is exactly Alpöge's line "length of
chain of theta lifts … was √X (and yet, what else could it be?)". Fewer lifts leave the trivial loss above √x. More lifts make
the lift costs exceed √x.

### 4.4 "qN = Mk + R" in the Lean (READ / INFERRED)
* The identity half: `AskMainBalance`, proved as `colour_balance` (`…/ChainSlopeWindows/StronglyIrreducible.lean:695`), says
  **ASK − MAIN = De**. De is a finite-Fourier expansion over the characters ψ of E0 × H of the twisted Riesz sums GPsi. It comes
  from `characterSum_expansion` (`AskExpansion`), `orthogonality_unramified` (`OrthoUnram`), `colour_sum` and
  `partner_bijection`, all PROVED [A]. This is the "iterate Parseval" step: Gauss-sum cancellation is replaced by exact
  orthogonality on a finite abelian group.
* `ThickPartTranslateMean` (READ; a stated Prop with no proof or use found anywhere in the tree, so the on-chain counterpart is inside `ch_core`'s proof, not traced, OPEN) says the thick part of De, averaged over u₀, equals
  (δ₂ b/2)·Gsig(A/b), where **Gsig = G_0 + G_ψC − G_χC − G_χψC**. INFERRED: this is the four-character combination with
  coefficients (1+ψ_C)(1−χ_C) ∈ {0, 4}, i.e. **non-negative**, and the trivial character's term is Gz.
* So "qN = Mk + R" becomes: [definite-sign arithmetic count (ASK ≥ 0, ≤ size by §3.1)] = [main term] + [Gz × mass] + [other
  characters]. The other characters are R, bounded by `ThinScales`/`ValuesAtX` at ≤ x^{1/2+ε}. "Opposite signs mean no
  cancellation" holds because ASK can only be **too small**, never too large, relative to size.

---

## 5. Translation into Daniel's language — **HEURISTIC throughout**

1. **Is the determinant a Gram matrix (‖Af‖² shape)?** *Partly.* The template (`GramZ`) is a Gram determinant of a
   positive-definite pairing, so yes, "x² is positive whichever sign x has" is literally there, as `PosDef ⇒ det > 0`. But the sign
   that does the work is **not** "a square is ≥ 0". It is "**a positive integer is ≥ 1**". Positivity plus **integrality**
   (quantization) is what gives a one-sided bound. A real Gram determinant could be 10⁻¹⁰⁰⁰. An integer one cannot.
   In Daniel's words: **the positivity comes from "counted once" (whole numbers, no fractional counts), not from squaring.**
2. **Is the amplifier "k copies of the boundary pairing"?** *Yes in structure, with one sharp correction.* The m lifts are m copies
   of the same divisibility pairing (ASK sums over j < m). The amplification works only because the **trivial error is shared and
   does not copy**, while each copy adds the same signal. The relevant ratio is x/m + m, not "error shrinks like √ each time". The
   ++/−− picture fits each copy (one copy = one sign-definite count). The win comes from stacking copies against a fixed noise
   floor, like averaging m coherent detectors against one shared background.
3. **Conservation law (t497u) = product formula.** Σ_{all places v} log|N|_v = 0 for a nonzero rational N. Each prime place can
   only have |N|_p ≤ 1 ("no negatives" at the finite places), so the archimedean place must carry the balance. This is the cleanest
   Daniel-language reading of Alpöge's sign: **a conservation law across places, with every finite place one-signed.**
4. **Can it be aimed at our D(L) to give F13's hypothesis?** The **exit matches**: F13's D(L) ≥ −C e^{εL} is the same direction
   as his Gz(x) ≤ x^{1/2+ε} (an upper bound on a smoothed prime sum), and his `FormB` says windows suffice. **The engine does not
   transfer.** D(L) is a real-analytic function of L built from Λ. It has **no integer attached whose p-adic divisibility is forced
   and whose size is controlled**. Every piece of Alpöge's sign needs such an integer: CM theta derivatives, inert-prime
   D^p ≡ 0, slope inequality. Nothing in our R13–R28 benches (Weil forms, screw function, zitter/spinor pictures, four-sector
   holography) supplies one. Our own Round 16 log already named this gap: "GAP = integer Stepanov lemma (square-root-sharp
   auxiliary integers)". Alpöge's tree is, as far as I can read it, exactly a (claimed) construction of such auxiliary integers,
   done automorphically.

**Blunt verdict.** What we can steal cheaply: (a) the **windowed Landau exit**; (b) the **shape** of the argument: a
one-sided bound comes from *integrality × multiplicity against a shared trivial bound*, with chain length √x by AM–GM; (c) the
**finite-Fourier (Parseval) explicit formula** in place of analytic cancellation. What we cannot steal without his whole
arithmetic engine: the integers themselves. A physics proof in Daniel's language would need a *physical reason the relevant
quantity is quantized* (an integer, an index, a count of states) and a *physical reason its local (prime) obstructions are
one-signed*. That is a precise target, and it is not a W ≥ 0 target.

---

## 6. OPEN and proposed next steps (none aimed at W ≥ 0)

* **OPEN-29a:** locate the step-4 determinant lemma on the main chain: which period determinant is shown ≠ 0, and how.
  Lead: `TopG`/`shapeWin` → `link_irr hGreen` → the `DivisorInequality`/`HarderNarasimhan` statements.
* **OPEN-29b:** confirm that `ChainTied`'s m is the theta-lift depth k in `IteratedLiftChains` (currently INFERRED).
* **Proposal R29-L (Lean, small, ours):** *F13-window*. Prove that D(L) ≥ −C e^{εL} on infinitely many windows
  [L₀, (3/2)L₀] suffices for ξ ≠ 0 on Re s > ½ + ε. Port the idea of `zeroFree_of_rpow_window` to PR #649's setting
  (Claude's lane; to be offered via the conversation file, not filed here).
* **Proposal R29-T (toy, numerical):** reproduce "signal ∝ m, trivial loss shared, cost ∝ m" with Stepanov over F_q
  (R16 bench), and check that the best multiplicity scales like √q.
* **Proposal R29-Q (physics, Daniel's):** "Where is our integer?" Find a quantity in the s = Δ / boundary picture that is
  **quantized** (an index, a degree, a state count) and whose prime-local contributions are one-signed. The product formula is the
  template conservation law.

---

*Sources: GPP-bridge `research/external/alpoge-rh/` (remark.pdf, key-files, lean-parts; sha256 above), and Claude Code (D1)'s notes
in GPP-bridge `CLAUDE_CODE_RESEARCH_NOTES.md`. Nothing in this file was kernel-checked by this seat. All Lean citations are to
Alpöge's source as released.*

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
