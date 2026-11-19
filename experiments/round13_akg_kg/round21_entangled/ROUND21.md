# Round 21: "Two prime bulks entangled through x = 1": the thermofield double

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). 2026-11-19 (Toronto).

**Frame.** d = 1, s = Δ, states in H = L²(ℝ>0, dx), shadow (Jf)(x) = x⁻¹f(1/x) (Δ ↔ 1 − Δ). In log coordinates x = eᵘ, (Vf)(u) = e^{u/2}f(eᵘ) is unitary onto L²(ℝ, du) and J becomes the reflection u ↦ −u. The two halves are H_R = L²(x > 1) (u > 0) and H_L = L²(x < 1) (u < 0). The boundary is x = 1 (u = 0), J pairs n with 1/n. Integer state F_ζ(x) = Σ_{n≥1} e^{−πn²x²} − 1/(2x), φ_ζ(u) = e^{u/2}F_ζ(eᵘ) is even, and ∫φ_ζ(u)e^{(s−½)u}du = M_ζ(s) := ½Γ_ℝ(s)ζ(s) on 0 < Re s < 1.
**RH: OPEN. Nothing here claims it.** Zeros are never proof inputs; the two zeros located in §3a (one ζ, one Epstein) are illustrations only, and every count in §3b is computed from M alone. Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN. "Cited" = published theorem stated, not re-derived.

**Daniel's idea.** "Two connected prime bulks, each 1/x with its partner x from the other half, by entanglement. The zero is the connections passing through the boundary."

## Short answer

1. **The state exists and is natural (PROVED).** There are two exact versions:
   - the **Bost–Connes TFD** Ψ_β = Z^{−1/2} Σ_n n^{−β/2} |1/n⟩_L |n⟩_R, with Schmidt weights n^{−β}/ζ(β), which is the primon gas at temperature 1/β;
   - the **theta TFD** Ψ(x, y) = F_ζ(x/y)/y for x < 1 < y, the cross-boundary block of the theta transfer operator. That operator is T x^{−s} = M_ζ(s) x^{−s}. In log coordinates the block is the Hankel operator (Γg)(u) = ∫₀^∞ φ_ζ(u+v) g(v) dv.
   
   The second is the first with the Hagedorn divergence at β = 1 (the pole) Poisson-subtracted, so it exists down to the critical line.
2. **Positivity (Task 2).** Each TFD fails differently:
   - BC TFD: positivity is automatic (a sum of squares), and the zeros never enter. Trivially not RH.
   - Theta TFD: its reflected form is **indefinite for every model**, ζ included (PROVED via Widom's theorem, plus a 25-digit CHECKED witness). Wrong level.
   - The reflected form that is RH-equivalent is the **Weil form**: the same TFD construction, but with the theta kernel replaced by the kernel of log ξ. That is 14B/Connes, not new.
3. **"Zero = connection through the boundary" is a theorem (PROVED, Riemann 1859 in new words).**
   - M(s) = A(s) + A(1−s), where A(s) = ∫₁^∞ x^{s−1}F(x)dx is the flux of the x > 1 bulk and A(1−s) is the flux of its J-partner bulk x < 1.
   - A zero is exact destructive interference at the boundary: A(ρ) = −A(1−ρ).
   - On the line the partner amplitude is automatically the complex conjugate, which is the TFD pairing. So an on-line zero is a sign flip of the real transmitted amplitude T(t) = 2 Re A(½+it).
   - An off-line zero is a cancellation against a **non-conjugate** partner.
   
   RH ⇔ every boundary cancellation is TFD-paired ⇔ #sign flips of T on (0, T] = N(T) for all T. That last form is Turing's method. NUMERICAL:
   - ζ: 29 = 29 up to T = 100;
   - EPS: 32 vs 42 up to T = 50 (5 ghost pairs);
   - DH: 48 vs 50 up to T = 95 (1 ghost pair).
4. **New, small, ζ-relevant diagnostic.** ζ's TFD is a tensor product over primes: the prime sectors are mutually unentangled (total correlation 0). That is the Euler product written in entanglement language. EPS's prime sectors are entangled (I(v₂; v₃) = 0.145 at β = 2), and DH is not even a state (10.5 % negative weight). **But F_θ passes both (positive and prime-unentangled) and has zeros off the line.** So "positive Schmidt weights + Euler product" ⇏ RH (PROVED by counterexample).
5. **GAP-21.** Prove that the integer flux A(s) can cancel its partner only when the partner is its conjugate. That is RH verbatim; the TFD gives it a physical name but no new inequality. **Link to Round 20:** a positive ladder is exactly a TFD purification with positive Schmidt weights at *real* energies (γ²). That positivity lives at the log-ξ (Weil) level, not at the ξ level where the theta TFD is indefinite. So the circuit should be built from ξ′/ξ, which is what Round 20's plan does.

Files: `r21_models.py` (states, kernels); `r21_hankel.py/.json/.log` (theta TFD spectra and witness); `r21_tfd.py/.json/.log` (Bost–Connes TFD); `r21_transmission.py/.json/.log` (boundary flux and Turing count).

---

## 1. The precise states (Task 1)

### 1.1 Theta transfer operator and its cross-boundary block (PROVED)

Define (Tf)(x) = ∫₀^∞ F_ζ(x/y) f(y) dy/y, with kernel K(x, y) = F_ζ(x/y)/y. Then T is multiplicative convolution, and for 0 < Re s < 1
  T x^{−s} = M_ζ(s) x^{−s}  (substitute x/y = w).
So every primary is a transfer eigenmode, and its eigenvalue is ½Γ_ℝζ. In log coordinates T is convolution by φ_ζ.

Split T into blocks over H_L ⊕ H_R. The cross block T_LR : H_R → H_L, composed with J : H_L → H_R, has log-coordinate kernel φ_ζ(u + v) on u, v > 0. That is the **Hankel operator Γ_φ** on L²(0, ∞). As a two-particle wavefunction it is
  **Ψ(x, y) = K(x, y)·1{x < 1 < y} = F_ζ(x/y)/y**,
which is Daniel's "bulk x > 1 connected to bulk x < 1 through x = 1". Its quadratic form ⟨f, Γ_φ f⟩ = ⟨Θf, T f⟩ with Θ = J∘(complex conjugation) is exactly the Osterwalder–Schrader reflected pairing of T across x = 1.

### 1.2 Bost–Connes TFD (PROVED, classical: Julia 1990, Bost–Connes 1995; cited)

  Ψ_β = ζ(β)^{−1/2} Σ_{n≥1} n^{−β/2} |1/n⟩_L ⊗ |n⟩_R,  β > 1.
The reduced state on R is the primon-gas Gibbs state ρ_R = Σ n^{−β}/ζ(β) |n⟩⟨n| with H = log n. The Euler product makes Ψ_β = ⊗_p (Σ_k p^{−kβ/2}|p^{−k}⟩|p^k⟩)·norm, a product of independent TFDs, one per prime (free bosonic primons). The return amplitude ⟨Ψ_β| e^{−itH}⊗1 |Ψ_β⟩ = ζ(β+it)/ζ(β) is the spectral form factor. The entanglement entropy S(β) = log ζ(β) − βζ′(β)/ζ(β) diverges at β = 1: **the pole is a Hagedorn point, where the boundary carries infinite entanglement.** That is Daniel's "the pole is the infinity" (Round 18/19) in TFD form. NUMERICAL, `r21_tfd.json`: S = 1.638, 3.218, 12.74, 105.0, 1007 at β = 2, 1.5, 1.1, 1.01, 1.001.

### 1.3 How the two are related (PROVED)

For Re s > 1, ∫x^{s−1}Σ_{n≥1}e^{−πn²x²}dx = ½Γ_ℝ(s)·Σn^{−s} = ½Γ_ℝ(s)·Z_BC(s). So the theta series is the BC partition function in position space, dressed by the Tate local factor. Subtracting the Poisson term −1/(2x) (the Hagedorn volume) continues it into 0 < Re s < 1. **The theta TFD is the renormalized BC TFD.** The BC TFD does not exist at β = ½, but the theta TFD does, which is why Daniel's picture has to be run with F_ζ.

## 2. Entanglement spectrum and signs (Task 2)

### 2.1 Theta TFD: Schmidt values and signs (NUMERICAL; `r21_hankel.json`)

Nyström, Gauss–Legendre on [0, L], N = 200 and N = 400, agreeing to all printed digits. The Schmidt coefficients are the |eigenvalues| (≥ 0 automatically, so they carry no information). The **signs** are the eigenvalues of the self-adjoint Γ_φ, i.e. of the reflected form.

| model | φ on (0, ∞) | top eigenvalues of Γ_φ (signed) | # pos / # neg (above 10⁻¹³·max) | trace (= ∫₀^∞φ(2u)du) |
|---|---|---|---|---|
| ZETA | < 0 everywhere | −0.4994, +2.43e−3, −1.37e−4, +1.20e−5, −1.32e−6, +1.68e−7, … | 7 / 7 | −0.49712 |
| EPS x²+5y² | < 0 everywhere | −0.9218, +7.19e−2, −1.07e−2, +1.67e−3, −2.83e−4, … | 10 / 10 | −0.85914 |
| DH | > 0 | +0.2082, −3.52e−2, +7.66e−3, −1.86e−3, +4.62e−4, … | 13 / 13 | +0.17915 |
| F_θ (θ = 0.2) | > 0 | +1.464, −6.09e−2, +3.97e−5, −1.34e−6, … | ≥ 2 / ≥ 2 (rest is grid noise ≲ 10⁻⁷) | +1.40355 |

**Every model, ζ included, has an indefinite reflected form with alternating, geometrically decaying signs.** The dominant eigenvalue in each case is the Poisson (Hagedorn) volume mode: −½e^{−u/2}e^{−v/2} for ζ, giving −0.4994 ≈ −½. ζ is not distinguished.

**PROVED (ζ, EPS; DH below).** Γ_φ is not ≥ 0, because φ < 0 on (0, ∞): any f > 0 gives a negative value. For ζ, φ(w) = e^{w/2}θ₊(e^{2w}) − ½e^{−w/2}. The first term is decreasing on w ≥ 0 from θ₊(1) = 0.0432 < ½. And −Γ_φ = Γ_{−φ} is not ≥ 0 either, by **Widom's theorem** (cited: a Hankel integral operator on L²(0, ∞) is ≥ 0 iff its kernel is the Laplace transform of a positive measure, hence completely monotone, hence analytic in Re w > 0). −φ is not analytic in Re w > 0: θ₊(e^{2w}) has the natural boundary Im w = ±π/4 (the real line of the modular variable). The same argument applies to EPS. For DH (kernel > 0) the kernel decays like e^{−c·e^{2w}}, faster than any nonzero completely monotone function, so Γ is not ≥ 0 (PROVED), and Γ is not ≤ 0 since φ > 0. For F_θ, indefiniteness is NUMERICAL (the −6.09e−2 eigenvalue is stable in N).
**CHECKED (25 digits), ζ witness.**
- ⟨e^{−u/2}, Γ_φ e^{−u/2}⟩ = −0.499418072132222.
- With f = e^{−3u} − (2/7)e^{−u/2} (f ⊥ e^{−u/2}), ⟨f, Γ_φ f⟩ = +0.000156648261960352.
- Opposite signs on the same operator.

**Verdict on the theta TFD.** Positivity is not automatic, and it is *false* for everyone. The theta TFD is at the **ξ level**: its symbol on the line is the transmission amplitude T(t) = M(½+it), which is real and **changes sign at every on-line zero**. A sign-changing symbol cannot give a positive form, so in this picture the *good* zeros are what destroy positivity. That is why it cannot be RH-equivalent.

### 2.2 The RH-equivalent reflected form is the Weil form (PROVED as an equivalence; Weil 1952, Round 14B; cited)

Replace φ by the Weil distribution W, the kernel of the explicit formula: W(g) = Σ_ρ M[g](ρ), built from Λ(n), the Γ-factor and the poles only. Then ⟨Θf, W ∗ f⟩ = Σ_ρ M[f](ρ)·conj M[f](1−ρ̄).
- On-line zeros give |M[f](ρ)|² ≥ 0.
- An off-line pair gives 2 Re(M[f](ρ) conj M[f](1−ρ̄)), which is indefinite.

So reflection positivity of the *Weil* TFD ⇔ RH (Weil; Round 14B "No-Ghosts ⇔ RH"). Relation to the theta TFD: W is the log level. Under RH its symbol on the line is Σ_γ δ(t − γ), the boundary value of (1/π)Re ξ′/ξ(½+0+it), i.e. the derivative of the phase of the transmission amplitude (a Wigner-type time delay). The theta TFD's symbol is the amplitude itself. Off-line zeros have no line symbol at all; they enter W only through analytic continuation of M[f], which is exactly where indefiniteness comes in. **Daniel's entangled state is reflection-positive exactly when it is built from log ξ. That is the Weil form, i.e. 14B/Connes. It is not a new positivity.**

### 2.3 BC TFD: automatic positivity, prime-sector entanglement (CHECKED/NUMERICAL; `r21_tfd.json`, n ≤ 10⁶)

TC = Σ_{p ≤ 13} H(v_p) + H(rest) − H(n) is the total correlation between the prime sectors. It is zero iff the TFD factorizes over primes.

| model | coefficients | Schmidt weights ≥ 0? | S(β=2) | I(v₂; v₃), β=2 | TC, β=2 | TC, β=3 | zeros on line? |
|---|---|---|---|---|---|---|---|
| ZETA | all ones | yes (automatic) | 1.6376 | 7e−12 | 3.9e−7 (truncation) | 2.3e−11 | OPEN |
| EPS | r_Q(n) ≥ 0 | yes | 0.9710 | **0.1451** | **0.2474** | **0.0568** | no (5 ghost pairs ≤ 50) |
| F_θ (θ=0.2) | Σ_{ab=n}(a/b)^θ > 0, multiplicative | yes | 2.7349 | 1e−9 | 1.6e−5 (truncation) | 3e−11 | **no** (zeros at ½ ± θ + iγ, PROVED) |
| DH | 1, κ, −κ, −1, 0 (mod 5) | **no**: 10.5 % (β=2), 2.7 % (β=3) of the weight is negative | — | — | — | — | no (1 ghost pair ≤ 95) |

Additivity check, ζ: S(3) = 0.678502221866 (exact formula) vs Σ_{p ≤ 2·10⁵} S_p(3) = 0.678502221828 (CHECKED).

- **PROVED (logic).** For real β > 1 the BC TFD is a positive state whenever c(n) ≥ 0, and ζ(σ+it) ≠ 0 for σ > 1 by the Euler product. The zeros enter only at complex β, where the "overlap" is a return amplitude and not a norm, and only for β < 1, below the Hagedorn point, where the state does not exist. The Rényi entropy S_α = log(ζ(αβ)/ζ(β)^α)/(1−α) is log-singular exactly at the complex Rényi indices α = ρ/β. So the positivity is automatic and trivially not RH.
- **PROVED (counterexample).** F_θ = ζ(s−θ)ζ(s+θ) has a positive, prime-factorized TFD and the shadow symmetry (Round 13/19). It vanishes at ρ ± θ for every ζ-zero ρ, and since those two real parts differ by 2θ ≠ 0, at least one is off Re s = ½. This is unconditional. So **positive Schmidt weights + Euler product + shadow symmetry ⇏ zeros on the line.** What F_θ breaks is "all ones" (degree 1), exactly as in Round 19.
- **Modestly new.** "Prime sectors unentangled (TC = 0)" is an entanglement reading of the Euler product. It cleanly flags EPS (TC = 0.247) and, with positivity, DH. It is a reformulation of multiplicativity, not a new constraint on zeros.

## 3. "The zero is the connection through the boundary" (Task 3)

### 3a. Boundary-flux decomposition (PROVED: Riemann 1859; CHECKED)

Split the Mellin integral at the boundary x = 1 and use J F = F:
  M(s) = A(s) + A(1−s), A(s) := ∫₁^∞ x^{s−1}F(x)dx (entire in s).
- ζ: A(s) = ∫₁^∞ x^{s−1}θ₊(x²)dx − 1/(2(1−s)).
- EPS: A(s) = ∫₁^∞ t^{s−1}(Θ_Q(t) − 1)dt − 1/(1−s).

A(s) is the right bulk's flux into the boundary, and A(1−s) = ∫₀¹x^{s−1}F dx is the partner bulk's flux. CHECKED: |A(s) + A(1−s) − M(s)| ≤ 1.1·10⁻²² at three test points for ζ and EPS.

**Zero ⇔ A(ρ) = −A(1−ρ)**: the two bulks' fluxes cancel exactly at the boundary. On the line 1 − s = s̄, so the partner flux is *automatically* conj A(s) (the thermofield pairing), and T(t) = M(½+it) = 2 Re A(½+it) is real. Hence (PROVED):
- an **on-line zero** ⇔ Re A = 0 ⇔ arg A(½+iγ) = ±π/2: the right flux is purely "imaginary" and its conjugate partner cancels it. The transmitted amplitude changes sign (for odd multiplicity).
- an **off-line zero** is a cancellation A(ρ) = −A(1−ρ) with A(1−ρ) ≠ conj A(ρ): the bulk cancels against a partner that is *not* its TFD mirror.

Illustration (zeros located only to illustrate; `r21_transmission.json`):

| | ρ | A_right = A(ρ) | A_left = A(1−ρ) | \|A_R + A_L\| | \|A_L − conj A_R\| | arg A_R / π |
|---|---|---|---|---|---|---|
| ZETA | ½ + 14.1347251417i | −0.0326391042 i | +0.0326391042 i | 0 | 0 | −0.5000000000 |
| EPS | 0.932969697485 + 15.6682495313i | 0.00087351 − 0.03176927 i | −0.00087351 + 0.03176927 i | 1e−22 | **1.75e−3** | −0.4912501734 |

So Daniel's sentence becomes the exact statement: **RH ⇔ every connection through the boundary is TFD-paired (partner = conjugate)**. This holds for EPS's on-line zeros too; the EPS ghosts are the unpaired connections.

### 3b. Counting connections: sign flips vs the argument principle (NUMERICAL; Turing 1953 / Backlund method)

Here "flips" means the sign changes of T(t) on a 0.01 grid, and N(T) is the argument-principle count of zeros of ξ_M in [−1, 2] × (0, T], from M alone. No zeros were used.

| model | T | flips (TFD-paired connections) | N(T) | deficit | zeros right of Re s = 0.52 |
|---|---|---|---|---|---|
| ZETA | 100 | 29 | 29.0 | 0 | 0.0 |
| EPS | 50 | 32 | 42.0 | 10 | 5.0 |
| DH | 95 | 48 | 50.0 | 2 | 1.0 |

The deficit is always 2 × (ghosts right of the line), as the ρ ↔ 1 − ρ̄ pairing requires. RH up to height T ⇔ deficit = 0. This is the standard rigorous verification method, here read as "every zero is a boundary connection that flips the transmitted amplitude".

## 4. Verdict, what is new, GAP

**Plain words.** Daniel's picture is exactly right as a *description*.
- ξ really is the overlap of two mirror-image bulks glued at x = 1.
- A zero really is their fluxes cancelling at the boundary.
- On the critical line the partner is automatically the mirror image (conjugate), so a zero there is a clean "connection" that flips the sign of what gets through.

What the picture does not do is *force* that pairing:
- The automatically positive version (Bost–Connes) never sees the zeros.
- The version that sees the zeros at the amplitude level (theta) is indefinite for everybody.
- The version whose positivity is equivalent to RH is the Weil form we already had.

The fakes show the boundary-cancellation statement is shared (EPS satisfies A(ρ) = −A(1−ρ) at its ghosts) and that positive, prime-factorized entanglement is shared (F_θ).

**New here (small, honest).**
- (i) The theta transfer operator's cross-boundary block is a Hankel operator, and Widom plus the theta natural boundary makes ξ-level reflection positivity impossible for every model. This is the clean reason the ξ-level TFD is the wrong level.
- (ii) The TFD reading of Riemann's split and of Turing's count: "zero = TFD-paired boundary cancellation, ghost = unpaired".
- (iii) Prime-sector total correlation as an entanglement form of the Euler product, plus the F_θ counterexample showing it cannot suffice.

**Repackaging.** Bost–Connes / primon gas, Riemann's 1859 split, Weil positivity (14B/Connes), and the Turing method.

**GAP-21.** Show that the integer flux A_ζ(s) = ∫₁^∞ x^{s−1}F_ζ(x)dx satisfies A_ζ(ρ) = −A_ζ(1−ρ) ⇒ 1 − ρ = ρ̄. This is RH verbatim. Any proof must use something that F_θ (all-ones fails), EPS (Euler fails) and DH (both fail) lack, and it must produce an inequality, the same positivity-transfer gap as GAP-19.

**Link to Round 20 (HEURISTIC/PROVED-classical).** A Stieltjes function ∫dμ(λ)/(λ+z) with μ ≥ 0 is ⟨Ψ|(H+z)⁻¹|Ψ⟩ for a TFD purification with Schmidt weights μ_k ≥ 0 at energies λ_k. In z = (s−½)², under RH ξ′/ξ = Σ 1/(z + γ²): unit Schmidt weights at real energies γ², and an off-line zero makes an energy complex. So "positive ladder" = "positive Schmidt weights at real energies", and it lives at the **log-ξ level**. Round 21 shows the ξ-level (theta) TFD is indefinite for every model. The circuit therefore has to be built from ξ′/ξ, not ξ, and Round 20's continued fraction of the log-derivative is the right object. Round 21 informs 20 more than the reverse.

## 5. Status ledger
- PROVED (classical/cited): T x^{−s} = M(s)x^{−s}; cross block = Hankel Γ_φ = OS form; BC TFD and Euler factorization; Hagedorn at β = 1; renormalized-BC = theta; Riemann's split M = A(s) + A(1−s); on-line zero ⇔ arg A = ±π/2; Weil RP ⇔ RH; Widom's criterion.
- PROVED (here): ξ-level (theta) reflection positivity fails for ζ, EPS (φ < 0, plus Widom with the natural boundary) and DH (super-exponential decay); NUMERICAL for F_θ; BC positivity + Euler product + shadow symmetry ⇏ RH (F_θ).
- CHECKED: flux split (10⁻²²); ζ Hankel witness of both signs (25 digits); S additivity.
- NUMERICAL: Hankel spectra (table 2.1); TC/I table (2.3); Turing counts (3b); EPS unpaired connection (3a).
- HEURISTIC: Round 20 link; "Hagedorn = pole = Daniel's infinity".
- GAP-21 as above. OPEN: RH.

## References
Riemann (1859) · Weil (1952) · Turing, Proc. LMS (1953) · Widom, "Hankel matrices", Trans. AMS 121 (1966) · Julia, "Statistical theory of numbers" (1990) · Bost–Connes, Selecta Math. 1 (1995) · Connes, Selecta Math. 5 (1999) · Davenport–Heilbronn (1936) · Rounds 14B, 17, 18, 19 (`ed0a255`).

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
