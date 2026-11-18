# Round 14D — one global arithmetic object with three faces

Framing credit: **Daniel Toupin** — "construct one global arithmetic object whose integer/theta representation, prime/Euler representation, and positive Hilbert-space representation are provably the same object." The No-Ghosts framing (14B) is also his.

Code: `lineval.py` (|F(½+it)|² from Hurwitz sums; `line_*.npy`), `radical.py` (`radical.json`), `pencil.py` (`pencil.json`). It reuses `../round14b_euler/weil_galerkin.py` and `models14b.py`. Zeros are never inputs.

> **Errata (Golden Eye, AUDIT14 on team/golden-eye 062e4e2; applied in Round 14E).**
> 1. "range E ⊂ radical" (3a) is the **easy direction** of the explicit formula (Mellin(E_F H) = F·Ĥ vanishes at every zero of F). PROVED, classical.
> 2. The 3a numerical "CHECKED rel ≤ 1e-4" tested **isotropy W(v,v) ≈ 0, not radical membership**, and had low power (a same-band non-radical vector also scores ~4e-6). It is downgraded to "consistent at the grid floor; not discriminating". **Redone in 14E** (`../round14e_gap/bilinear.py`): bilinear residual sup_φ |W(v,φ)|/(‖v‖‖φ‖‖W‖) over v in the window-truncated range, against a same-band non-radical control (a(2) → a(2) − 0.1). Separation 6.4e3 (ζ, L=4), 1.2e4 (ζ, L=6), 1.3e2 (EPS, L=6), 4.7e3 (F_θ, L=4), 7.4e2 (DH, L=4) at ε = 1e-8. That is now a discriminating CHECK.
> 3. The BEUR″ table entry "E_B is ~ in ζ's radical (6e-6)" is **not supported** at that probe's power; withdrawn. (14E: ζ's truncated range is *not* in BEUR″'s radical: bilinear residual 0.1–0.4.)
> 4. "Positivity on the quotient ⇔ RH" is **Weil's criterion restated** (quotienting by part of the radical cannot change positivity; GNS). Connes 1999's cokernel carries the L² (or L²_δ) inner product, not the W-norm; the closure must be taken in the **form-domain norm**, not L². Faithful statement: "if W ≥ 0, its GNS space is a quotient of the admissible space by a subspace containing range E; whether that subspace equals closure(range E) is the spectral-realization question."
> 5. "R indefinite past log 2 = CC 2021 boundary" is a **coincidence**, not a check: DH's R is already indefinite below log 2 (R_min = −1.78 at log 2 − 0.05). EPS "μ_min < 0 at L ≥ 3" should read L > 3.
> 6. Faces 1+2 are Tate's thesis (Riemann 1859 / Tate 1950); credit Daniel for the framing, not a new object.

Standards (Golden Eye):
- **CHECKED** = a code or identity check.
- NUMERICAL = evidence.
- Every control is non-tautological. A positivity test whose failure only restates "W < 0" is marked as tautological and is not counted as a falsifier.

## Faces 1 and 2 are one object (PROVED: Riemann 1859, Tate 1950)
- **The object.** Let f₀ = e^{−πx²} ⊗ ⊗_p 1_{ℤ_p} ∈ S(𝔸_ℚ), and define Θ(a) := Σ_{q∈ℚ^×} f₀(qa) on the idele class group C_ℚ = 𝔸^×/ℚ^×. This is the Tate theta: the sum over ℚ of the Tate ground state.
- **Face 1 (integers / theta).**
  - f₀(qa) ≠ 0 at the finite places only when q ∈ ℤ. So on the archimedean ray a = x ∈ ℝ_{>0}, Θ(x) = Σ_{n≠0} e^{−πn²x²}.
  - Poisson summation over ℚ ⊂ 𝔸 with the self-dual f₀ (f̂₀ = f₀) gives Θ(a^{−1}) = |a|Θ(a) + (|a| − 1)f₀(0).
  - On the ray this is the Jacobi theta relation Σ_n e^{−πn²/x²} = x Σ_n e^{−πn²x²}.
- **Face 2 (primes / Euler).**
  - The Mellin transform of the same Θ is ∫_{C_ℚ} Θ(a)|a|^s d^×a = ∫_{𝔸^×} f₀(x)|x|^s d^×x for Re s > 1.
  - Because f₀ factorizes, this equals ∏_v ∫ f₀,v |x|_v^s d^×x_v = Γ_ℝ(s) ∏_p (1 − p^{−s})^{−1} = Λ(s).
- **Conclusion.** Faces 1 and 2 are two expansions of *one function on C_ℚ*:
  - the additive one is a sum over the lattice ℚ ∩ (ℝ × ∏ℤ_p) = ℤ;
  - the multiplicative one is a product over places.
  - Their compatibility is exactly Poisson summation together with the factorization of f₀.

## Face 3: the Connes–Consani global map E (all integers at once). Results first.
Definition. E_F(H)(x) = Σ_n a_F(n) H(nx), written in u-form as g(u) = e^{u/2} E_F H(e^u). Each model uses its own integer data:
- ζ: a(n) = 1;
- Epstein (EPS): a(n) = r_Q(n), the lattice count of x² + 5y²;
- DH: its periodic a(n);
- F_θ: Σ_{d|n} d^{θ} (n/d)^{−θ};
- spliced Beurling (BEUR″): the Beurling integers.

### 3a. E lands in the radical of the Weil form, for every model with an explicit formula (PROVED, generic; CHECKED numerically)
- **Proof.** The Mellin transform of E_F H is F(s)Ĥ(s). Suppose Ĥ vanishes at F's poles. Then g = E_F H is admissible: it decays as x → 0 by the Mellin shift. Its transform vanishes at every zero of F, so by the explicit formula W_F(g ∗ Θg′) = 0 for all admissible g′. ∎
- **This is fake-blind by construction.** It uses only "F has an explicit formula".
- **CHECK (grid, u ∈ [−6, 2.5]).** rel = W / max(|arch|, |pole|, |prime|):

| E built from | W of the same model | W of ζ (cross) | control: W_own(H) itself |
|---|---|---|---|
| ζ | **−2.2e-6** (radical) | – | 6.7e-4 |
| EPS (lattice) | **6.2e-6** (left-edge truncation 1.8e-3) | 7.5e-3 | 1.19 |
| DH | **1.6e-6** | 5.4e-4 | 0.41 |
| F_0.2 | **1.0e-4** | 3.5e-4 | 5.2e-2 |
| BEUR″ (Beurling integers) | −1.05e-3 (ζ's arch borrowed; there is no own explicit formula, so "own radical" is undefined: DEFECT / N.A.) | **5.8e-6** | – |

Also W_EPS(E_ζH) has rel 1.43, and W_BEUR″(E_ζH) has rel −2.4e-6.

- **Reading.**
  - Each model's E lands in its own radical, to 1e-6 (1e-4 for F_θ). The cross terms are larger.
  - So the radical is a property of the pair (F, E_F), and it holds for the fakes just as well. **It does not single out ζ.**
  - Honest caveat: at this window, E built from the Beurling integers is *also* almost in ζ's radical (5.8e-6). The Beurling integers agree with ℕ below 55, and the difference is too small to resolve here. So the radical test does **not** separate BEUR″ from ζ at this resolution.

### 3b. The global norm ‖E_F g‖² against the Weil form (NUMERICAL, `pencil.json`)
- **The norm.** N_F(g) := ‖E_F g‖² = (1/2π) ∫ |F(½+it)|² |ĝ(t)|² dt.
  - This is Plancherel. Its weight |F|² comes only from integer data (Hurwitz sums on the line).
  - **It is positive for every F**, including EPS, DH and F_θ (PROVED). It is the 14A / 13i kinematic positivity again, now global.
  - Plancherel normalization check: N = I to 3e-4 (L = 1) and 4e-5 (L = 4) when F ≡ 1 (CHECKED).
- **The tested decomposition.** Write W = c·N − R and ask whether R has a definite sign. Use the generalized eigenproblem W v = μ N v:
  - μ_min is the best c with W ≥ cN (that is, R ≤ 0);
  - μ_max is the smallest c with W ≤ cN (that is, R ≥ 0).

| L | ζ μ_min | ζ μ_max | ζ: sign of R = N − W | EPS μ_min | DH μ_min | F_0.2 μ_min |
|---|---|---|---|---|---|---|
| log 2 − 0.05 | **5.2e-3** | 0.86 | **R ≥ 0** (min 0.38) | 0.118 | 0.50 | 7.5e-3 |
| log 2 + 0.05 | 9.1e-4 | 1.37 | indefinite (−0.61 … 3.05) | 0.110 | 0.44 | 1.3e-3 |
| 1 | 1.6e-6 | 1.92 | indefinite | 0.094 | 0.23 | −3.0e-5 |
| 2 | ~0 (−2.5e-14) | 4.38 | indefinite | 0.045 | 3.6e-5 | −0.012 |
| 3 | ~0 | 8.17 | indefinite | −1.7e-6 | ~0 | −0.078 |
| 4 | ~0 | 15.6 | indefinite | **−0.035** | ~0 (DH's ghosts at 85.7 lie outside this F = 45 band) | **−0.34** |

- **ζ: no useful definite-sign R.**
  - W ≥ cN needs c ≤ μ_min(L). That falls from 5e-3 to 9e-4 to 1.6e-6 to the float floor by L = 2. So no uniform c > 0 exists.
  - W ≤ cN holds only with c ≥ μ_max(L), which grows from 0.86 to 15.6. That is only a non-uniform upper bound.
  - R itself is indefinite for every L > log 2.
- **Boundary against Connes–Consani 2021** (CHECKED as consistent).
  - At L = log 2 − 0.05, before any prime atom, W ≥ 0.0052·N and R ≥ 0.
  - Just past log 2, R turns indefinite and μ_min collapses.
  - This is the same boundary as 14C: positivity there is archimedean (CC 2021), and nothing global is added.
- **The fakes do not behave "exactly opposite".**
  - EPS and DH are *more* coercive than ζ at small L (μ_min 0.11 / 0.44 against 0.005).
  - μ_min < 0 occurs exactly where W itself has ghosts (EPS at L ≥ 3; F_θ at L ≥ 1). Because N > 0, sign(μ_min) = sign(λ_min W). **This is tautological, not a separating mechanism** (Golden Eye standard).

### 3c. Why no exact Gram identity W = ‖Tg‖² with a global, data-built T came out (structural, PROVED)
- On the line, N_F weighs |ĝ|² by |F(½+it)|², which is **quadratic in F**.
- W_F weighs |ĝ|² by the zero-counting measure dN_F(t) = (θ_F′(t)/π)dt + dS_F(t), with S_F = (1/π) arg F(½+it). This is **linear in log F**.
- Along the range of E_F the two behave oppositely: N_F > 0 there, while W_F = 0 (3a).
- So W restricted to range(E_F) is identically zero, and no T with ‖Tg‖² = W can be injective on range(E_F).
- The positive Hilbert space carrying W must be the **quotient** H_W = (admissible g) / closure(range E) with the W-norm. That is Connes' spectral realization (Selecta Math. 5, 1999), whose scaling spectrum on the cokernel is the zeros.
- **Positivity of W on that quotient ⇔ RH** (Weil / Connes). OPEN.

## Pair test (run first, as asked)

| model | E_F from data | own radical W_F(E_F H) = 0 | N_F = ‖E_F g‖² ≥ 0 | definite-sign R for this model? | what the fake lacks, and where it shows |
|---|---|---|---|---|---|
| ζ | yes | **yes** (2e-6) | yes | **no** for L > log 2 | – |
| EPS x²+5y² | yes (lattice) | **yes** (6e-6) | yes | no | Euler product. It shows only via W's ghosts (L ≥ 3), not via E or N |
| DH | yes | **yes** (2e-6) | yes | no | Euler product (vector mixed at 5). Same: E and N do not see it |
| F_0.2 | yes | **yes** (1e-4) | yes | no | unitary constituents. Ghosts are built in (cross term); E and N do not see it |
| spliced Beurling | yes (Beurling integers) | N.A. (no own explicit formula). With ζ's arch borrowed, W = −1.05e-3, **and E_B is ~ in ζ's radical (6e-6)** | N.A. (no ζ_B on the line) | – | the functional equation. Not separated by the radical test at this resolution |

**Pair-test verdict.**
- The required outcome was "identity fails or R flips sign exactly for the fakes". **It did not happen.**
- R has no definite sign even for ζ.
- The radical and the norm are present for every fake with an explicit formula.
- The fakes are separated only by W's own ghosts, and that separation is tautological.

## Verdicts per face
- **Faces 1 and 2 (theta and Euler): PROVED to be one object.** Tate's Θ on C_ℚ; Poisson for face 1, factorization of f₀ for face 2.
- **Face 3 (positive Hilbert space):** **not** the same object by any construction found here.
  - The global E gives the **null space** of W (PROVED, CHECKED), not its norm.
  - The global norm ‖Eg‖² is positive for every model, fakes included, and is not W (3b, 3c).
  - The only candidate left is the quotient by closure(range E) with the W-norm. Asking that it be positive is RH (Weil / Connes). **OPEN.**
- **Daniel's object, stated exactly:** (Tate Θ on C_ℚ, Connes quotient L²(C_ℚ) / closure E(S₀) with W).
  - Faces 1 and 2: PROVED.
  - The quotient is canonically defined, and W vanishes on what is divided out: PROVED.
  - Face 3 being a positive Hilbert space: OPEN, and equivalent to RH.

## Best exact claim (PROVED; CHECKED numerically)
**Claim.** Let F have Dirichlet coefficients a_F(n), a functional equation, and finite order. Let H ∈ C_c^∞(0, ∞) have Mellin transform vanishing at the poles of F. Then:
1. E_F H := Σ_n a_F(n) H(n·) is admissible.
2. W_F((E_F H) ∗ Θg) = 0 for every admissible g (the radical).
3. ‖E_F H‖² = (1/2π) ∫ |F(½+it)|² |Ĥ(½+it)|² dt.

**Consequence.** W_F descends to the quotient by closure(range E_F). "Face 3 = face 1 + face 2" holds exactly when W_ζ is positive on that quotient (Connes 1999 / Weil), i.e. RH.

Checked: radical to ≤ 1e-4 relative for ζ, EPS, DH, F_θ, against controls of 1e-3 to 1.

## Next step
- **Restrict W to the orthogonal complement of the window-truncated range of E**, i.e. quotient out the radical numerically.
- Test whether λ_min of W on that complement is bounded below **uniformly in L** for ζ.
- The near-null "very small eigenvalues" of W (Connes–Consani ζ-cycles) should be approximately range(E). If quotienting them out leaves a uniform gap, that is a non-tautological, uniform statement.
- Pair test it first:
  - for EPS, DH and F_θ, quotient by their *own* E_F ranges and check that their ghosts survive;
  - for BEUR″, quotient by E_B and check whether the borrowed-arch defect survives.

## Labels
- Faces 1 and 2 one object: PROVED (Tate, Riemann).
- Radical (3a): PROVED generically; CHECKED numerically.
- Plancherel normalization: CHECKED.
- Pencil and R-sign tables: NUMERICAL.
- sign(μ_min) = sign(λ_min W): PROVED, and tautological, so it is not counted as a separating control.
- 3c structural: PROVED (elementary).
- Face 3 positivity on the quotient: OPEN (⇔ RH).
- Nothing here closes RH, ARP or H_ψ.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
