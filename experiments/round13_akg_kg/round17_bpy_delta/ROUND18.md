# Round 18: "s IS Δ" in d = 1. Five zero-to-state candidates, with Epstein

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-19 (Toronto). Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN.
**RH: OPEN.** Zeros are never proof inputs. All states are built from theta series or Dirichlet coefficients. Off-line Epstein zeros are found by Newton from Round 17's printed label, and only to exhibit failures.

**Frame (Daniel, 2026-11-19).** We work in d = 1. The Pasterski–Shao line is Re Δ = d/2 = ½, so **s = Δ literally**. States live in H = L²(ℝ_{>0}, dx), primaries are ψ_Δ(x) = x^{−Δ}, and the shadow (Jf)(x) = x^{−1}f(1/x) sends Δ ↔ 1 − Δ. Requirement: **each zero ρ must yield a normalizable 1D primary/state.** This is Round 17's setting and GPPVerify's `XiZeroIsNormalizable1dState` (OPEN; equivalent to RH, PROVED in Lean). The d = 2 graviton version (Δ = 2s) is kept only as a cross-check in §6 and in ROUND17B.md.

Note: I chose candidates (a)–(d) myself because no explicit list reached this seat: the three ζ-objects named in the brief (BPY state, Tate vector, theta) plus Nyman–Beurling. (e) is the arithmetic Hardy-space candidate from the GPPVerify2 survey.

## Summary

| | Candidate (built from integers) | What a zero ρ gives | Normalizable when? | ζ-only? | Epstein |
|---|---|---|---|---|---|
| (a) | BPY state f_ζ (Round 17) | a primary invisible to f_ζ: (ψ_ρ, f_ζ) = ξ(ρ) = 0 | ψ_ρ δ-normalizable ⇔ Re ρ = ½. Not implied (GAP = RH) | no | same construction; invisible primary at Re Δ = 0.933, not normalizable |
| (b) | Tate local vector e^{−πx²} | nothing: Mellin = Γ_R(s)/2 ≠ 0 | — | no (local) | — |
| (c) | **integer lattice state** F_ζ(x) = Σ_{n≥1}e^{−πn²x²} − 1/(2x) ∈ L²(dx) (new) | invisible primary: ∫x^{ρ−1}F_ζ = ½Γ_R(ρ)ζ(ρ) = 0 | same as (a), GAP | no | F_EPS ∈ L², blind to non-normalizable Δ = 1 − ρ_E (overlap 1e−27) |
| (d) | Nyman–Beurling one-sided primary x^{ρ̄−1}1_{(0,1)} ⊥ floor-function span | an honest L² vector | **⇔ Re ρ > ½**, i.e. only for ghosts (PROVED) | uses floor functions of ℤ; general-Dirichlet-series version not checked here | — |
| (e) | arithmetic Hardy coherent states e_s(n) = n^{−s}, joint eigenvectors of A_p | e_ρ, A_p-eigenvalue p^{−ρ} | **ℓ² ⇔ Re ρ > ½** (ghosts only). At σ ↓ ½ *every* point of the line becomes a generalized vector, zero or not | **yes** (complete multiplicativity; EPS fails, CHECKED) | not an A_p-eigenvector |

**Bottom line (d = 1).** No candidate turns "ρ is a zero" into "ψ_ρ is δ-normalizable" without assuming RH. The candidates split into two classes:
- **Two-sided, positive or normalizable states (a), (c).** Every zero becomes a primary the state cannot see, but normalizability of that primary is a property of the label alone. The step "invisible ⇒ normalizable" is exactly RH, and it fails for Epstein, whose states pass every other check.
- **One-sided, Hardy-type states (d), (e).** Here zeros do produce honest normalizable vectors, but only off the line to the right, Re ρ > ½. These are the ghosts. RH says no such vectors exist (NB density, PROVED as an equivalence). On the line they degenerate to generalized vectors for every t, so the boundary limit cannot single out the zeros.

So in d = 1 the requirement "each zero yields a normalizable state" points the opposite way for the one-sided, ζ-specific constructions. The usable form is the **Nyman–Beurling / Hardy form: no zero yields a normalizable one-sided state.** The Euler product (A_p-joint eigenvectors, (e)) is the ζ-only input. The GAP is to make (e)'s multiplicativity force NB density, with no Weil form.

---

## 1. The d = 1 Pasterski–Shao facts used

- PS 2017 (arXiv:1705.01027) eq. (4.7) at general d: the massless KG product ∝ ∫dω ω^{2c−d−1+iν} is δ-normalizable only at c = d/2. **At d = 1 the principal line is Re Δ = ½.** The Mellin form of the wavefunction is (4.4), and completeness is (4.6) and (4.12).
- GPPVerify (main `7ebc350`): `radialDeltaNormalisable_one_iff_eq_conjShadow` (Δ = 1 − Δ̄), `normalizable_1d_implies_principal`, `mellin_half_plancherel`, `shadow_dilate`. The premise `XiZeroIsNormalizable1dState` (OPEN) ⇒ RH via `rh_of_xi_zero_is_normalizable_1d_state` (inventory in ROUND17.md §0).
- In d = 1, (i) "Re Δ = ½ ⇔ ψ_Δ is δ-normalizable" is a statement about labels (PROVED). (ii) RH needs every zero to be such a label. Every candidate is judged on whether it supplies (ii).

## 2. Candidates (a)–(c): two-sided states whose overlaps vanish at zeros

**(a) BPY state** (Round 17, PROVED; CHECKED 1e−41). f_ζ(x) = x^{−1/2}Φ(log x) ∈ L²(dx), positive, J f_ζ = f_ζ, and (ψ_Δ, f_ζ) = ξ(Δ) for every Δ ∈ ℂ. A zero gives a primary that f_ζ cannot see. Nothing forces that primary onto Re Δ = ½. The Epstein law is not positive (Φ_EPS(0) < 0, certified in Round 17), but F_θ is positive with every zero off the line, so positivity cannot be the missing input.

**(b) Tate local vector** (CHECKED). ∫x^{s−1}e^{−πx²}dx = Γ_R(s)/2, which never vanishes. The local Archimedean vector sees no zeros. Zeros appear only after the global sum over n ∈ ℤ, which is candidate (c). This places the zero-to-state problem in the global (integer) step, not the local one.

**(c) The integer lattice state ("plane waves of ℤ") — new, CHECKED** (`r17b_d1.py` → `r17b_d1.json`).
- ζ: F_ζ(x) = Σ_{n≥1}e^{−πn²x²} − 1/(2x). By Poisson summation F_ζ → −½ as x → 0, and xF_ζ → −½ as x → ∞. So **F_ζ ∈ L²(ℝ_{>0}, dx)**, and for 0 < Re s < 1 (Riemann's strip formula)
  ∫₀^∞ x^{s−1}F_ζ(x)dx = ½Γ_R(s)ζ(s),  i.e. (ψ_Δ, F_ζ) = ½Γ_R(1−Δ)ζ(1−Δ).
  The identity is CHECKED to 1e−26 at four points in the strip. ‖F_ζ‖² = 0.489371176033634 equals (1/2π)∫|½Γ_R(½+it)ζ(½+it)|²dt, Mellin–Plancherel on the d = 1 PS line, to relative 0.
- The Poisson pole term −1/(2x) is exactly what makes the raw theta sum normalizable. This is the additive half of Round 16 in Hilbert-space form.
- EPS: F_E(t) = Θ_Q(t) − 1 − 1/t, with Θ_Q(t) = θ(t/√5)θ(√5t). Then ∫t^{s−1}F_E = (√5/π)^sΓ(s)Z_Q(s) (CHECKED 1e−27), and ‖F_E‖² = 1.24552106630555 = Plancherel (relative 4e−24). The off-line zero ρ_E = 0.932969697 + 15.668249531i gives overlap 1.3e−27 against a scale of 9.3e−11 at the same height on the line. **A normalizable integer-built EPS state is blind to the non-normalizable primary Δ = 1 − ρ_E (Re Δ = 0.067).**
- Verdict: (c) is the cleanest d = 1 "integers as plane waves" object. It is normalizable, sees every zero, and is built only from ℤ and Poisson. But every step also holds for Epstein, so it cannot supply (ii).

## 3. Candidate (d): Nyman–Beurling, one-sided primaries

- With ρ_θ(x) = {θ/x} − θ{1/x} on (0, 1), ∫₀¹ρ_θ(x)x^{s−1}dx = −(ζ(s)/s)(θ^s − θ) for 0 < Re s < 1. So a zero ρ makes e_ρ(x) := x^{ρ̄−1}1_{(0,1)} orthogonal to every ρ_θ.
- ‖e_ρ‖² = ∫₀¹x^{2Re ρ−2}dx = 1/(2Re ρ − 1), **finite iff Re ρ > ½** (PROVED). Off-line zeros come in pairs ρ, 1 − ρ̄, so any off-line zero gives such a vector. For the Epstein zero, 1/(2·0.93297 − 1) = 1.1548 (`r17b_d1.json`). That is shown only to illustrate the norm; the Epstein floor-function analogue was not built.
- Nyman (1950), Beurling (1955), Báez-Duarte (2003) (PROVED): RH ⇔ 1_{(0,1)} lies in the closed span of the ρ_θ (or, with Báez-Duarte, of the integer dilates ρ_{1/k}).
- In PS language: off-line zeros are exactly the L²-normalizable truncated primaries ("bound states") orthogonal to the integer floor span. On the line the truncated primary is only a generalized vector. So the d = 1 requirement "zero ⇒ normalizable" holds here **only for ghosts**. The RH-equivalent statement is **"no zero yields a normalizable one-sided state"**, which is density, a Weil-free form.

## 4. Candidate (e): arithmetic Hardy space, coherent states and σ ↓ ½ (already in GPPVerify2)

Survey clone `/workspace/gppv2-survey`, main `d6cbaa2`, read-only:
- `GppVerify/RiemannHypothesis/ArithmeticHardySpace.lean` (PROVED; note the folder is `RiemannHypothesis/`, not `CelestialHolography/`). `coh z n = n^{−½−z}`, i.e. e_s with s = ½ + z. In ℓ² only for Re z > 0 (`summable_norm_sq`, l.72). Gram `coherent_inner` (l.117): ⟨e_w, e_z⟩ = ζ(1 + z + w̄). Joint eigenvectors: `annihilate_coherent` (l.195), (A_p f)(n) = f(pn), eigenvalue p^{−½−z} = p^{−s}. Unitary dilation `dilate_coherent` (l.208): z ↦ z + it.
- `docs/critical_zeta_coherent_state_boundary_2026-09-13.md`: K_σ(t,u) = ζ(2σ + i(u−t))/ζ(2σ) → δ_{t,u} as σ ↓ ½ (an orthogonality catastrophe; §3), and the vacuum fidelity 1/ζ(2σ) → 0 (§4–5).
- `GppVerify/CelestialHolography/ArithmeticTemperedPrincipalSeries.lean`: `CenteredTempered ρ` (l.38) ⇔ Re ρ = ½ (`centeredTempered_iff_principal`, l.46). `all_strip_zeros_principal_of_tempered_bridge` (l.75) takes **`hbridge`: zero ⇒ CenteredTempered** as a hypothesis. That hypothesis is OPEN and is RH restated.
- The prime-local OS and factor-map files (`ArithmeticPrimeLocalOS.lean`, `ArithmeticPrimeFactorMap.lean`, `ArithmeticOSFactorization.lean`) take `p : ℝ` with `1 ≤ p` and never require primality. They are positivity facts about any base, not ζ-specific.

**Can a σ ↓ ½ boundary limit make each zero a δ-normalizable generalized eigenvector? No, not without RH.** Four reasons (PROVED, with CHECKED illustrations in `r18_e_hardy.py` → `r18_e_hardy.json`):
1. **The Gram matrix never sees a nontrivial zero.** For e_s, e_w in ℓ², or at the boundary, ⟨e_w, e_s⟩ = ζ(s + w̄) with Re(s + w̄) ≥ 1, the zero-free region (Hadamard–de la Vallée Poussin). Reaching a zero ρ needs Re s + Re w = Re ρ < 1, so at least one vector has Re < ½. That vector is neither in ℓ² nor a boundary vector.
2. **Zero-ness is a pairing with the non-normalizable vacuum column.** ζ(ρ) = ⟨𝟙, e_ρ⟩_reg with 𝟙 = (1, 1, 1, …) ∉ ℓ². The A_p-eigenvalue p^{−ρ} of e_ρ is the same whether or not ζ(ρ) = 0. The eigenvector structure exists for every ρ, so it cannot by itself mark the zeros.
3. **The boundary limit is indiscriminate.** As σ ↓ ½ every e_{½+it} becomes a generalized vector, zero or not, and no off-line e_ρ does. Unnormalized, ζ(1 + ε + iτ) → πδ(τ) − i·PV(1/τ) + regular part. The CHECKED integral ∫_{|τ|<1} Re ζ(1+ε+iτ)dτ = 4.1145, 4.2808, 4.2974 for ε = 0.1, 0.01, 0.001 approaches π + 2γ ≈ 4.296. So the boundary vectors are δ-normalizable only up to a Hilbert-transform (PV) term: a one-sided Hardy structure, not an orthonormal δ-basis. Normalized by ζ(2σ), K_σ(t, t+1) = 0.4799, 0.2128, 0.02187, 0.002189, 0.0002189 at σ = 0.75, 0.6, 0.51, 0.501, 0.5001 (the catastrophe).
4. **Normalizability points the wrong way, as in (d).** e_ρ ∈ ℓ² ⇔ Re ρ > ½. Indeed e_ρ(n) is the sampling of x^{−ρ}1_{x≥1} at the integers, the discrete twin of NB's one-sided primary. A zero gives an honest ℓ² eigenvector only if it is a ghost.

**ζ-specificity (CHECKED).** For x² + 5y², v_s(n) = r_Q(n)n^{−s} is not an A_p-eigenvector:
- p = 2 and p = 3: r(pn) = 0 on the support of r, yet r(pn) ≠ 0 for n = 2, 3, 7, 8, 10, where r(n) = 0. The genus swap moves between the two forms.
- p = 29: r(29n)/r(n) takes the values {1.5, 2}.

So (e) genuinely uses the Euler product, which (a)–(c) do not.

**Comparison of (d) and (e).** Both are one-sided (Hardy) pictures in which off-line zeros are the normalizable "bound states":
- (d) has the **PROVED** RH-equivalence (density) but no multiplicative structure in the statement.
- (e) has the **Euler product** (joint A_p eigenvectors) but no density theorem.
- **GAP-18:** prove NB density for the integer dilates using (e)'s joint-eigenvector structure, i.e. show that the completely multiplicative coherent family leaves no room for an ℓ²/L² bound state orthogonal to the floor span. That is a Weil-free formulation in which the ζ-only input (complete multiplicativity, equivalently c(n) = Λ(n) ≥ 0 from Round 16) actually enters. It is RH-equivalent, so OPEN.

## 5. Epstein verdict (d = 1)

| Step | ζ | EPS |
|---|---|---|
| (a)/(c) normalizable integer-built state, Plancherel on Re Δ = ½ | ✓ | ✓ (F_E: 1.2455 = 1.2455) |
| shadow Δ ↔ 1 − Δ | ✓ | ✓ |
| overlap vanishes at zeros | ✓ | ✓ |
| invisible primary is normalizable | **OPEN (= RH)** | ✗ (Re Δ = 0.067; overlap 1e−27) |
| (e) joint A_p eigenvectors (Euler product) | ✓ (Lean) | ✗ (CHECKED) |

The only discriminating structure among the candidates is (e)'s complete multiplicativity, and it enters on the one-sided (Hardy/NB) side.

## 6. Cross-check: the graviton frame (d = 2, Δ = 2s), see ROUND17B.md

- The unitary map (Ub)(ω) = √2 b(ω²) from L²(ℝ₊, dy) to L²(ℝ₊, ω dω) sends the d = 1 primary of weight s to the graviton energy profile of weight Δ = 2s. The Lorentz-invariant measure d³p/2ω is flat in ω². Every d = 1 statement above transports verbatim (PROVED).
- The graviton adds kinematics only:
  - it forces the factor 2;
  - it rules out the Planck reading s = Δ_{d=2}, which would put the critical line at Re Δ = ½ (CHECKED);
  - its soft poles at Δ = 2, 0 come from the same Poisson term −1/(2x) as in candidate (c).
- It adds no constraint on zeros. The d = 2 GAP is GPPVerify's `GravitonZeroIdentification` (OPEN), the same premise.

## Files (this folder)
- `r17b_d1.py` → `r17b_d1.json`, `r17b_d1.log`: candidates (b), (c) (ζ and EPS: strip identity, L² norm = Plancherel), the EPS invisible primary, and the (d) norm. 25 digits, about 90 s.
- `r18_e_hardy.py` → `r18_e_hardy.json`: candidate (e). Orthogonality catastrophe, the boundary δ + PV check, and the Epstein non-eigenvector check.
- ROUND17B.md and `r17b_graviton.py`, `r17b_spin2_sympy.py`: the d = 2 cross-check.

## References
Pasterski–Shao, arXiv:1705.01027, eqs. (4.4), (4.6)–(4.12) · Pasterski–Shao–Strominger, arXiv:1701.00049 · Riemann (1859) strip formula · Tate (1950) · Nyman (1950), Beurling, PNAS 41 (1955), Báez-Duarte, Rend. Lincei (2003) · Biane–Pitman–Yor, Bull. AMS 38 (2001) · GPPVerify main `7ebc350` · GPPVerify2 main `d6cbaa2` (ArithmeticHardySpace.lean, ArithmeticTemperedPrincipalSeries.lean, ArithmeticPrimeLocalOS / PrimeFactorMap / OSFactorization, docs/critical_zeta_coherent_state_boundary_2026-09-13.md) · Rounds 16 and 17 (`9e40820`, `3383dd6`).

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
