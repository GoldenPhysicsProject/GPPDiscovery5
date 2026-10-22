# Queen of Mathematics, session 1 (2026-10-22): what a non-circular positivity mechanism for the Weil form has to be, and one uniform-in-L statement that kills DH. RH OPEN.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

**Labels.** PROVED means a complete argument here or a citation. CERTIFIED would mean Arb; nothing in this note is CERTIFIED. NUMERICAL means float64 or mpmath, not interval arithmetic. CONJECTURE means unproved. All citations previously marked *(recall; verify)* were checked against the source on 2026-10-22 (fix commit). Unverifiable ones were dropped. The remaining unverified-from-memory items are general background (Weil, Yoshida, Connes, Suzuki, Burnol, Lagarias) and carry no claim here beyond what the project already cites.

**L convention (binding throughout).** L is the **half-width** of supp f in log scale. The window is [−L, L], the full support width is 2L, and the prime atoms seen are log n < 2L.
- Team certificates are quoted by FULL width: "support 2.8" means L = 1.4, and "support 2.84" means L = 1.42 (S72's L = 71/50).
- S75's DH onset L* ∈ (1.70, 1.7105] is a half-width (support ≈ 3.421).
- Golden Eye's sign-change formula 2Lγ/π ≈ 93 at L = 1.71 is also half-width (§5, sign-change note).
- I have not seen Golden Eye's L ≈ 1.61 derivation, so its convention is unconfirmed. 1.61 is not 1.7105/2 or 2·1.7105, so a factor-2 convention mix does not explain the gap.

**Files** (all in `team/queen/`):

| Script | Output | What it does |
|---|---|---|
| `pf_ground_state.py` | `.log` | Spectrum, nodes and pole overlap of the non-pole operator A − P from the s91 dps-90 cache. Run it from `team/queen/`. It needs `../../s91/forms_cache.pkl`. |
| `pf_nodes.py` | `.log` | Sign changes of the first three eigenvectors per parity. |
| `dh_coeff.py` | `.log` / `.json` | Lévy–Khintchine check; the DH log-derivative coefficients; Riesz-mean sups. |
| `dh_riesz2.py` | `.log` | DH Riesz mean at the zero 0.8085 + 85.70i compared with the Perron residue. |
| `dh_sinh.py` | `.log` / `.json` | An exact, Galerkin-free Rayleigh quotient of Q_DH and Q_ζ on the off-line-pair vector. |

---

## 0. Summary

1. **The decomposition (PROVED, classical).** On L²[−L, L] the window Weil form is exactly
   **Q_L = Pole + E_M − κ_L·Id**.
   - E_M is the Beurling–Deny Dirichlet form of a symmetric Lévy process on the log-line. Its Lévy measure M is the archimedean density ν(u) = e^{u/2}/sinh u plus atoms 2Λ(n)n^{−1/2} at u = log n.
   - κ_L is a constant, ≈ 4e^L.
   
   So a genuine Markov semigroup built from the Euler product and Γ exists *unconditionally*. That is the non-circular positive object the brief asked for. Its covariance is E_M, though, not W. W differs from it by the indefinite rank-2 pole term and the huge negative constant −κ_L.
2. **Perron–Frobenius (PROVED given Λ ≥ 0).**
   - The non-pole operator R_L = A − P = H_M − κ_L has a simple, strictly positive, even ground state φ_L.
   - The pole's positive functional cosh(x/2) > 0 therefore always overlaps φ_L.
   - NUMERICAL (ζ, L = 0.5–1.5, N = 24/32): R_L has **exactly one** negative eigenvalue, in the even sector, and none in the odd sector. φ_L is nodeless, and |cos(φ_L, pole)| = 0.995–0.9997.
   - This is the structural explanation of S84–S90's "nodeless pole bump" and of S98's nodeless soft mode.
3. **The S98 constraint, made exact.** If v ≥ 0, then g_v = v * ṽ ≥ 0 pointwise, so every positive atom *lowers* Q(v). On the PF ground state the Euler product's positivity works against positivity.
   - Window-RH is decided by the **second** even level λ₂(R_e) and the odd sector. NUMERICAL: λ₂(R_e)/‖R‖ runs from 2.4e−3 down to 7e−27. Its eigenvector has exactly one sign change on (0, L).
   - So any forcing property has to act on **sign-changing** vectors, through cancellation in Σ Λ(n)n^{−1/2} g_v(log n). That is a phase or correlation property of the primes, not their sign.
4. **The specification.** A non-circular mechanism needs properties P1–P7 (§3). Two pairs of examples pin down what it must use:
   - **{F_θ, G_θ′}**: F_θ has a positive (Markov) Lévy measure and fails; G_θ′ has a *signed* Lévy measure and passes. So "Λ ≥ 0" is neither necessary nor sufficient. The correct local property is **unitarity of the Satake parameters** (Herglotz positivity of k ↦ Λ_D(p^k)/log p), not sign.
   - **{DH, DMV}**: DH has the functional equation (FE) but no Euler product; DMV has an Euler product but no FE or lattice. Both have off-line zeros. So the mechanism must use the Euler product **and** the additive (Poisson / FE) structure *of the same lattice* ℤ ⊂ ℝ (ℚ ⊂ 𝔸) jointly. Neither alone suffices.
   - On a fixed window this joint input is invisible (hybrid lemma). It must enter **uniformly in L**.
   - I know of no construction meeting P1–P7. OPINION: the closest known framework is the Connes adèle-class trace formula, which appears to have P1, P3 and P5. Its missing piece is the uniform cutoff estimate, which is open.
5. **One uniform-in-L statement (§5).**
   - **U:** λ_min(Q_L) ≥ −(5.38 + 2(sinh L − L) + 4.16e^L) for all L. This is PROVED for ζ here. It uses Chebyshev ψ(x) < 1.0389x, i.e. the Euler product through absolute convergence.
   - It is fatal to DH *modulo a Landau-type lemma* (CONJECTURE, provable in principle; sketch in §5). DH has zeros with Re s > 1 (Davenport–Heilbronn 1936), and an off-line pair at real part ½ + β forces λ_min ≲ −e^{2βL}/(2β).
   - NUMERICAL support for the pair mechanism: on the explicit vector e^{−iγ₀x} sinh(β₀x) at DH's zero ρ₀ = 0.8085 + 85.70i, Q_DH/‖f‖² = −0.24, −1.05, −3.09, −10.5, −24.4, −53.0 at L = 2, 2.5, 3, 4, 5, 6. ζ on the same vector stays at +0.9 to +3.8.
   - U is **not** RH-equivalent; it is weak (Euler-product convergence strength). Strengthening it to λ_min ≥ −Ce^{δL} with δ < 1 is equivalent (heuristically) to a zero-free half-plane, which is open. δ → 0 is RH.
   - It separates ζ from DH but not from F_θ, H₋ₜ or DMV.

---

## 1. Setup and the exact Lévy–Khintchine / Beurling–Deny decomposition

**Explicit formula.** For g ∈ C_c^∞(ℝ), with ĝ(z) = ∫ g(x)e^{izx} dx (Weil; Bombieri 2000 normalisation):

  Σ_ρ ĝ(γ_ρ) = ĝ(i/2) + ĝ(−i/2) − Σ_n Λ(n)n^{−1/2}[g(log n) + g(−log n)] + (1/2π)∫ ĝ(t) Ω(t) dt,

where ρ = ½ + iγ_ρ and Ω(t) = Re ψ(¼ + it/2) − log π.

**The window form.** For f ∈ L²[−L, L], let g = f * f̃, so g(u) = ⟨f, τ_u f⟩ and supp g ⊂ [−2L, 2L]. Then Q_L(f) := W(g).

**The pole term** is P(f) = 2Re(a b̄), where a = ∫ f e^{−x/2} and b = ∫ f e^{x/2}. Equivalently P(f) = ½|⟨f, 2cosh(x/2)⟩|² − ½|⟨f, 2sinh(x/2)⟩|². This has rank 2 and signature (1, 1): it is positive on the even sector and negative on the odd sector.

**(1.1) Lévy–Khintchine (PROVED: Gauss's digamma integral; NUMERICAL check `dh_coeff.log` to 3e−15).**

  Ω(t) = Ω(0) + ∫₀^∞ (1 − cos tu) ν(u) du, where ν(u) = e^{u/2}/sinh u > 0 and Ω(0) = ψ(¼) − log π = −5.37218.

**(1.2) Window decomposition (PROVED; elementary from (1.1)).** Two facts give it:
- (1/2π)∫|f̂|²(1 − cos tu) dt = g(0) − Re g(u) = ½‖f − τ_u f‖².
- g(u) = 0 for |u| ≥ 2L.

With them,

  Q_L(f) = P(f) + E_M(f) − κ_L‖f‖², where E_M(f) = ½∫_{(0,∞)} ‖f − τ_u f‖² dM(u),

  M = ν(u)du + Σ_n 2Λ(n)n^{−1/2} δ_{log n}, and κ_L = −Ω(0) + Σ_{n<e^{2L}} 2Λ(n)n^{−1/2} (≈ 4e^L).

(The jump of size u = 2L or more contributes exactly 2‖f‖², and that is absorbed into κ_L. The only point is that only atoms with log n < 2L see the window.)

**Interpretation.**
- Restricted to f supported in [−L, L], E_M is the regular Dirichlet form of the symmetric Lévy process with Lévy measure M_sym, killed on leaving [−L, L].
- Equivalently E_M is the jump form on the window plus the killing potential k(x) = M_sym(ℝ ∖ ([−L, L] − x)).
- H_M is its generator, with symbol ψ_M(t) = ∫(1 − cos tu) dM(u) ≥ 0, and e^{−sH_M} is sub-Markovian (Beurling–Deny; Fukushima–Oshima–Takeda Thm 1.4.1).
- M is **not** a Lévy measure on all of ℝ, because Σ Λ(n)n^{−1/2} = ∞. Only the truncations are. This is the probabilistic face of "κ_L ≈ 4e^L grows with L".

**(1.3) The same for other data.**
- DH uses c_DH(n) from −DH′/DH, which is **signed** and not supported on prime powers: 74% of its support is off prime powers and 50% of the signs are negative (NUMERICAL, n ≤ 4·10⁵, `dh_coeff.log`). DH uses ν_odd(u) = e^{−u/2}/sinh u, Ω_DH(0) = ψ(¾) − log(π/5), and has **no pole**.
- F_θ: the atoms are Λ(n)n^{−1/2}(n^θ + n^{−θ}) > 0, and there are two pole pairs.
- G_θ′: the atoms are 2Λ(n)n^{−1/2}cos(θ′ log n), which are **signed**.

## 2. Perron–Frobenius structure of the non-pole operator; the S98 constraint

Let R_L := Q_L − P = E_M − κ_L (in the s91 cache, R = A − P).

**Lemma 2.1 (PROVED, standard).** Assume M ≥ 0 (i.e. Λ ≥ 0 and ν ≥ 0). Then R_L, as a closed form on L²[−L, L], has compact resolvent, and its lowest eigenvalue is simple, with a strictly positive eigenfunction φ_L. Since M is symmetric, φ_L is even.

*Proof.*
1. ψ_ν(t) ~ log|t| → ∞, so the form domain embeds compactly in L²[−L, L] (Kolmogorov–Riesz).
2. The first Beurling–Deny criterion holds: E(|f|) ≤ E(f), because ||f(x)| − |f(y)|| ≤ |f(x) − f(y)|. So e^{−sR} is positivity preserving.
3. ν > 0 on (0, ∞) makes the jump kernel strictly positive, so the semigroup is positivity improving.
4. Then apply Reed–Simon IV, Thm XIII.44. ∎

**Corollary 2.2 (PROVED).** ⟨φ_L, cosh(x/2)⟩ > 0. So the pole's single positive direction always "sees" the bound state.

**Lemma 2.3 (PROVED; rank-one interlacing).**
- Q_e = R_e + pp^T and Q_o = R_o − qq^T, where p and q are the pole vectors.
- So Q_L ⪰ 0 implies n₋(R_e) ≤ 1 and R_o ⪰ qq^T ⪰ 0.
- In particular, window positivity forces **λ₂(R_e) ≥ 0** and **λ₁(R_o) ≥ 0**. That is, the second even level and the first odd level of the arithmetic Lévy generator H_M must both sit above the threshold κ_L.

**NUMERICAL (dps 90, s91 cache, `pf_ground_state.log`, `pf_nodes.log`).**

| Form | n₋(R_e) | n₋(R_o) | Nodes of φ_L on (0, L) | \|cos(φ_L, p)\| | λ₂(R_e)/‖R‖ |
|---|---|---|---|---|---|
| ζ, L 0.5–1.5 (N 24; N 32 at 0.6, 1.0, 1.2, 1.4) | **1** at every L | **0** | 0 | 0.9949–0.9997 | 2.4e−3, 1.4e−5, 1.2e−8, 1.4e−12, 4e−17, 1.1e−22 (L 1.0, N 32), 2.4e−27 (L 1.2, N 32), 5e−31 (L 1.4, N 32) |
| F_0.1 (2 pole pairs) | 2, 3, 6, 7 at L 0.6, 0.8, 1.0, 1.2 | 0, 2, 3, 5 | 0 | ≈ 0.997 | negative |
| DH (no pole) | 0 | 0 | 0 | n/a | PD |

Further observations from the same runs:
- The second even eigenvector of ζ's R has exactly one sign change on (0, L) at every L.
- F_θ is Markov (M_θ ≥ 0), so its ground state is nodeless too. From L 0.8 it fails **by count**: n₋(R_e) = 3 > 2 positive pole directions.
- DH's ground state is also nodeless (PF is not guaranteed for signed M, but it holds here), so nodelessness is **not** a separator.
- For L ≥ 1, N 24 is unconverged at the λ₂ level. The values are upper bounds, as in S91.

**N-convergence caveat (lead directive, binding).** At L ≥ 1 the N = 24 basis is NOT converged at the λ₂ level. Where both exist, N 24 and N 32 differ by orders of magnitude:

| L | λ₂(R_e)/‖R‖ at N 24 | at N 32 |
|---|---|---|
| 1.0 | 1.6e−19 | 1.1e−22 |
| 1.2 | 8.7e−23 | 2.4e−27 |
| 1.4 | 3.2e−24 | 5.0e−31 |

So **no λ₂ or λ_min value at L ≥ 1 in this note is a converged number.** They are NUMERICAL upper bounds only.

What *is* stable between N 24 and N 32 at L = 0.6, 1.0, 1.2 and 1.4:
- the counts n₋(R_e) = 1 and n₋(R_o) = 0;
- the nodelessness of φ_L;
- |cos(φ_L, p)| (to 4 digits);
- the one sign change of the second even eigenvector.

These are the claims I make (NUMERICAL). The converged value is at L = 0.6, where λ₂(R_e)/‖R‖ = 1.42e−5 (N 24) against 1.26e−5 (N 32), and the counts agree. A converged λ₂(L) curve, at N = 24/32/40/48 until 2 digits are stable, is next step (c).

**S98, made exact.**
- If v ≥ 0, then g_v(u) = ∫ v(x)v(x − u) dx ≥ 0, so each atom contributes −2a_n n^{−1/2} g_v(log n) ≤ 0 when a_n > 0. This is S98's H97a observation, and it is automatic once the soft mode is nodeless; Lemma 2.1 explains why it is nodeless.
- So positivity of prime weights is the *wrong* sign on the PF direction. The pole has to rescue that direction, and does (Cor 2.2).
- What remains, the RH-tight part, is λ₂(R_e) and λ₁(R_o). Both live on sign-changing vectors (Courant: their eigenvectors are orthogonal to φ_L > 0, so they must change sign).
- **Consequence for the spec:** the property that does the work cannot be "Λ ≥ 0" in any monotone way. It has to control the *oscillating* sums Σ Λ(n)n^{−1/2} g_v(log n) for sign-changing v, which amounts to cancellation in ψ(x) − x. G_θ′ (§4) confirms independently that sign is not the right invariant.

## 3. The specification: what a non-circular mechanism must be

**Definition.** A *positivity mechanism* is a construction Φ that takes arithmetic data D and outputs a Hilbert space 𝓗_D and a linear map J_D : C_c^∞(ℝ) → 𝓗_D. The data D are: local factors {L_v(s)} at all places v ≤ ∞, the global set (lattice or monoid) on which the Dirichlet series lives, and its functional equation if it has one. Φ must satisfy:

- **P1. Unconditional existence.** 𝓗_D is a genuine Hilbert space. Its positivity is manifest (a measure, a state, a Markov semigroup, a Fock space), proved without reference to zeros.
- **P2. Covariance identity.** Q_D(f) = ‖J_D f‖² + R_D(f) for all f ∈ C_c^∞(ℝ), with R_D ≥ 0 manifestly (R_D = 0 is allowed). Lemma 2.3 shows how the pole can be included: the term must produce the odd negative pole direction q and beat it.
- **P3. Non-circularity.** Φ is defined by a rule applied to D. It is not a GNS or OS reconstruction of W_D, and it uses neither zeros nor window positivity.
  - *Why this is needed (PROVED; Bochner–Schwartz / GNS):* a pair (𝓗, J) with ‖Jf‖² = Q(f) on C_c^∞ exists iff W is a positive-definite distribution, which happens iff RH holds. Existence alone is RH; only the rule is content.
  - Log-time OS (x ↔ 1/x) is false for ζ. Heat-time OS is equivalent to RH. The S93–S97 Schur, Kreĭn and doubling reductions are all P2 with J = GNS of W, so they fail P3.
- **P4. Local unitarity (Satake/Herglotz).** For each p, the sequence k ↦ Λ_D(p^{|k|})/log p, with value deg D at k = 0, must be positive definite on ℤ. Equivalently, the local Lévy measure at p is the push-forward of a positive measure on the unit circle (Satake parameters with |α_j(p)| = 1). The analogous condition at ∞ is that the Γ-shifts are tempered, so that ν_D ≥ 0 comes from Re μ_j ≥ 0 unitary data.
  - *Why the sign of Λ is the wrong invariant (PROVED by identity).* Q_{G_θ′}(f) = Q_ζ(e^{iθ′x}f) + Q_ζ(e^{−iθ′x}f) (unitary twist), whereas Q_{F_θ}(f) = 2Re B_ζ(e^{θx}f, e^{−θx}f), where B_ζ is the sesquilinear Weil form. That second expression is a polarisation between two non-unitary twists, and so it is indefinite.
  - Both identities are linearity of the explicit formula: G_θ′ = ξ(s + iθ′)ξ(s − iθ′) and F_θ = ξ(s − θ)ξ(s + θ), with zeros at γ ∓ θ′ and γ ∓ iθ respectively. The rep-theory seat's twist identity is the same statement.
  - So Φ must be functorial under unitary twists, Φ(D ⊗ χ) = Φ(D) ⊗ χ, and undefined (or indefinite) for non-unitary ones.
- **P5. Joint multiplicative–additive structure.** Φ must use the Euler product **and** the Poisson/FE structure of the *same* lattice, i.e. Tate's thesis: ζ is the Mellin transform of a theta series over ℚ^× acting on 𝔸, and the FE is Poisson summation for ℚ ⊂ 𝔸.
  - DH has the FE (each L(s, χ mod 5) is Tate) but the combination has no Euler product.
  - DMV has an Euler product over Beurling primes, but its integers are not a lattice: N(x) = Ax + O(x^θ) with θ > ½, and there is no FE.
  - Since both have off-line zeros, neither input alone can be sufficient.
- **P6. Uniformity in L.** The identity in P2 must hold with no L-dependent cutoff left over, or the cutoff error must be controlled uniformly in L.
  - *Why (PROVED; Goddess hybrid lemma given DMV, Math. Ann. 334 (2006)):* for every X there is a Beurling system with Λ ≥ 0, an Euler product, de la Vallée Poussin PNT and off-line zeros whose window forms (in ζ's arch and pole slot) **equal** ζ's for every 2L < log X.
- **P7. Correct pole and bound-state accounting.** The construction must reproduce, not assume, the facts that:
  - R_L has exactly one negative even level, absorbed by the pole (Lemma 2.3); and
  - the RH content sits on the sign-changing levels λ₂(R_e) and λ₁(R_o) (§2).
  
  Any mechanism whose positivity comes from monotonicity in Λ ≥ 0 is refuted by S98.

**Status of candidate shapes.**

| Shape | P1 | P2 | P3 | P4 | P5 | P6 | Verdict |
|---|---|---|---|---|---|---|---|
| Gaussian field / OS / KMS with covariance W | iff RH | ✓ | ✗ | – | – | – | Circular (target, not mechanism) |
| Lévy/Markov semigroup H_M (§1, this note) | ✓ | ✗ (covariance E_M; the remainder P − κ_L is indefinite) | ✓ | ✗ (Markov ⇔ Λ ≥ 0: F_θ is Markov and fails, G_θ′ is not and passes) | ✗ | n/a | Right object for the ground state (PF); says nothing on sign-changing levels. **DEAD END** as the mechanism. |
| Bost–Connes KMS_β; Riemann-zeta distribution (ζ(β+it)/ζ(β) infinitely divisible, Lévy measure Σ p^{−kβ}/k δ_{k log p}) | ✓ for β > 1 only | ✗ (PD functions in t at β > 1, not W) | ✓ | ✓ | ✗ (DMV has the same KMS structure) | – | Positivity lives on β > 1. Crossing the pole at β = 1 (the phase transition) is exactly the problem. DH fails already: DH(β+it)/DH(β) is not a characteristic function. **DEAD END** as stated. |
| Connes adèle class space / trace formula (Connes 1999; CCM 2007; CC 2021–23 prolates) | ✓ | ✓ up to a cutoff error | ✓ | ✓ | ✓ (uses ℚ ⊂ 𝔸) | ✗: RH ⇔ uniform control of the cutoff error (Connes 1999) | OPINION: **the closest known framework**; it appears to have P1, P3, P4 and P5. The open step is P6. S91–S95 numerics are its finite shadow. |
| de Branges / Kreĭn canonical systems | needs PSD Hamiltonian | ✓ | ✗ (S97) | – | – | – | Restates nested positivity (S96–S97). The de Branges positivity condition is false for ζ (Conrey–Li 2000). |

**Honest conclusion.** P4 + P5 + P6 together are, in substance, "ζ is an automorphic L-function of GL(1)/ℚ, used uniformly in the cutoff". No theorem is known that turns those properties into positivity. That would be GRH for the class. What this note adds is a sharpening of the *shape* the mechanism must have:
- it is not a sign condition (S98 + G_θ′);
- it acts on sign-changing levels λ₂(R_e) and λ₁(R_o) of an explicit Markov generator; and
- it must couple multiplicative and additive structure uniformly in L.

## 4. The fakes: which property breaks

| Fake | P4 local unitarity | P5 joint Euler + lattice | P6 | Where it fails in this note |
|---|---|---|---|---|
| **DH** | ✗: no local factors; the c(n) are signed and 74% supported off prime powers (NUMERICAL) | ✗: FE yes, Euler product no | – | Zeros in σ > 1: Davenport–Heilbronn, J. London Math. Soc. 11 (1936). Verified via two sources. Karatsuba, Math. USSR Izv. 36 (1991) states DH has ≫ T zeros in Re s > 1, 0 < Im s < T, citing DH 1936 and Titchmarsh §10.25. Saias–Weingartner, Acta Arith. 140 (2009) 335–344, doi:10.4064/aa140-4-4, Corollary: DH is a combination of L(s, χ) and L(s, χ̄) mod 5, so it is not P·L_χ and vanishes in Re s > 1. The supremum of the real parts is σ_DH ≈ 1.120362 (Bombieri–Ghosh, Russian Math. Surveys 66:2 (2011), as quoted by Odžak–Smajlović, LMS J. Comput. Math. 2016, doi:10.1112/s1461157016000115). So U (§5) fails, by Lemma 5.2 (now PROVED). Off-line pair 0.8085 + 85.70i: explicit negative vector for L ≥ 2 (NUMERICAL). |
| **F_θ** | ✗: Satake p^{±θ} is non-unitary; Γ((s ± θ)/2) is non-tempered | Euler product ✓, FE ✓ | – | Polarisation identity (P4); R_e has more bound states than pole directions from L 0.8 (NUMERICAL). Passes U: its zeros have Re < 1. |
| **H₋ₜ** (de Bruijn–Newman, t < Λ_dBN) | ✗: no local factors | ✗: no Euler product; the Weil distribution is ζ's plus an entire (non-discrete) correction (S93–S94 "dd" piece) | – | Fails P5 at input level. Passes U. The S93 zero-swap model has finitely many off-line zeros, so no density statement sees it either; only full positivity does. |
| **DMV** (in ζ's arch and pole slot) | ✓: Λ_P ≥ 0, Satake 1 | ✗: no lattice, no FE, N(x) − Ax ≍ x^θ with θ > ½ | equal to ζ on every window 2L < log X (hybrid lemma) | Passes U and every fixed-L test. Excluding it needs **additive** input: the FE/Poisson for the same integers, or at least integer regularity β < ½ (open whether that suffices). Broucke–Debruyne–Révész, *Some examples of well-behaved Beurling number systems*, Trans. AMS, doi:10.1090/tran/9274, construct [α, β]-systems for α ∈ [0, 1), β ∈ [½, 1), and certain β < ½ families only assuming RH. (The earlier Hilberdink citation is dropped as unverified.). |
| **G_θ′** (must pass) | ✓: Satake p^{±iθ′} is unitary, but Λ_G is **signed** | ✓: Eisenstein on GL(2) | ✓ | Passes by the unitary twist identity. This is the reason "Λ ≥ 0 / Markov" cannot be the invariant. |

## 5. The uniform-in-L statement

**Statement U (proposal).** For all L > 0,

  λ_min(Q_L) ≥ −B(L), where B(L) := 5.373 + 2(sinh L − L) + 4.156 e^L.

Equivalently, the *negativity exponent* η := limsup_{L→∞} L^{−1} log⁺(−λ_min(Q_L)) is at most 1.

**Theorem 5.1 (PROVED for ζ).** U holds for ζ.

*Proof.* Take ‖f‖ = 1, supp f ⊂ [−L, L].
1. **Arch.** By (1.1), Ω(t) ≥ Ω(0), so arch ≥ Ω(0) = −5.3722.
2. **Pole.** P(f) ≥ −½‖2sinh(x/2)‖² = −2(sinh L − L).
3. **Primes.** |g(u)| ≤ 1, so |prime part| ≤ 2Σ_{n<x} Λ(n)n^{−1/2} with x = e^{2L}.
4. By partial summation, Σ_{n≤x} Λ(n)n^{−1/2} = ψ(x)x^{−1/2} + ½∫₁^x ψ(t)t^{−3/2} dt ≤ 2·1.03883·√x, using ψ(x) < 1.03883x (Rosser–Schoenfeld 1962).
5. So the prime part is at most 4.156e^L. Summing the three bounds gives U. ∎

The only arithmetic input is Chebyshev's bound ψ(x) ≪ x. Since Λ ≥ 0, Landau's theorem makes that equivalent to absolute convergence of the Euler product on σ > 1.

**Class 𝒟.** D(s) = Σ a(n)n^{−s}, with real a(n) and a(1) = 1, and these properties:
- D continues meromorphically with finitely many poles;
- Λ_D(s) = γ_D(s)D(s) = εΛ_D(1 − s), where γ_D is a product of Γ_ℝ-factors and Λ_D has order 1;
- −D′/D = Σ c(n)n^{−s} converges absolutely on some line Re s = σ₀.

For such D and real even g ∈ C_c^∞(ℝ), the explicit formula reads

  W_D(g) := Σ_ρ m_ρ ĝ(γ_ρ) = [poles] − Σ c(n)n^{−1/2}(g(log n) + g(−log n)) + arch_D(g),

summed over the non-trivial zeros ρ = ½ + iγ_ρ with multiplicity m_ρ, and absolutely convergent. This is PROVED, standard: Weil (1952). The contour proof shifts from Re s = σ₀ to 1 − σ₀ and uses only absolute convergence of −D′/D on Re s = σ₀, the FE and order 1. For DH, σ₀ = 2 works.

Q_{D,L}(f) := W_D(f * f̃) for f ∈ L²[−L, L] (half-width convention). ζ and DH are in 𝒟.

**Lemma 5.2 (PROVED, qualitative; Landau).** Let D ∈ 𝒟. If there are b ≥ 0, C and L₀ with λ_min(Q_{D,L}) ≥ −Ce^{2bL} for all L ≥ L₀, then every non-trivial zero satisfies Re ρ ≤ ½ + b.

Equivalently, with η_D := limsup_{L→∞} L^{−1} log⁺(−λ_min(Q_{D,L})) and Θ_D := sup Re ρ:

  **η_D ≥ 2Θ_D − 1.**

*Proof.*
0. **The bound extends to all L.** Each Q_{D,L} is bounded below: arch ≥ min Ω_D, the prime part is a finite sum, and the poles are bounded. λ_min(Q_{D,L}) is non-increasing in L (nested windows). So, enlarging C, the bound holds for all L ≥ ℓ/2, for any fixed ℓ > 0.
1. **Test vectors.** Take real φ ∈ C_c^∞(0, ℓ) with φ ≥ 0, g = φ * φ̃, and f_t = φ − τ_tφ for t ≥ 0. Then f_t * f̃_t = 2g − τ_tg − τ_{−t}g. W_D is even (the zeros are symmetric under ρ ↦ 1 − ρ by the FE, and ĝ is even), so Q(f_t) = 2C(0) − 2C(t), where C(t) := W_D(τ_tg).
2. **One-sided bound.** f_t is supported in an interval of length ℓ + t, and W_D is translation invariant (it depends on f only through f * f̃). So Q(f_t) ≥ −Ce^{b(ℓ+t)}‖f_t‖² ≥ −4Ce^{bℓ}‖φ‖²e^{bt}. Hence h(t) := C(0) + Me^{bt} − C(t) ≥ 0, with M = 2Ce^{bℓ}‖φ‖².
3. **Zero expansion.** (τ_tg)^(z) = e^{izt}ĝ(z) and e^{iγ_ρt} = e^{(ρ−½)t}, so C(t) = Σ_ρ m_ρ ĝ(γ_ρ)e^{(ρ−½)t}.
   - Here ĝ(z) = φ̂(z)φ̂(−z) and |φ̂(z)| ≤ C_k e^{ℓ|Im z|}(1 + |z|)^{−k}.
   - The zeros lie in a vertical strip with N(T) ≪ T log T.
   - So Σ m_ρ|ĝ(γ_ρ)| < ∞, the series converges absolutely, and C(t) = O(e^{Bt}), where B = Θ_D − ½.
4. **Laplace transform.** For Re s > max(b, B),

     H(s) := ∫₀^∞ h(t)e^{−st} dt = C(0)/s + M/(s − b) − Σ_ρ m_ρ ĝ(γ_ρ)/(s − (ρ − ½)).

   Termwise integration is justified by absolute convergence. The right side is meromorphic on ℂ, with simple poles at s = ρ − ½ of residue −m_ρ ĝ(γ_ρ), and at 0 and b. The series converges locally uniformly off its poles, since Σ m_ρ|ĝ(γ_ρ)| < ∞.
5. **Choice of φ.** Suppose a zero ρ₀ = ½ + β₀ + iγ₀ has β₀ > b. Choose φ (still ≥ 0) with ĝ(γ_{ρ₀}) = φ̂(γ_{ρ₀})φ̂(−γ_{ρ₀}) ≠ 0. This is possible because, for fixed w, {φ : φ̂(w) = 0} is a closed hyperplane, and φ₀ + εψ with φ₀ > 0 on supp ψ stays ≥ 0. Then H has a pole with real part β₀ > b, so the abscissa of convergence satisfies σ_c ≥ β₀ > max(b, 0).
6. **Landau gives the contradiction.** Since h ≥ 0, Landau's theorem (Widder, *The Laplace Transform*, Ch. II §5) says the real point σ_c is a singularity of H. The only real singularities of H in (b, ∞) are poles s = ρ₁ − ½ coming from **real** zeros ρ₁ > ½ + b.
   - Their residue is −m ĝ(−i(ρ₁ − ½)) = −m (∫φe^{β₁x})(∫φe^{−β₁x}) < 0, because φ ≥ 0.
   - A negative-residue pole at σ_c forces H(σ) → −∞ as σ ↓ σ_c. That contradicts H(σ) = ∫ h e^{−σt} ≥ 0.
   - If σ_c is not such a point, H is analytic at σ_c, which contradicts Landau directly.

   Hence no zero has β₀ > b. ∎

**What is and is not proved.**
- PROVED: the qualitative growth η_D ≥ 2Θ_D − 1, for every D ∈ 𝒟. This is the window-form analogue of the Fixed-Window paper's Landau step (RH map R1) at exponential scale.
- NOT proved: the quantitative form −λ_min(Q_{D,L}) ≥ c·e^{2β₀L} with an explicit c, on the vector e^{−iγ₀x}sinh(β₀x). The gap is precise. Other zeros with real part ≥ β₀ contribute terms of the same or larger exponential order with phases e^{iΔL}, and for DH such zeros accumulate up to σ_DH (Saias–Weingartner density, Bombieri–Ghosh supremum). So a single zero cannot be isolated by this vector without a Turán-type power-sum lower bound, which I have not supplied. Status: CONJECTURE.
- Multiplicity m enters only through the residue; m ≥ 1 changes nothing.

**Corollary 5.3 (PROVED, given Lemma 5.2 and the verified citations).** DH violates U.
- Θ_DH > 1: Davenport–Heilbronn 1936, via Karatsuba 1991 and Saias–Weingartner 2009 Cor. Indeed Θ_DH = σ_DH ≈ 1.120362 (Bombieri–Ghosh 2011 via Odžak–Smajlović 2016).
- So η_DH ≥ 2Θ_DH − 1 > 1, with ≈ 1.2407 at the quoted σ_DH.
- U gives η_ζ ≤ 1.

Hence ζ and DH are separated uniformly in L by the exponent of the window negativity. NUMERICAL consistency check: DH's coefficient mass Σ_{n≤x}|c(n)|n^{−1/2} has local log-slope 0.62–0.63 in x, i.e. exponent 1.24–1.26 in L, matching 2σ_DH − 1 = 1.2407. ζ's slope is ½, i.e. exponent 1 in L.

More generally, among Dirichlet series with periodic coefficients, those with no zeros in σ > 1 are of the form L(s, χ)·(Dirichlet polynomial). This is Saias–Weingartner (2009), Corollary, verified; their theorem also gives c₁T ≤ N(σ₁, σ₂, T) ≤ c₂T for ½ < σ₁ < σ₂ ≤ 1 + η. So among periodic-coefficient fakes in 𝒟, U plus Lemma 5.2 singles out the genuine Dirichlet L-functions.

**Explicit DH zero with Re s > 1: BLOCKER (not computable here).**
- The computed DH zero lists (Spira 1994; Balanzario–Sánchez-Ortiz 2007) contain off-line zeros with Re < 1 only.
- By Bombieri–Ghosh's criterion, a zero in σ > 1 requires log(L(s, χ̄)/L(s, χ)) = 2iθ + (2k+1)πi. That needs the phases p^{−it} of essentially all primes up to ~10³ to align simultaneously, so the first such zeros sit at astronomically large heights (Kronecker, log t ≫ 10²). They cannot be evaluated with mpmath Hurwitz at those heights.
- I therefore did **not** locate one, and claim nothing NUMERICAL or CERTIFIED about a specific zero.

What was run instead (`dh_torus.py`, NUMERICAL):
- **Torus surrogate.** For f = e^{−itx}sinh(0.6x) on [−L, L], the prime part depends on t only through the torus point (p^{−it})_{p<e^{2L}}, and Kronecker realises any such point. Maximising the prime part over that torus by coordinate ascent:
  - DH, P*/‖f‖²: 23.2, 57.3, 111.3, 401.9, 636.3 at L = 3, 3.5, 4, 5, 5.5 (fitted exponent 1.32, against 2σ_DH − 1 = 1.24; the optimiser is heuristic, and L 4.5 underperforms);
  - ζ: 18.6 → 242.7 (exponent 1.03, i.e. e^L, the Chebyshev mass).
  
  **Caveat:** the Kronecker heights needed to align all primes below e^{2L} have log t ≈ π(e^{2L})·log(1/ε) ≫ e^{1.3L}. There the archimedean term Ω_D(t) ≈ log t dominates, so the surrogate does not exhibit negativity. It confirms only the exponent of the prime part. Not a proof of anything.
- **Lemma 5.2's own vectors** (`dh_lemma_vectors.py`; φ a C^∞ bump on (0, 0.5)): R(t) = Q_DH(f_t)/‖f_t‖² stays ≥ 0.889 for t ≤ 12.4 (L ≤ 6.45).
  - This is expected. These low-frequency vectors see the zero at 85.7i only through |ĝ(γ₀)|, which is tiny for a smooth bump. Their negativity is purely asymptotic.
  - The modulated vectors of §5 check 1 (e^{−iγ₀x}sinh(β₀x)) are the efficient ones, negative from L ∈ (1.5, 2].

**Growth-exponent fit for the 0.8085 + 85.70i pair** (`dh_sinh.log`, NUMERICAL):
- −Q/‖f‖² = 3.09, 10.50, 24.44, 52.95 at L = 3, 4, 5, 6, giving local exponents 1.22, 0.85, 0.77 per unit L.
- The prime part's local exponent is 0.67–0.69.
- Both are drifting toward 2β₀ = 0.617, with lower-order (1/L, arch) corrections still visible at L ≤ 6. This pair has 2β₀ < 1, so it does not by itself violate U, consistent with Lemma 5.2.

**Sign-change note** (half-width convention; vector e^{−iγ₀x}sinh(β₀x), γ₀ = 85.699348, β₀ = 0.308517; counted on a 2·10⁶-point grid of [−L, L]):

| L | Re sign changes | Im sign changes | 2Lγ₀/π |
|---|---|---|---|
| 0.5 | 29 | 26 | 27.3 |
| 1.0 | 55 | 54 | 54.6 |
| 1.4 | 77 | 76 | 76.4 |
| 1.42 | 79 | 76 | 77.5 |
| 1.61 | 89 | 86 | 87.8 |
| 1.7105 | 95 | 92 | 93.3 |
| 2.0 | 111 | 108 | 109.1 |

- Exact formulas: Re = 2⌊Lγ₀/π + ½⌋ + 1 (the cos zeros plus the sinh sign change at 0); Im = 2⌊Lγ₀/π⌋ (the sin zero at 0 is a double zero with sinh, so no sign change there).
- Golden Eye's 2Lγ/π ≈ 93 at L = 1.71 matches with **L = half-width**.
- The failing DH vector is thus a modulated "edge" vector carrying ≈ 2Lγ₀/π nodes, sign-changing as S98 requires. It is the opposite of ζ's nodeless soft mode, where Λ ≥ 0 hurts.

**Why U is not RH, and what it does not do (honest).**
- U is equivalent to η ≤ 1. Lemma 5.2 gives η ≥ 2Θ − 1, so U implies "no zeros with Re s > 1" for any D ∈ 𝒟. The converse direction is heuristic. That is classical for ζ.
- It does **not** see DH's strip zeros (0.8085 + 85.70i has 2β₀ = 0.617 < 1), nor F_θ (Re < 1), H₋ₜ or DMV.
- Its value is that it is the first statement on the team's books that is (a) uniform in L, (b) proved for ζ, (c) fatal to DH, and (d) uses infinitely many primes (Chebyshev for all x). It also shows the hybrid-lemma obstruction does not bind once statements are asymptotic in L.

**The ladder (what strengthening costs).**

| Statement | Equivalent to (≥ direction PROVED by Lemma 5.2; ≤ direction heuristic) | Status for ζ | Kills |
|---|---|---|---|
| η ≤ 1 (U) | no zeros in σ > 1 | PROVED | DH (given 5.2) |
| −λ_min = o(e^L) | *not* PNT: zeros at the de la Vallée Poussin boundary 1 − c/log γ with log γ ≈ L² would give ≈ e^L(1 − O(1/L)) | OPEN (**DEAD END** as "PNT-provable", logged) | DH |
| η ≤ δ < 1 | zero-free half-plane Re s > ½ + δ/2 | OPEN (quasi-RH) | DH (δ < 0.617), F_θ (θ > δ/2), DMV (zeros → σ = 1) |
| η = 0 / λ_min ≥ 0 | RH | OPEN | all |

**NUMERICAL sanity checks (cheap; float64/mpmath; not interval).**
1. *DH pair negativity on an explicit vector* (`dh_sinh.log`). Exact Lévy–Khintchine evaluation, no Galerkin. Columns are Q/‖f‖² at L = 1, 1.5, 2, 2.5, 3, 4, 5, 6:

   | Form | L 1 | 1.5 | 2 | 2.5 | 3 | 4 | 5 | 6 |
   |---|---|---|---|---|---|---|---|---|
   | DH | +3.22 | +1.90 | **−0.24** | −1.05 | −3.09 | −10.50 | −24.44 | −52.95 |
   | ζ (same vector) | +3.69 | +3.27 | +3.81 | +3.06 | +3.10 | +1.04 | +0.95 | +1.79 |

   - DH's prime part on this vector grows from 1.0 to 57.2, with local exponent 0.67–0.69 against 2β₀ = 0.617 plus lower-order terms.
   - DH turns negative between L 1.5 and 2, consistent with the Galerkin onset L ≈ 1.71 (S75). The vector is not optimal.
   - **DH onset reconciliation (lead input).**
     - Golden Eye's density heuristic predicts DH failure at L ≈ 1.61. The recorded onset (S75) is 1.7105, with DH numerically PD below it.
     - The pair heuristic of Lemma 5.2 shows why a heuristic cannot settle a 0.1 difference. The onset is where e^{2β₀L}/(2β₀) ≈ C·log(γ₀/2π) with β₀ = 0.3085 and γ₀ = 85.70, and the O(1) constant C (on-line neighbours, arch) moves the predicted onset from L ≈ 0.8 (C = 1) to ≈ 2.6 (C = π). C ≈ 1.8 reproduces 1.71.
     - So 1.61 vs 1.7105 is a constant-level discrepancy, not a contradiction of mechanism. The Galerkin/S75 value should be preferred, provided its PD side at 1.61 ≤ L < 1.71 is N-converged. Given the N caveat in §2, that needs checking: run DH at L = 1.60, 1.65, 1.70 at N = 24/32/40 with LDLᵀ.
     - The explicit vector here (not optimal) only brackets the onset in (1.5, 2.0], which is consistent with both.
2. *Riesz mean at γ₀* (`dh_riesz2.log`). R_L(γ₀) = Σ_{n<e^{2L}} c(n)n^{−1/2}cos(γ₀ log n)(1 − log n/2L) against the Perron residue −e^{2β₀L}/(2Lβ₀²):

   | L | 3 | 4 | 5 | 5.5 | 6 |
   |---|---|---|---|---|---|
   | ratio | 0.53 | 0.71 | 0.82 | 0.86 | 0.89 |

   The ratio tends to 1, and the minimum over t ∈ γ₀ ± 3 sits exactly at γ₀. ζ's Riesz means on t ∈ [10, 200] stay in [−5.8, 1.7] (no growth beyond log-size).
3. *Coefficient mass* (`dh_coeff.log`).
   - Σ_{n≤x}|c_DH(n)|n^{−1/2} = 4953 against ζ's 1262 ≈ 2√x at x = 4·10⁵.
   - The local log-slope is 0.62–0.63, above ½. That is consistent with DH's absolute-convergence abscissa exceeding 1, but not a proof, and the σ > 1 zeros themselves were not located.

## 6. Prior art (verified items marked ✔ with source; the rest are standard background)

- **Weil (1952), Bombieri (2000, "Remarks on Weil's quadratic functional")**: the criterion and the pole terms.
- **Yoshida (1992)**: positivity on small support.
- **Connes (1999, Selecta)**: RH ⇔ the global trace formula with uniform cutoff, i.e. P6.
- **Connes–Consani–Marcolli (2007)**: Weil-proof analogy.
- **Connes–Consani**: archimedean Weil positivity 2021; prolate near-radical 2106.01715; 2310.18423.
- **Burnol**: Sonine spaces; Nyman–Beurling.
- **Li (1997); Bombieri–Lagarias (1999)**: Li coefficients (dead per project list).
- **de Branges (1992); Conrey–Li (2000)**: de Branges' positivity condition fails for ζ.
- **Lagarias (2006)**: de Branges spaces for L-functions.
- **Suzuki**: 1204.1827, 2206.03682 and 2209.04658, screw functions.
- **Mori**: window positivity to log 15, as cited by the project.
- **Kreĭn–Langer**: κ negative squares.
- **Probability**: Beurling–Deny; Fukushima–Oshima–Takeda; Reed–Simon IV §XIII.12; infinite divisibility of the zeta distribution (Khinchine; Lin–Hu 2001); ✔ Nakamura, arXiv:2209.13257, *Zeta distributions generated by Dirichlet series and their (quasi) infinite divisibility* (A(n) > 0 gives infinitely divisible; signed coefficients give at most quasi-infinite divisibility); Bost–Connes (1995).
- **Number theory**:
  - ✔ Davenport–Heilbronn (1936), via Karatsuba (Math. USSR Izv. 36 (1991)) and Bombieri–Ghosh (Russian Math. Surveys 66:2 (2011)). The latter gives the team's DH, with ξ = τ₊ and FE (π/5)^{−s/2}Γ((1+s)/2), identical to `dh_coeff.py`'s coefficients (1, ξ, −ξ, −1, 0).
  - ✔ Saias–Weingartner (2009): ≍ T zeros of DH in every strip ½ < σ₁ < σ < σ₂ ≤ 1 + η. This replaces the dropped Voronin strip citation.
  - ✔ Balanzario–Sánchez-Ortiz, Math. Comp. 76 (2007), doi:10.1090/s0025-5718-07-01999-0, and Spira, Math. Comp. (1994), doi:10.1090/s0025-5718-1994-1254148-8: computed off-line DH zeros (0.808517 + 85.699348i, …). None of the zeros they report has Re > 1.
  - ✔ Conrey–Li, IMRN 2000 no. 18, 929–940: de Branges' positivity fails for ζ (34th zero).
  - ✔ Rosser–Schoenfeld, Illinois J. Math. 6 (1962) 64–94, Thm 12: ψ(x) < 1.03883x for all x > 0 (max of ψ(x)/x at x = 113; sharp constant 1.0388205…, OEIS A206431). U's constant 4.156 ≥ 4·1.03883 = 4.1553 stands.
  - Bombieri–Hejhal (1995);
  - Saias–Weingartner (2009);
  - Rosser–Schoenfeld (1962);
  - Carlson's density theorem;
  - DMV, Math. Ann. 334 (2006); ✔ Broucke–Debruyne–Vindas, Adv. Math. 370 (2020), and Broucke–Debruyne–Révész, Trans. AMS (doi:10.1090/tran/9274); Rodgers–Tao (2020).
- **Novelty, honestly.**
  - The decomposition (1.2) is classical: it is the Lévy–Khintchine reading already in the project, papers_review item 3. Using Perron–Frobenius on the killed arithmetic Lévy generator to explain the nodeless soft mode, and the λ₂/odd-level localisation of window-RH, I have not seen in the literature. Treat it as possibly known.
  - U is classical in content.
  - The twist and polarisation identities for G_θ′ and F_θ are linearity of the explicit formula, and are also the rep-theory seat's.

## 7. Dead ends logged this session

1. **Markov/Lévy positivity (Λ ≥ 0) as the mechanism.** It is neither necessary (G_θ′) nor sufficient (F_θ). DEAD END.
2. **PF ground-state positivity as the forcing property.** Λ ≥ 0 lowers the form on the nodeless soft mode (S98, Lemma 2.1). DEAD END.
3. **Gaussian / OS / KMS field with covariance W.** GNS circularity (P3). DEAD END as a mechanism.
4. **Bost–Connes KMS / zeta-distribution infinite divisibility.** Positive only for β > 1, and DMV shares the structure. DEAD END as stated.
5. **"−λ_min = o(e^L) follows from PNT".** False heuristic: dlVP-boundary zeros at height e^{L²} would saturate e^L. DEAD END.
6. **Counting negative eigenvalues of Q_L over all frequencies as a density detector.** An off-line pair at height γ is visible at window L once e^{2βL}/(2β) ≳ log γ, i.e. over a doubly exponential height range. This only makes sense with frequency localisation. Reformulated as next-step CONJECTURE N below.

## 8. Next

- **(a) Write Lemma 5.2 rigorously.** The same-real-part clustering step via a Turán power sum.
  - Then locate one DH zero with σ > 1 (cheap with mpmath; DH is a combination of two L(s, χ mod 5)).
  - Then extend `dh_sinh.py` to it, to *see* η_DH > 1.
- **(b) CONJECTURE N (frequency-localised density dictionary).** For δ > 0, H → ∞ and L ≥ L₀(δ), let ν(L, T, H; δ) be the number of eigenvalues below −e^{δL} of Q_L compressed to {f ∈ L²[−L, L] : f̂ essentially supported in [T, T + H]}.
  - Conjecture: ν(L, T, H; δ) ≍ #{ρ : Re ρ > ½ + δ/2, γ ∈ [T, T + H]} + O(L).
  - If so, Carlson's density theorem (N(σ, T) = o(T), which uses the Euler product via mollifiers) gives ν_ζ = o(T) unconditionally. Saias–Weingartner give ν_DH ≫ H, and F_θ gives ν ≍ H log T.
  - That would be a PROVABLE, uniform, DH- and F_θ-fatal, non-RH statement that touches DH's *strip* zeros. This is the step after U.
- **(c) λ₂(R_e) and λ₁(R_o) as functions of L.** Fit their decay against the Zhu / Landau–Widom law. Test whether the DH onset (L ≈ 1.71) is where DH's analogue of "λ₂ ≥ 0" breaks. DH has no pole, so it needs λ₁ ≥ 0; that should equally be the *first* level for DH. This is the "one extra bound state" picture, now precise.

RH is OPEN. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
