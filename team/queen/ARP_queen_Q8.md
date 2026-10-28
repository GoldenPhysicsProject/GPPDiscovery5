# §Q.8 The off-line residual and the exact identity R6 would need (Queen of Mathematics; appendix to ARP_queen.md)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Queen of Mathematics (RH sub-team), 2026-10-28. Against GPPDiscovery5 main `8079b30`, which includes:
- G.4d `0244ce4`/`7328e39`;
- §S `9184613` (no-slack lemma S(c)4);
- §P `team/principal-series/ARP_principal-series.md`.

**RH is OPEN. §Q.8 proves nothing about RH.** It does four things:
- (1) It writes the off-line residual as a theorem, with hypotheses stated exactly. The identity on the zero side is PROVED for ℚ (elementary). Its realization as a limit of Connes' cutoff traces is PROVED for function fields (Connes Lemma 3, with one proof step repaired here). For ℚ it is a CONJECTURE, and the missing steps are named.
- (2) It settles the P_Λ vs P̂_Λ placement question as far as δ_Λ is concerned: the placements are equivalent on real f.
- (3) It states the exact identity (R6) that would give Q_L ≥ 0 at every Λ, and lists what that identity is forced to satisfy.
- (4) It checks that F_θ and DH must violate the identity, and at which ingredient.

No numerics were run. The only computer check was a symbolic `sympy` check of four elementary identities (marked ✓sym).

**Corrections found on the way** (these do not change any team label):
- **(c1)** In the proof of Connes 1999 §VIII Lemma 3, the vectors η_z are said to be "asymptotically orthogonal". That is false once two zeros lie on the same side of the line: the limiting Gram entries are Szegő-kernel entries, which are nonzero. The conclusion of the lemma still holds by a basis-free trace argument (Q.8.2).
- **(c2)** For ℚ, Connes' inclusion (23), Q′_{Λ,0} ≤ S_Λ, is *not exact* with the prolate B_Λ. At ∞ no nonzero f has both f and f̂ vanishing beyond Λ (Connes §VIII, the paragraph before (31)). So in §Q.4 item 2, "free positivity of S_Λ − Q′" is exact for function fields only.
  - For ℚ, the exact positive-type object is the compression Δ^♮_Λ (Q.8.3).
  - ARP G.4d's parenthesis "(cokernel deficit, S_Λ − Q′, which is free)" should carry this qualifier.
- **(c3)** §Q.2 calls F_θ an "OS2a GL(1)-constituent" failure; main G.4e/§P.3 now say "polarization at G.4, not a local temperedness filter". These are the same fact. Q.8.5 shows the constituent statement is exactly the adjoint relation V_θ(f)^* = V_{−θ}(f^*), and that relation is what turns the G.4 norm into a polarization. Proposed reading for §Q.2: "OS2a, read as unitarity of the GL(1) representation through which E acts (not a local sign filter)."

Notation is §G/§Q:
- u = log|g| on the ζ-sector L²(C_ℚ)^{Ẑ^*} ≅ L²(ℝ, du).
- f is real in 𝒯_L = C_c^∞((−L, L)), g = f ∗ f̃, f̂(t) = ∫ f e^{itu} du. Then ĝ = |f̂|² on ℝ.
- V is the regular representation, V(f) = ∫ f(a) V(a) da, and f^*(a) = f(−a).
- Q_L(f) = Σ_ρ m_ρ ĝ(γ_ρ), where a zero is written ρ = ½ + iγ_ρ, γ_ρ = γ − iβ, β = Re ρ − ½ ∈ (−½, ½).

---

## Q.8.1 Definitions and the exact zero-side identity (ℚ; PROVED)

**Harmonic measure.** For a zero ρ = ½ + β + iγ:
- μ_ρ := (1/π) |β| / ((t − γ)² + β²) dt on ℝ if β ≠ 0, and μ_ρ := δ_γ if β = 0.
- μ_ρ is the harmonic (Poisson) measure of the point β + iγ with respect to the line iℝ. This is Connes' normalization: his ρ_C = ρ − ½, and his line iℝ is our critical line.
- It is a probability measure. Its Fourier transform is ∫ e^{itu} dμ_ρ(t) = e^{iγu − |β||u|} (✓sym).

**Lemma Q.8.1 (absolute convergence; PROVED).** Let g ∈ C_c²(ℝ) with supp g ⊂ [−2L, 2L]. Then:
- Σ_ρ m_ρ |∫ ĝ dμ_ρ| < ∞ and Σ_ρ m_ρ |ĝ(γ_ρ)| < ∞.
- Explicitly, ∫ ĝ dμ_ρ = ∫ g(u) e^{iγu − |β||u|} du and ĝ(γ_ρ) = ∫ g(u) e^{iγu + βu} du.

*Proof.*
- Integrating by parts twice gives |ĝ(y)| ≤ C_g (1+|y|)^{−2} on ℝ, and |ĝ(γ − iβ)| ≤ C_g e^{2L|β|} (1+|γ|)^{−2}.
- For the smear, split at |y − γ| = |γ|/2:
  - On |y − γ| < |γ|/2, (1+|y|)^{−2} ≤ (1+|γ|/2)^{−2}, and the measure has mass ≤ 1.
  - On |y − γ| ≥ |γ|/2, the Poisson density is ≤ 4|β|/(πγ²) ≤ 2/(πγ²).
  - So |∫ ĝ dμ_ρ| ≤ C′_g (1+|γ|)^{−2}, uniformly in |β| < ½.
- Σ_ρ (1+γ²)^{−1} < ∞ by Riemann–von Mangoldt. The two closed forms are Fubini plus the Fourier transform of μ_ρ. ∎

R-line: R1, per-window | R2, n/a | R3, n/a | R4, n/a | R5, holds for every fake with a zero set of order 1 (all six) | R6, exact.

**Definitions.**
- Smeared zero form: 𝒮(f) := Σ_ρ m_ρ ∫ ĝ dμ_ρ = Σ_ρ m_ρ ∫ |f̂|² dμ_ρ.
- One-zero residual: R_ρ(g) := ∫ ĝ dμ_ρ − ĝ(γ_ρ) = ∫ g(u) e^{iγu} (e^{−|β||u|} − e^{βu}) du. It is 0 when β = 0.
- Total residual: Res(f) := Σ_{ρ: β≠0} m_ρ R_ρ(g).

**Theorem Q.8.A (exact zero-side identity, ℚ; PROVED (restatement; no new RH content), unconditional).** "RH ⇔ Res ≡ 0" below is Weil's criterion re-expressed. For every L > 0 and every real f ∈ 𝒯_L:

  Q_L(f) = 𝒮(f) − Res(f),  with 𝒮(f) ≥ 0.

All three series converge absolutely. Moreover RH ⇔ Res ≡ 0 on 𝒯_L for all L ⇔ Q_L = 𝒮 for all L.

*Proof.*
- The identity is termwise, using Lemma Q.8.1.
- 𝒮 ≥ 0 because ĝ = |f̂|² ≥ 0 on ℝ and μ_ρ ≥ 0.
- RH ⇒ every β = 0 ⇒ Res ≡ 0.
- Conversely, suppose Res ≡ 0. Then Q_L = 𝒮 ≥ 0 for all L, which gives RH by Weil's criterion (Weil 1952; Bombieri 2000). ∎

R-line: R1, the equivalence needs all L | R2, enters only through the explicit formula | R3, the positivity of 𝒮 is the trivial ĝ = |f̂|², not a weight sign | R4, the poles do not appear on the zero side (they live in the explicit-formula expression of Q_L) | R5, 𝒮 ≥ 0 holds for every fake, so it is fake-blind and RH-silent, the zero-side twin of OS2′_MF | R6, exact.

> Consequence for the team: the "residual formula" is not itself the hard part. On the zero side it is an elementary identity. What is non-trivial is to realize 𝒮 as lim_Λ Δ_Λ (Q.8.2/Q.8.3); and even that realization is unconditional, so it is **not** RH-equivalent. RH remains the statement Res ≡ 0, equivalently G.4d.

---

## Q.8.2 Function fields: Connes' Lemma 3 with the Gram step repaired (PROVED)

**Theorem Q.8.B.**
- *Setting.* k is a global field of characteristic p, with C_k ≃ C_{k,1} × q^ℤ. Λ = q^N with N → ∞ through the integers. h ∈ S(C_k) has compact support, so only finitely many characters χ of C_{k,1} occur.
- *Cutoffs.* B_{Λ,0} and Q′_{Λ,0} are as in Connes (21)–(23). Since 1_{O_v} is self-dual, P_Λ and P̂_Λ commute exactly, so E(B_{Λ,0}) ⊂ S_Λ and there is no placement issue.
- *Definitions.* Δ_Λ(h) := Tr((S_Λ − Q′_{Λ,0}) V(h)), and δ_Λ(h) := Δ_Λ(h) − Δ_W(h), with Δ_W Weil's distribution (17).
- *Conclusion.* For each fixed h,

  lim_{N→∞} δ_Λ(h) = Σ_{χ} Σ_{ρ ∈ B/N^⊥, Re ρ ≠ 0, L(χ̃, ½+ρ) = 0} N(χ̃, ½+ρ) [∫_{iℝ} ĥ(χ̃, z) dμ_ρ(z) − ĥ(χ̃, ρ)].

  The sum is finite. Both sides are 0 here, because RH holds for function fields (Weil). The content of the theorem is the structure of the limit.
- *Rate.* O(1/N) = O(1/log Λ), coming from the on-line zeros (Fejér). Off-line terms, if there were any, would converge geometrically.

*Proof.* Connes reduces the claim to a finite model. For each χ, E(B_{Λ,χ}) ⊂ S_{Λ,χ} ≅ ℓ²([−N, N]) is cut out exactly by the conditions (28), one for each zero. The codimension count uses Weil [W1] Thm VII.6. We accept this as Connes states it. For χ = 1 Connes writes "assume first χ ≠ 1"; the χ = 1 sector is as stated by Connes and not re-checked here.

What remains is the finite-model claim: for a finite multiset F ⊂ ℂ^*, E_N := {ξ ∈ ℓ²[−N, N] : Σ_n ξ(n) z^n = 0 for all z ∈ F}, and V the shift,

  Tr((S_N − E_N) V^k) → Σ_{z∈F} c_k(z),

where c_k(z) := ∫_{|u|=1} u^{−k} P(u) dm(u) is the Poisson integral at the point determined by z (as in Connes' (29)).

*Repaired step.*
- S_N − E_N is the orthogonal projection P_N onto W_N = span{w_z}, where w_z(n) := z̄^n restricted to [−N, N]. Then Tr(P_N V^k) = Tr(G_N^{−1} M_N^{(k)}), where (G_N)_{zz′} = ⟨ŵ_z, ŵ_{z′}⟩, (M_N^{(k)})_{zz′} = ⟨V^k ŵ_z, ŵ_{z′}⟩, and ŵ = w/‖w‖. This expression is basis-free; no orthogonality is needed.
- Case |z|, |z′| > 1. The mass of ŵ concentrates at n ≈ +N, and
  - G_N does **not** converge entrywise when zeros have different arguments: the off-diagonal entries carry phases conj(u_z) u_{z′}, u_z := (z/|z|)^N, which oscillate in N. What holds is convergence **up to diagonal unitary conjugation**: G_N = D_N^* G D_N + o(1), with D_N := diag(u_z) unitary and
    G_{zz′} = (1−|z|^{−2})^{1/2} (1−|z′|^{−2})^{1/2} / (1 − (z̄ z′)^{−1}).
  - G is a Szegő-kernel Gram matrix. It is nonzero off the diagonal, so (c1), and positive definite for distinct points.
  - Likewise, up to diagonal unitary conjugation, M_N^{(k)} = D_N^* M^{(k)} D_N + o(1), with M^{(k)} = diag(z̄^{−k}) G for k ≥ 0 and M^{(k)} = G diag(z′^{k}) for k < 0 (only the truncated top strip matters).
  - Tr(G_N^{−1} M_N^{(k)}) is invariant under the simultaneous conjugation by D_N, and G is invertible, so in both cases Tr(G_N^{−1} M_N^{(k)}) = Tr(G^{−1} M^{(k)}) + o(1) → Σ_z c_k(z).
- Case |z|, |z′| < 1 is the mirror image (mass at −N).
- Cross terms between the two sides tend to 0, since their mass sits at opposite ends.
- On-line |z| = 1: ⟨V^k ŵ_z, ŵ_{z′}⟩ = z̄^{−k} δ_{zz′} (2N+1−|k|)/(2N+1) + O(1/N) (Fejér). The cross terms with off-line z tend to 0.
- So G_N = D_N^* G D_N + o(1) up to diagonal unitary conjugation (D_N diagonal unitary over all zeros), with G block-diagonal and invertible; the trace is conjugation-invariant, so it converges.
- Multiple zeros: use w_{z,j}(n) = n^j z̄^n. The limiting Gram matrix is a confluent Szegő kernel, still positive definite, and V^k acts triangularly with diagonal c_k(z). This gives the multiplicity N(·). That case is PROVED in sketch.
- Finally, sum over k with h's Fourier coefficients (finitely many k) and over χ. ∎

R-line: R1, a limit Λ → ∞ for fixed h, not uniform in h | R2, (28) encodes the Euler product, and (22) encodes the FE | R3, n/a | R4, the χ = 1 pole sector is as stated by Connes | R5, n/a (function fields) | R6, an exact identity at each N, with an exact limit.

---

## Q.8.3 ℚ: placement of the cutoff, the exact positive object, and the conjectural realization

**Placements.**
- (P): B^{P}_{Λ,0} is the span of the even prolates ψ_n ⊗ 1_Ẑ with n ≤ 4Λ² (Connes (31)–(33)), so ran P_Λ means support exactly in |x| ≤ Λ, intersected with the two boundary conditions φ(0) = φ̂(0) = 0.
- (P̂): B^{P̂}_{Λ,0} := F_𝔸(B^{P}_{Λ,0}), i.e. exactly band-limited; again both conditions hold, since F swaps them.

Let Q′^{•}_Λ := E Q^{•}_Λ E^{−1} for • ∈ {P, P̂}.

**Lemma Q.8.2 (where the boundary data sit; PROVED, Poisson).** For f ∈ 𝒮(𝔸):

  E(f)(g) = E(f̂)(g^{−1}) + f̂(0) |g|^{−1/2} − f(0) |g|^{1/2}.

Consequences:
- Both pole functionals live at the small-|g| end (u → −∞): f̂(0) e^{−u/2}, which is not in L², and −f(0) e^{u/2}.
- In (P), E(f) vanishes exactly for |g| > Λ, because |qg| = |g| for q ∈ ℚ^* and f vanishes beyond Λ. The leak out of S_Λ is at |g| < Λ^{−1}, the same end as the pole data.
- In (P̂) the leak is at |g| > Λ.
- A one-sided placement (only φ(0) = 0 or only φ̂(0) = 0) is excluded. If f̂(0) ≠ 0 then E(f) ∉ L²(C). If f̂(0) = 0 but f(0) ≠ 0, then (22) E∘F = inv∘E fails, so Θ-symmetry is lost.
- So both conditions are forced.

*Proof.* Poisson summation over ℚ ⊂ 𝔸, Σ_{q∈ℚ} f(qg) = |g|^{−1} Σ_{q∈ℚ} f̂(q g^{−1}), with the q = 0 terms separated out. ∎

**Proposition Q.8.3 (placement independence; PROVED).**
- Define Δ^{•}_Λ(h) := Tr((S_Λ − Q′^{•}_Λ) V(h)). This is well defined: Q′ has finite rank, and S_Λ V(h) is trace class for smooth compactly supported h.
- Then Δ^{P̂}_Λ(h) = Δ^{P}_Λ(h∘inv).
- For h = f ∗ f^* with f real (every case the Weil form needs), h∘inv = h. Hence **Δ^{P̂}_Λ = Δ^{P}_Λ and δ^{P̂}_Λ = δ^{P}_Λ on real f.**
- For complex f the two placements give complex-conjugate values.
- The same holds for the compressed Δ^♮ below.

*Proof.*
- By (22), E∘F_𝔸 = Θ_C∘E on 𝒮(𝔸)_0, so Q′^{P̂} = Θ_C Q′^{P} Θ_C.
- Θ_C S_Λ Θ_C = S_Λ, since the window is symmetric.
- Θ_C V(h) Θ_C = V(h∘inv).
- By trace cyclicity, Tr((S − Θ Q′ Θ) V(h)) = Tr((S − Q′) Θ V(h) Θ).
- h = h^* gives h(g^{−1}) = conj h(g), which equals h(g) for real f. ∎

So the open question "pole at P_Λ or P̂_Λ" **does not affect δ_Λ** on the forms the GAP uses. It does affect where the leak (1 − S_Λ)Q′ sits, which matters for step M1 below.

R-line: R1, every Λ | R2, (22) is the FE half | R3, n/a | R4, both pole conditions are forced (Lemma Q.8.2) | R5, n/a | R6, exact.

**Lemma Q.8.4 (the exact positive object for ℚ; PROVED).** Let A_Λ := S_Λ − S_Λ Q′_Λ S_Λ, for either placement and any closed subspace in place of E(B). Then 0 ≤ A_Λ ≤ S_Λ, and

  Δ^♮_Λ(f ∗ f^*) := Tr(A_Λ V(f) V(f)^*) = ‖A_Λ^{1/2} S_Λ V(f)‖²_{HS} ≥ 0 for every Λ.

In function fields A_Λ = S_Λ − Q′_{Λ,0}, so Δ^♮_Λ = Δ_Λ. For ℚ, Δ_Λ − Δ^♮_Λ = −Tr((Q′ − S Q′ S) V(h)), which is controlled by the leak (1 − S_Λ) Q′_Λ only.

*Proof.* S Q′ S ≤ S. Also S_Λ V(f) has a smooth compactly supported kernel, hence is Hilbert–Schmidt. ∎

R-line: R1, every Λ | R2, it uses E | R3, n/a | R4, n/a | R5, positive type for any datum with an E (so fake-blind) | R6, exact.

**Conjecture Q.8.C (Connes' Lemma 3 for ℚ, smeared form; CONJECTURE, unconditional and not RH-equivalent).** For either placement (equivalent by Q.8.3), every fixed real f ∈ C_c^∞(ℝ_u) (ζ-sector, h = f ∗ f^*), and Λ → ∞ through ℝ:

  Δ^♮_Λ(h) → 𝒮(f) = Σ_ρ m_ρ ∫ |f̂|² dμ_ρ, and hence δ^♮_Λ(f) := Δ^♮_Λ(h) − Q(f) → Res(f).

The mode is pointwise in f; no uniformity in f is claimed. Equivalently:

  Tr(S_Λ Q′_Λ S_Λ V(h)) = 2‖f‖² log Λ − 𝒮(f) + o(1).

Status: Connes states this as "the analogue of Lemma 1, Theorem 5 and Lemma 3" (§VIII, after (33)); no proof is cited. Missing steps:
- **M1 (leak).** Show Tr((Q′_Λ − S_Λ Q′_Λ S_Λ) V(h)) → 0. This needs trace-norm control of (1 − S_Λ)Q′_Λ against V(h): the Landau–Pollak–Slepian transition band has about log Λ states with eigenvalue strictly between 0 and 1 (Connes, §VIII after (32)). In (P) the leak is at the pole end (Lemma Q.8.2).
- **M2 (infinitely many zero vectors).** Over ℚ, S_Λ ⊖ S_Λ E(B) is infinite-dimensional. Show that its trace against V(h) splits as Σ over zeros with |γ| ≤ T(Λ) of the Gram/Poisson limit (the Q.8.2 argument, now with |F| → ∞), plus a remainder equal to 2‖f‖² log Λ minus the archimedean principal value, up to o(1). Zeros with log|γ| ≳ L are not resolved by the window, so they must be absorbed into the counting term (Connes (34)). The Q.8.2 argument does not cover this: in function fields F is finite.
- **M3 (uniform Gram).** Resolved zeros need the limiting Gram matrices (Szegő kernels on accumulating points) to be uniformly invertible as their number grows.

Even if proved, Conjecture Q.8.C would leave RH exactly as open as before: RH ⇔ Res ≡ 0 (Thm Q.8.A).

R-line: R1, the Λ-limit is exactly Connes' "uniformity" problem (§0.2) | R2, E carries both halves | R3, n/a | R4, the pole conditions are forced and Θ-symmetric | R5, n/a (it is a statement about ζ's own cutoff traces) | R6, an exact identity plus o(1).

**Theorem Q.8.D (structure on an off-line quartet; PROVED, exact, ✓sym).** Let ρ₀ = ½ + β + iγ₀ with 0 < β < ½, ρ₀′ = 1 − ρ̄₀ = ½ − β + iγ₀, and multiplicity m. Put χ_L(x) = e^{−iγ₀x} cosh βx and φ_L(x) = e^{−iγ₀x} sinh βx, with ⟨f, ψ⟩ = ∫ f ψ̄. Then for real f ∈ 𝒯_L:
1. **Weil (point) part.** ĝ(γ_{ρ₀}) + ĝ(γ_{ρ₀′}) = 2∫ g(u) e^{iγ₀u} cosh βu du = 2|⟨f, χ_L⟩|² − 2|⟨f, φ_L⟩|². This is the spec's quartet form: complex signature (1,1).
2. **Smeared part.** ∫ ĝ dμ_{ρ₀} + ∫ ĝ dμ_{ρ₀′} = 2∫ g(u) e^{iγ₀u} e^{−β|u|} du = 2∫ |f̂|² dμ_{ρ₀} ≥ 0. Both members of the pair have the same harmonic measure.
3. **Residual.** R_{ρ₀} + R_{ρ₀′} = −2∫ g(u) e^{iγ₀u} sinh(β|u|) du, using cosh(βu) − e^{−β|u|} = sinh(β|u|) (✓sym). Equivalently:

   R_{ρ₀} + R_{ρ₀′} = 2∫ |f̂|² dμ_{ρ₀} − 2|⟨f, χ_L⟩|² + 2|⟨f, φ_L⟩|².

   So **the residual is the quartet (1,1) form with its sign reversed, plus a positive Poisson smear.** Adding the conjugate pair (γ₀ → −γ₀) gives the real quartet residual 2Re(·) = −4m ∫ g(u) cos(γ₀u) sinh(β|u|) du.
4. **Sign structure.**
   - On the quartet's negative direction (⟨f, χ_L⟩ = 0), the residual is ≥ 2|⟨f, φ_L⟩|² > 0. This is exactly the amount by which the Weil quartet is negative there: the off-line zero enters Δ_∞ as a positive resonance and Q as an indefinite point pair.
   - As a real form on real f, the quartet residual has negative index ≤ 2, matching the real (2,2) signature of the Weil quartet (`QuartetSignature.Q_neg_plane`). It is genuinely indefinite for every L > 0: take f ≥ 0, f ≢ 0, supported in [−ε, ε] with 4ε|γ₀| < π. Then g ≥ 0 on |u| ≤ 2ε, cos(γ₀u) > 0 and sinh(β|u|) > 0 there, so the quartet residual is < 0.
   - Growth along φ_L: |⟨f, φ_L⟩|² for normalized f reaches order e^{2βL}/(2β), the Turán Theorem A rate.
5. **Comparison with the pole (orientation only, not an identification).** The pole pair splits the other way: 2cosh(u/2) = e^{|u|/2} + e^{−|u|/2} (✓sym), i.e. the PNT mean plus a unit-mass smear of height ½ (G.3a: MF_L(f) = ∫ g e^{−|u|/2}). The pole is not in Connes' open strip, so it is not a Lemma 3 term. The two decompositions differ in mass (one smear for two poles, versus one smear per zero), and no structural claim is made.

*Proof.*
- Items 1–3: write g(u) = ∫ f(u+v) f(v) dv, so ∫ g K = ∫∫ f(w) f(v) K(w − v). Then cosh(β(w−v)) = cosh βw cosh βv − sinh βw sinh βv.
- The Fourier transform of μ_ρ is e^{iγu − |β||u|} (✓sym), and ĝ = |f̂|² on ℝ.
- Item 4: Cauchy–Schwarz and min–max, using the rank-two negative term −2|⟨f, χ_L⟩|² (real and imaginary parts). The bump argument is pointwise. ∎

R-line: R1, the indefiniteness holds at every L; the growth is the L → ∞ statement | R2, n/a | R3, the residual's sign does not depend on the prime weights | R4, the comparison with the pole is orientation only | R5, it locates the residuals of DH, Epstein, F_θ (spec quartet numerics, NUMERICAL, not re-run) | R6, exact.

---

## Q.8.4 The exact identity R6 would need (GAP/CONJECTURE)

**Lemma Q.8.5 (operator form of the mean field; PROVED, ✓sym).** Let η(u) := e^{−u/2} 1_{u≥0}. Then ‖η‖ = 1 (the soft mode, i.e. the half-line pole state of OS2′_MF), and

  MF_L(f) = ⟨V(f) V(f)^* η, η⟩ = ‖V(f^*) η‖² = Tr(|η⟩⟨η| V(f) V(f)^*).

*Proof.* ∫ η(x) η(x+r) dx = e^{−|r|/2} (✓sym), and ‖V(f^*)η‖² = ∫∫ f(b) f(b′) (η ⋆ η)(b − b′). ∎

So G.3a in operator form reads Q_L(f) = Tr(|η⟩⟨η| V V^*) + Fl_L(f), and Fl_L is not PSD (§S S(b)4).

**Conjecture Q.8.E (R6 identity; this is GAP G.4d in operator form).** There is one explicit rule L ↦ D_L, built from E, S_Λ, P_Λ, P̂_Λ and Q′_Λ with Λ = e^L, producing bounded operators D_L on L²(ℝ_u) such that, for every L > 0:
- **(I1) exact identity:** Q_L(f) = Tr(D_L V(f) V(f)^*) for all real f ∈ 𝒯_L, with D_L S_Λ V(f) Hilbert–Schmidt;
- **(I2) positivity:** D_L ≥ 0;
- **(I3) strong form:** D_L ≤ A_Λ = S_Λ − S_Λ Q′_Λ S_Λ.

Labels:
- (I1)+(I2) for all L ⇒ Q_L ≥ 0 for all L ⇒ RH. PROVED implication; the existence is the GAP.
- (I1)+(I3) gives δ^♮_Λ(f) = Tr((A_Λ − D_L) V V^*) ≥ 0 **at every Λ = e^L**, which is the "δ_Λ ≥ 0 at every Λ" form. Together with (I2) it gives 0 ≤ δ^♮_Λ(f) ≤ Δ^♮_Λ(f ∗ f^*). It is CONJECTURE (strong). (I3)+(I2) ⇒ RH (via (I1)). Whether RH ⇒ (I3) is **unknown**: under RH, δ^♮_Λ → 0 by Conjecture Q.8.C, but its sign at finite Λ is not determined (same direction issue as H1 in §C). No claim that (I3) is strictly stronger than RH is made.
- **Uniformity (PROVED, trivial).** Q_{L′} restricted to 𝒯_L equals Q_L, so (I1) at L′ implies (I1) at L with the compression of D_{L′}. "Uniform in Λ" therefore means exactly "one rule for all L" (R1).
- Under RH, some positive D with (I1) exists abstractly (GNS: PROVED-LEAN `gns_from_positive_definite`), but it is not explicit. R6 asks for the explicit rule.

**Proposition Q.8.F (what (I1) is forced to satisfy; PROVED, conditional on (I1)+(I2)).**
- **(N1) No slack on the soft mode.** By §S S(c)4, r(L) = inf_f Q_L/MF_L is nonincreasing with lim r(L) ≤ 0. So for every c > 0 there are L and f ∈ 𝒯_L with Tr((D_L − c|η⟩⟨η|) V V^*) < 0. No positive multiple of |η⟩⟨η| splits off D_L with a positive remainder, uniformly in L. The identity must hold with **zero margin, tightest on η** (dilated nodeless bumps; §S).
- **(N2) No "pole square + positive remainder."** D_L − |η⟩⟨η| represents Fl_L, which is not PSD (S(b)4). So the pole square must be spent exactly, and the negative part of the Fl-operator must be absorbed by it in the limit. This is the operator form of OS4 being fully consumed (S103).
- **(N3) Zero side.** Tr(D_L V V^*) = 𝒮(f) − Res(f) (Thm Q.8.A), so (I1)+(I2) for all L forces Res ≡ 0. Hence D_L cannot be "𝒮 plus a positive correction" unless that correction vanishes in the limit.

*Proof.* Direct from S(c)4, S(b)4, Lemma Q.8.5 and Thm Q.8.A. ∎

R-line for Q.8.E/F: R1, all L, one rule | R2, E is the only ingredient carrying both halves (P.5–P.6: the ℚ^*-sum, Meyer Lemma 5.3, is the only non-isometric step) | R3, the positivity is (I2), an operator statement, not a weight sign | R4, η is in the positive space and has zero margin (N1) | R5, see Q.8.5 | R6, it *is* the exact identity.

---

## Q.8.5 Fakes against the identity: F_θ and DH must fail, and where

**F_θ = ξ(s−θ)ξ(s+θ), 0 < θ < ½ (PROVED).**
- E commutes with |x|^{±θ}. The C-representation through which E acts on the constituents is V_{±θ}(a) := |a|^{±θ} V(a), and V_{±θ}(f) = V(e^{±θu} f).
- **Adjoint relation (this is the constituent statement).** V_θ(f)^* = ∫ f(a) |a|^θ V(a^{−1}) da = V_{−θ}(f^*). So the reflection f ↦ f^* *swaps* the constituents rather than preserving each one unitarily. That is the precise content of "OS2a fails on GL(1) constituents" (§Q.2), and it is what main G.4e/§P.4 call norm → polarization at the gluing (c3).
- If ζ satisfied (I1) with D_L, polarization gives B_ζ(f₁, f₂) = Tr(D_L V(f₁) V(f₂)^*), and §Q.2 (PROVED) gives

  Q_{F,L}(f) = 2Re B_ζ(e^{θu}f, e^{−θu}f) = 2Re Tr(D_L V_θ(f) V_{−θ}(f)^*).

  This is a polarization. It would be a trace of the (I1) form only if V_{−θ}(f) = V_θ(f), i.e. θ = 0.
- **No positive D′ represents Q_{F,L} in any form Tr(D′ W(f) W(f)^*):** Q_{F,L} is indefinite for large L (zeros at ½ ± θ + iγ; Lemma 5.2).
- The GL(2) product is unitarizable: the Knapp–Stein form is positive definite on the complementary series (§P.3). That does not help, because (I1) lives in L²(C_ℚ), the GL(1) object, where the only involution available is the adjoint above.
- **Failing ingredient:** OS2a in its GL(1)-adjoint form, consistent with §Q.2 and §P.3–P.4.
- R-line: R1, the failure holds for all large L | R2, F_θ has both halves | R3, F_θ's weights are positive and it still fails | R4, MF ≥ 0 holds | R5, ✓ | R6, exact.

**DH (PROVED, with the uniqueness step standard).**
- Any D_L built from E restricts, in the conductor-5 sectors χ and χ̄, to cokernel traces whose spectral sides (in the limit) are supported on the zeros of L(χ̃, ·): (28) in Connes' setting; Meyer 2005 spectral realization.
- **Unique decomposition.** A Weil-type functional ĝ ↦ Σ m_ρ ĝ(γ_ρ) on all ĝ = (f₁ ∗ f̃₂)^, f_i ∈ C_c^∞, determines the zero multiset. The reason is that exponentials e^{iγ_ρ u} with distinct complex frequencies are linearly independent as distributions (Laplace-transform uniqueness). This is standard and is used as PROVED.
- So a combination of the sector forms Q_χ, Q_χ̄ has spectral support inside Z(L(χ)) ∪ Z(L(χ̄)).
- DH has zeros in σ > 1 (Davenport–Heilbronn 1936), where L(χ, s) L(χ̄, s) ≠ 0 because the Euler product converges absolutely. So **Q_DH is not any combination of E-sector traces, and (I1) fails for DH.**
- The product L(χ) L(χ̄) (a degree-2 Euler product, form Q_χ + Q_χ̄) *passes* this step. So the failing ingredient is **OS3**: a single ⊗′ Euler datum, as opposed to a sum of functions. This is consistent with §Q.2 and with §P.1/P.5 ("no single ⊗′").
- In addition DH has no pole, so η is absent and (N1) is vacuous (R4 fails, §Q.2). Its form would have to be Fl-only, and that form is indefinite.
- R-line: R1, the failure already holds for fixed L once DH's σ > 1 zero is resolved (Lemma 5.2) | R2, DH has the FE but not the Euler product | R3, n/a | R4, fails (no pole) | R5, ✓ | R6, exact.

**G_θ′ (control).** The unitary twists |·|^{±iθ′} satisfy V_{iθ′}(f)^* = V_{−iθ′}(f^*), and V_{±iθ′} are unitary. Q_G = 2Q_ζ(f cos θ′u) + 2Q_ζ(f sin θ′u) is a *sum of two (I1)-forms*. So (I1)+(I2) for ζ ⇒ (I1)+(I2) for G_θ′. **It passes**, as the control must.

---

## Q.8.6 Summary of labels

| item | label |
|---|---|
| Lemma Q.8.1 absolute convergence; Thm Q.8.A Q = 𝒮 − Res, 𝒮 ≥ 0, RH ⇔ Res ≡ 0 | PROVED (restatement; no new RH content) — Weil's criterion re-expressed |
| Thm Q.8.B function-field residual limit (Connes Lemma 3, Gram step repaired) | PROVED (simple zeros); multiple zeros PROVED in sketch; χ = 1 as stated by Connes |
| Lemma Q.8.2 boundary data at the small-\|g\| end; one-sided placements excluded | PROVED |
| Prop Q.8.3 P_Λ ↔ P̂_Λ placement equivalence on real f | PROVED |
| Lemma Q.8.4 Δ^♮_Λ positive type for ℚ at every Λ | PROVED |
| Conj Q.8.C ℚ realization lim Δ^♮_Λ = 𝒮 (missing M1–M3) | CONJECTURE (unconditional; not RH-equivalent) |
| Thm Q.8.D quartet structure of the residual | PROVED (exact) |
| Lemma Q.8.5 MF_L = ‖V(f^*)η‖² | PROVED |
| Conj Q.8.E (I1)+(I2) explicit identity | GAP (⇒ RH; equivalent to G.4d with an explicit D_L) |
| Conj Q.8.E strong form (I3), δ^♮_Λ ≥ 0 at every Λ | CONJECTURE ((I3)+(I2) ⇒ RH; RH ⇒ (I3) unknown) |
| Prop Q.8.F (N1)–(N3) | PROVED, conditional on (I1)+(I2) |
| F_θ fails (I1) at OS2a-adjoint; DH fails at OS3; G_θ′ passes | PROVED |

No new Lean names are cited. The names used are those of §Q (ARP_queen.md "Q.L"), plus `QuartetSignature.Q_neg_plane` and `PositiveDefiniteGNS.gns_from_positive_definite`.

References:
- Connes 1999, Selecta Math. 5, 29–106 (math/9811068): §VIII (16)–(33), Lemma 3 and its proof (28)–(29).
- Meyer 2005, Duke Math. J. 127, 519–595: §5.4, Lemma 5.3 (via §P).
- Weil 1952; Bombieri 2000; Weil, *Basic Number Theory* Thm VII.6.
- Landau–Pollak–Slepian 1961–62.
- Davenport–Heilbronn 1936.
- Knapp–Stein 1971 (via §P.3).
- Team: ARP §G/G.4d; §S S(b)4, S(c)4; §P.3–P.6; spec §5; `2026-10-22_turan_lemma52.md`.
