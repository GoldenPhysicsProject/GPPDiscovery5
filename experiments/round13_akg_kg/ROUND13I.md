# Round 13i — celestial bulk from the integers: primon gas / Bost–Connes Fock space with Pasterski–Shao

Code: `primon.py`. Data: `primon.json`. Bulk: ℋ = ℓ²(ℕ) = Sym(ℓ²(primes)), basis |n⟩, H|n⟩ = log n |n⟩, Z(β) = Tr e^{−βH} = ζ(β) (Julia 1990; Spector, CMP 127 (1990); Bost–Connes 1995).

## 1. Pasterski–Shao transform
- PS take the Mellin transform of plane waves in energy, φ_Δ = ∫ dω ω^{Δ−1} e^{iωq·X}. Here the energy ladder is the scale n itself, so the Mellin transform of the eigenbasis is
  |Δ⟩ := Σ_n n^{−Δ}|n⟩ = e^{−ΔH} Σ_n |n⟩.
  This is time evolution e^{−itH} continued to complex time, applied to the uniform vector.
- **Inner product inherited from Fock unitarity:** ⟨Δ|Δ′⟩ = Σ n^{−(\bar Δ + Δ′)} = ζ(\bar Δ + Δ′). (PROVED, one line.)
  - Normalizable iff Re Δ > 1/2, with ‖|Δ⟩‖² = ζ(2 Re Δ).
  - On **Re Δ = 1/2** it is delta-normalizable but **not diagonal** (PROVED):
    ⟨1/2+it | 1/2+it′⟩ = lim_{ε↓0} ζ(1+ε+i(t′−t)) = π·Res·δ(t′−t) − i·PV 1/(t′−t) + [ζ(1+iu) − 1/(iu)]_{u=t′−t}.
    The last bracket is smooth and non-zero (= γ at u = 0).
- So the line Re Δ = 1/2 does come out of unitarity. But it is **kinematic**: it is half the pole at 1, the Mellin–Plancherel line of any system whose state count is N(x) ~ c·x. It is not tied to the zeros.
- **The zeros never appear in this inner product.** The kernel samples ζ on Re = 1, never on the critical line.

## 2. Is the Weil form this inner product? No — wrong side of the explicit formula (PROVED, structural)
- A Fock norm Σ_n |⟨n|Ψ_g⟩|² of a vector linear in g is a translation-invariant (Toeplitz) form only if ⟨n|Ψ_g⟩ = c_n ĝ(log n). In that case the form is Σ_n |c_n|² |ĝ(log n)|²: positive point masses at the frequencies log n.
- The Weil form is ∫∫ g(x)g(y)k(x−y) with k̂ = Σ_ρ δ_{γ_ρ}. Its positive point masses sit at the zero ordinates. The log n (in fact log p^k with weights Λ(n)n^{−1/2}) appear on the other side, as δ's in x, with a minus sign, alongside the arch and pole terms.
- Fock unitarity makes the arithmetic-side frequencies positive. Weil positivity needs the zero-side frequencies positive. Making the primon gas carry the Weil form would require a Hilbert space whose dilation spectrum is {γ_ρ}: the Hilbert–Pólya problem, not the primon gas.
- **The zeros enter as resonances.**
  - They are Fisher / Lee–Yang zeros of Z(β) = ζ(β) continued past the Hagedorn point β = 1.
  - Equivalently, they are poles of the supersymmetric primon gas's Witten index Tr(−1)^F n^{−β} = 1/ζ(β) (Spector).
  - The continuation into Re β < 1 requires the functional equation and arch, which the Fock space does not contain (no Γ factor, no archimedean place).

## 3. Counter-models: what distinguishes them?
- **Spliced Beurling:** a Fock space Sym(ℓ²(P_B)) with H = dΓ(log p_B). It is unitarily the same kind of object.
  - Its Gram kernel ζ_B(\bar Δ + Δ′) is delta-normalizable on Re Δ = 1/2 with weight π·Res_B. Res_B = 1.003 / 0.998 / 0.998 from (s−1)ζ_B(s) at s = 1 + 10⁻² / 10⁻³ / 10⁻⁴ (NUMERICAL).
  - Nothing in the Fock structure tells log p_B from log p. "Is p_B an integer?" is an additive question, and the Fock space only knows the multiplicative monoid.
- **Epstein:** ℓ²((ℤ²∖0)/±) with H|v⟩ = log Q(v)|v⟩ and Z = Σ a(n) n^{−β}. It is delta-normalizable on Re Δ = 1/2 with weight π·π/(2√5) = π·0.7025.
  - It is **not** a Fock space: no tensor factorization, because there is no Euler product.
  - It does have ℤ² Poisson summation, hence the functional equation.
- **Witten-index field** F_μ(u) = Σ μ_F(n) n^{−1/2} g(log n − u), L = 2, u ≤ 13, Fock / coefficient data only (NUMERICAL):
  - ζ: rms 0.0019, flat (growth rate −0.003/u).
  - Epstein: rms ≈ 1.3, flat at this width. It is about 700× ζ's because the Dirichlet inverse of r_Q/2 has unbounded, non-multiplicative coefficients with no Möbius cancellation. The off-line growth e^{0.433u} is not yet resolved at u ≤ 13, the same as in 13e at L = 2.
  - Spliced Beurling: about ζ's size up to u ≈ 4. Then a large transient decays like e^{−u/2} (rms 0.078 at u = 5–7 down to 0.0033 at u = 11–13). This comes from the continuous comparison model, whose 1/ζ_C is large left of Re s = 1. The quartet term, about 1e-5 at this width, is not resolved.
  - **No separation from Fock data at accessible u.** As with the 13e field, separation is a large-u (tail) statement.

## 4. Verdict
- **Right:** the PS construction on the primon Fock space does single out Re Δ = 1/2 as the delta-normalizable line, exactly from unitarity.
- **Wrong as a source:**
  - That line is kinematic: it comes from the pole at 1 and is shared by Beurling and Epstein.
  - The inner product is ζ on the 1-line, not on the critical line, and it is not diagonal (there is a PV off-diagonal part).
  - The Weil form is not this inner product: it sits on the dual side of the explicit formula.
  - The zeros appear only as resonances (Lee–Yang zeros of Z, poles of the Witten index 1/ζ), reached by analytic continuation that needs the functional equation.
- **Precisely where it fails:** the Fock space carries the Euler product (the multiplicative ℕ^×) and nothing else. The functional equation / arch, i.e. the additive ℤ structure, is absent. So Beurling shares everything, and nothing distinguishes ζ inside the Fock space. That matches the 13e sentence: the missing ingredient is an additive operator acting on the same ℓ²(ℕ), e.g. the shift |n⟩ ↦ |n+1⟩ or Poisson / θ(1/x), coupled to the Fock grading.
- Labels: §1, §2 PROVED (elementary); residues and fields NUMERICAL; resonance picture CHECKED (Spector / Julia); RH OPEN.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
