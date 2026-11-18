# Round 13h — compactifying the scale direction ℝ_{>0}

Daniel's framing: integers are fundamental, infinity is not real, and 0 and ∞ should be glued by x ↦ 1/x. Question: does the Weil form become a norm on some compactified space? Code: `connes_cutoff.py`. Data: `connes_cutoff.json`.

## Test: Connes cutoff Λ = Weil form on test functions supported in [Λ^{−1/2}, Λ^{1/2}], i.e. our width L = log Λ (NUMERICAL)
Smallest eigenvalue λ_min, with the number of negative eigenvalues in brackets, at Galerkin sizes M = 10L / 14L:

| L | ζ | Epstein (own arch) | spliced Beurling (ζ's arch) |
|---|---|---|---|
| 4.0 | −1e-13 (0) | −1.56 / −1.60 (3 / 4) | identical to ζ |
| 4.5 | −1e-13 (0) | −2.88 / −2.94 | **−0.44 / −0.57 (6 / 8)** |
| 5 | −1e-13 (0) | −4.42 | −1.00 / −1.12 (11 / 15) |
| 6 | −2e-13 (0) | −8.61 | −2.30 (21) |
| 7 | −2e-13 (0) | −16.2 | −3.93 / −4.05 (29 / 39) |
| 8 | −2e-13 / −4e-13 (0) | −28.9 / −29.0 (4 / 8) | −5.54 / −5.58 (36 / 46) |

- ζ stays at the float floor with 0 negative eigenvalues up to L = 8. These are the "very small eigenvalues" of Connes–Consani, *Spectral triples and ζ-cycles* (Enseign. Math. 69, 2023).
- The spliced Beurling model is identical to ζ for L ≤ log 55 = 4.007, and **fails immediately after**. The failure is broad (many negative directions), not quartet-localized: once its own primes enter, ζ's arch is the wrong partner, and the model has no arch of its own because it has no functional equation.
- Epstein fails from L ≈ 3.1 (13e) and gets steadily worse.
- The signed bound λ_max(prime − pole) depends on M: for ζ it is 1.55–1.60 at M = 10L and 1.86–1.93 at M = 14L.

## Survey: where each compactification stands and where it fails on the two counter-models
1. **Connes cutoff Λ / semilocal trace formula.** References: Connes, *Selecta Math.* 5 (1999); Connes–Consani, "Weil positivity and trace formula, the archimedean place", *Selecta Math.* 27 (2021); Connes–Consani–Moscovici, "Zeta zeros and prolate wave operators" (2024).
   - Weil positivity for every Λ ⇔ RH. Connes–Consani prove the archimedean-only case: support in [2^{−1/2}, 2^{1/2}], where no primes enter. Their bound is W_∞(g∗g*) ≥ Tr(ϑ(g) S ϑ(g)*), a Sonin-space trace that is positive by construction.
   - **Is it a norm on the compactified space?** For each fixed Λ, yes for ζ numerically (above, and 13b). But fixed Λ is a finite-width statement.
   - **Spliced Beurling:** for any Λ there is a splice that matches ζ exactly up to Λ (PROVED, 13e). It fails only after the splice point, because it lacks the functional equation and its arch.
   - **Epstein:** fails at a finite Λ ≈ e^{3.1} through the signed prime − pole channel (13e).
   - **Content:** uniformity in Λ. Compactifying by cutoff does not remove the limit; it re-indexes it.
2. **Circle of length L = 2 log S (ζ-cycles) and the scaling site [0,∞) ⋊ ℕ^×.** References: Connes–Consani, *Spectral triples and ζ-cycles* (2023); "Hochschild homology, trace map and ζ-cycles" (PSPUM 105, 2023); "Geometry of the scaling site" (2017).
   - This is literally gluing the scale line into a circle. The spectrum of scaling on the ζ-cycle spaces is **contained in the critical zeros** (their Thm 6.4 / 1.1(i)). Stability under the ℕ^× coverings makes the scaling site their parameter space.
   - Not a norm statement: it realizes zeros on the line, and off-line zeros are invisible to it. It neither proves nor refutes RH.
   - **Spliced Beurling:** the ℕ^× covering action needs the monoid of integers acting by n-fold covers. Beurling "integers" give the same abstract free monoid, but the periods log p_B are not logs of integers, and there is no Frobenius on ℤ_max. The construction can be copied only formally.
   - **Epstein:** there is no monoid at all (no Euler product), so no scaling-site sheaf.
3. **Bost–Connes / primon gas** (Bost–Connes, *Selecta* 1995).
   - The ℕ^× semigroup of isometries, H = log n, and Z(β) = ζ(β), with a phase transition at β = 1. It sees only β > 1 and never the critical strip.
   - No Weil norm. The archimedean place, the Γ factor and the functional equation are absent.
   - **Spliced Beurling:** has an identical Bost–Connes-type system (any free abelian monoid with Hamiltonian log n). It **cannot separate** the two.
   - **Epstein:** has no Bost–Connes structure (no multiplicativity). See Round 13i.
4. **Daniel's 0 ~ ∞ gluing by x ↦ 1/x.**
   - In u = log x this is u ↦ −u, and the glued space is the one-point compactification of ℝ, a circle.
   - The Weil form's pole terms carry weights e^{±u/2}, unbounded at both ends. The flip exchanges them, which is exactly the functional equation (Poisson / θ(1/x) = √x θ(x)).
   - **Epstein has this flip** (Poisson on ℤ²), and so does any off-line quartet (13c: functional-equation reversal held for the quartet). So the flip alone **cannot** be the source. **Spliced Beurling lacks it**, which is why it fails right after its splice.
   - Verdict: the flip is necessary (it kills the Beurling splice), not sufficient (Epstein).

## Verdict
None of the four compactifications turns the Weil form into a provable norm.
- **Cutoff:** a norm at each finite Λ (numerically), but the content is uniformity in Λ.
- **ζ-cycles / scaling site:** sees only critical zeros, so it is not a positivity statement.
- **Bost–Connes:** blind to the strip and to the functional equation.
- **x ↦ 1/x gluing:** the functional equation itself, which Epstein has.

This again matches 13e: the multiplicative ℕ^× structure (Bost–Connes, scaling site) and the flip (functional equation) must act together on the same space. The adele class space ℚ^× \ 𝔸_ℚ is where they do, and Weil positivity there is RH. Labels: table NUMERICAL; splice statements PROVED (13e); survey claims LITERATURE; RH OPEN.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
