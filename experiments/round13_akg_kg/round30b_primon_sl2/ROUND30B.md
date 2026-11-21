# Round 30-b — the primon gas and SL(2,ℝ): is the gas unitary in Δ = s?

Approved by Daniel Toupin in his own words: "Keep going don't stop find the way" (t551u, 2026-11-20).
Seat: Heart of the Cosmos / Grok (D5). No zeta zeros are used as inputs anywhere. Mpmath zeros are not called in this round.
The Mazáč–Paulos papers (AFB I, arXiv:1803.10233; AFB II, arXiv:1811.10646) are used for their window mechanism only, not for positivity.
Weil positivity is never assumed. It appears only as the statement whose truth would be the output.

**Verdict.** The primon gas gives a genuinely positive Hilbert space: each positive integer is one basis vector, and nothing is negative. But that positivity lives at temperatures β = Re Δ > 1, where the zeros are invisible.
The gas Hamiltonian log N cannot carry any SL(2,ℝ) action (PROVED no-go, §2). So the conformal weight Δ = s has to be the complex inverse temperature, and the zeros are Lee–Yang/Fisher zeros of the gas.
The SL(2) pieces that do exist are these:
- the dilation flow, intrinsic to the gas;
- the Weyl/shadow element, which is the Planck↔Matsubara (KMS) duality, i.e. the functional equation (PROVED, §3);
- a per-prime metaplectic, weight-½ "spinor" oscillator whose two sectors Δ = ¼, ¾ are shadows (PROVED, §2).

Unitarity of the resulting principal-series Casimir is exactly κ = 0, where κ is the Pontryagin (ghost) index of the Weil form. That is RH (Weil, Bombieri; §4). So this route reduces to Hilbert–Pólya/Weil at a precisely named point. **RH is not closed.**

Two new, honest outputs came out of it.
- (i) **A four-principle separation table (§5).** Each of the four standard fakes drops exactly one of Daniel's principles. F_θ keeps positivity, the Euler product, the FE and a unitary shadow intertwiner, and still fails. So no argument built from those four soft properties alone can prove RH.
- (ii) **The MP window mechanism on the principal series (§6).** It is the Hardy–Littlewood sign-change method. It kills F_θ outright and unconditionally, but DH and Epstein pass it, so it is strictly weaker than RH.

The sharpest next GAP is the one-sided, diagonal form of unitarity that F13S (GPPVerify #653, merged `c81c031`) already reduces to: the integer floor D(4 log k) ≥ −B (§7).

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

---

## 1. The candidate Hilbert space and where positivity sits

**H_F = Γ_s(ℓ²(primes)) ≅ ℓ²(ℕ_{≥1})** is the symmetric Fock space of one boson species per prime p.
- The occupation numbers (a_p) correspond to n = ∏ p^{a_p} by unique factorization, so the orthonormal basis is |n⟩, with **each positive integer counted once**.
- The inner product is ⟨m|n⟩ = δ_{mn}. It is positive by construction, with no negatives (PROVED, trivial).
- The Hamiltonian is H|n⟩ = log n |n⟩, and Z(β) = Tr e^{−βH} = Σ n^{−β} = ζ(β) for β > 1. This is the Julia/Spector/Bost–Connes primon gas.

**Coherent (Mellin) vectors.** Set v_s = Σ n^{−s̄}|n⟩ ∈ H_F for Re s > ½.
- Their Gram is ⟨v_s|v_{s'}⟩ = ζ(s̄ + s'). This kernel is positive definite on the half-plane Re s > ½ because it is a Gram matrix (PROVED).
- CHECKED, `r30b_fock.py` A1: random 12-point sets with Re s ∈ [0.51, 0.91] give min eigenvalues 0.269, 0.036 and 1.1·10⁻⁵, all > 0.
- On the boundary family Re s = ½ + ε the Gram becomes the Toeplitz kernel ζ(1+2ε+i(t'−t)). It is the Fourier transform of the positive measure Σ n^{−1−2ε} δ_{log n}, so it is PSD by Bochner (PROVED). Its diagonal blows up like 1/(2ε): the pole of ζ at 1 is the Plancherel delta of L²(ℝ₊).
  - CHECKED (A1b): min eigenvalue > 0 for ε = 0.25, 0.1, 0.03, 0.01.

**Where positivity sits (the precise point).** Fock positivity is a property of ζ(w) for Re w > 1, i.e. of Re s > ½ in each slot. Every zero has Re w < 1, so **this positivity never sees a zero** (PROVED: ζ ≠ 0 on Re w ≥ 1).

The inner product that would see the zeros is the Weil form on the critical line, Re s = ½. In it, the primes enter with the **opposite sign**:

W(g) = (pole + archimedean terms) − 2 Σ_n Λ(n) n^{−½} g(log n).

The gas's own positive objects all appear with a minus sign there relative to the vacuum terms:
- the free energy log ζ(β) = Σ Λ(n)/log n · n^{−β}, with coefficients ≥ 0;
- the energy −ζ'/ζ = Σ Λ(n) n^{−β}, with coefficients ≥ 0.

So "nothing goes negative" holds exactly in the gas (β > 1) and fails as a mechanism on the line. There, unitarity would be a competition between the vacuum (pole + Γ) and a positive prime sum, not a sum of positives. This agrees with Round 22 §1: no finite family of primes is passive on its own.

## 2. SL(2,ℝ) generators: what can and cannot exist

**No-go (PROVED).** Suppose L₀ = log N is the L₀ of an sl(2) triple on H_F, with [L₀, L_{±1}] = ∓L_{±1}.
1. Take matrix elements: ⟨m|L₁|n⟩(log m − log n + 1) = 0.
2. That needs m/n = e^{−1}, which is impossible because e is irrational. So L₁ = L₋₁ = 0.
3. Then [L₁, L₋₁] = 2L₀ forces L₀ = 0, which contradicts L₀|n⟩ = log n |n⟩ ≠ 0 for n ≥ 2.

Illustration (A2): the closest |log(m/n) − 1| over n ≤ 2000 is 4.05·10⁻⁸, at 2721/1001, and it is never 0.

The same argument works for L₀ = c·log N. Integer-spaced ladders exist only when e^{1/c} is a rational q, and then only along chains n, qn, q²n, ….
- **Consequence:** the conformal weight Δ cannot be an energy eigenvalue of the gas.
- With s = Δ, the weight is the **complex inverse temperature** β = Δ, and the zeros of ζ are the Lee–Yang/Fisher zeros of Z(β). This reading is known (Julia 1990; Bost–Connes 1995; Knauf's number-theoretic spin chain).

**What does exist.**
- (a) **Dilations** U_a: the ax+b torus. On the Mellin side they act by n^{−s} ↦ a^{s}n^{−s}, with generator D. The principal series Re Δ = ½ is where D is self-adjoint on L²(ℝ₊).
- (b) **Per-prime metaplectic oscillator ("primes as spinors"; PROVED).** On the Fock space of a single prime (occupation N = v_p(n)), set K₀ = (N+½)/2, K₊ = a†²/2, K₋ = a²/2.
  - These close sl(2) (A3: commutator errors 2·10⁻¹³ and 3·10⁻¹³).
  - The rep is unitary, with lowest weights **¼ (even N, the squares p^{2k}) and ¾ (odd N)**.
  - These weights are not in ½ℤ, so this is a genuine representation of the double cover Mp(2,ℝ) and not of SL(2,ℝ). It is the weight-½ (Jacobi θ) representation: spinor-like.
  - **¼ + ¾ = 1: the two sectors are shadows Δ ↔ 1−Δ with the same Casimir Δ(Δ−1) = −3/16** (A3: −0.1875 on both sectors to 2·10⁻¹³).
  - The Euler factor is the trace of the oscillator semigroup, split by sector: (1−x)^{−1} = (1−x²)^{−1} + x(1−x²)^{−1}. Globally, ζ(β) = ζ(2β)·Σ_{n squarefree} n^{−β} (CHECKED A3 at β = 2.3, agreeing to 6·10⁻⁸ with the n ≤ 2·10⁵ truncation).
  - Riemann's θ (n²) is the all-even (¼) sector. This is where the classical "factor 2" (s/2 in Γ(s/2)) lives; Δ = s with no factor 2 is the whole gas. That part is HEURISTIC as a dictionary; the identities are PROVED.
  - The weights log p are linearly independent over ℚ, so the per-prime Mp(2,ℝ)'s do **not** assemble into a diagonal SL(2) with L₀ = log N. This is the no-go again (PROVED).
- (c) **The Weyl element (shadow).** This is not intrinsic to the gas: see §3 and the Beurling control in §5.

## 3. The bridge: Planck ↔ Matsubara = the Weyl element (PROVED, classical; CHECKED)

The Planck factor has the partial-fraction (Matsubara/KMS) form
1/(e^ω−1) = 1/ω − ½ + Σ_{k≥1} 2ω/(ω² + (2πk)²).

Take the Mellin transform over −1 < Re s < 0:
- the **thermal side** (occupation n ≥ 1) gives Γ(s)ζ(s);
- the **Matsubara side** (frequencies k ≥ 1) gives π(2π)^{s−1}ζ(1−s)/cos(πs/2).

Setting them equal is the functional equation. CHECKED (A4) at s = −0.5, −0.3+4i, −0.7+11i: the Mellin quadrature, the thermal closed form and the Matsubara closed form agree to ≤ 10⁻¹⁵ relative in the first two cases and 10⁻¹⁰ in the third.

**Reading.** The Weyl element Δ ↔ 1−Δ is the KMS duality. It swaps two copies of the positive integers: occupation numbers n and Matsubara numbers k. The −½ = ζ(0) is the zero-point term, ½coth(ω/2) − 1/(e^ω−1).

**The shadow intertwiner is unitary exactly on the principal series.** Write ζ(Δ) = χ(Δ)ζ(1−Δ). Then |χ(½+it)| = 1, while off the line |χ(σ+it)| ≈ (t/2π)^{½−σ} (A5, matching to 10⁻⁵ at t = 14, 50, 200; e.g. |χ(0.75+200i)| = 0.421).
- Hence: zero quartet {ρ, 1−ρ, ρ̄, 1−ρ̄} collapses to a pair ⟺ ρ̄ = 1−ρ ⟺ Re ρ = ½ ⟺ the zero sits where the intertwiner is an isometry. This is Claude's 1D-CFT reading, made precise (PROVED, elementary).
- Caution (CHECKED A5): for F_θ = ξ(s−θ)ξ(s+θ), χ_θ(Δ) = χ(Δ−θ)χ(Δ+θ) is **also** exactly unitary on the line, yet every zero of F_θ is off the line. A unitary intertwiner on the line does not put the zeros there.

## 4. Unitarity of the Casimir: equivalent to RH, with the exact obstruction

**Dictionary (PROVED, elementary).** A zero ρ = ½ + iγ is a weight Δ = ρ with Casimir C_ρ = Δ(1−Δ) = ¼ + γ².
- It is unitary principal series iff C_ρ ≥ ¼ is real, iff γ is real.
- The complementary series (real Δ ∈ (0,1)) is excluded for ζ by GPPVerify F9 `xiE_real_ne_zero`.
- An off-line quartet gives a conjugate pair of non-real Casimirs, i.e. a **ghost pair** in a Krein space. By Pontryagin (1944), a self-adjoint operator in Π_κ has at most κ eigenvalues in the open upper half-plane.

**The unitarity Gram of the scaling flow.** Let D(L) = Σ_ρ (1 − cos γL)/γ² = Σ_ρ (1 − cos γL)/(C_ρ − ¼).
- D is computed from the integers alone by the explicit formula. It is the same D as GPPVerify `GppHeartLaplace.Dfun`, whose Laplace transform is ξ'/ξ(½+z)/z².
- G_ij = D(L_i) + D(L_j) − D(|L_i−L_j|) is the Gram of the vectors (e^{iγL}−1)/γ in L²(zeros).
- G is PSD on every grid ⟺ every C_ρ ≥ ¼ ⟺ RH (Schoenberg; Round 23).
- κ(L) = the number of negative eigenvalues on [0, L] is the resolved ghost content. Bombieri (2000) shows that, if only finitely many zeros are off the line, large truncations of the Weil form have exactly one negative eigenvalue per off-line conjugate pair.

**Results** (`r30b_casimir.py`, X = 10⁶, L ≤ 13.82, grids h = 0.25, 0.1, 0.05; NUMERICAL):

| model | κ at L=6 / 9 / 13.8 (h=0.05) | min eig (h=0.05) | D(12) |
|---|---|---|---|
| ZETA | 0 / 0 / 0 (all grids, all L) | +0.0029 | +0.029 |
| FTH (θ=¼) | 13 / 37 / 87 | −18.95 | −0.249 |
| EPS (x²+5y²) | 4 / 11 / 32 | −122.0 | −0.305 |
| DH | 0 / 2 / 6 | −0.237 | +0.199 |

The κ counts grow as the grid is refined, so the absolute numbers depend on the grid. They are lower bounds on ghost content, not literal pair counts. Each off-line pair contributes rank ≤ 4 to G.

New observation (NUMERICAL): for DH the diagonal D(L) stays positive up to L = 13.8 (D ≈ 0.18–0.26). So the full Gram detects the DH ghost (κ > 0 from L ≈ 8) well before the diagonal does. The DH off-line zero 0.808+85.70i contributes only ~e^{0.31L}/γ² ≈ 0.01 to D at L = 13.8.

**Exactly where this reduces to Hilbert–Pólya/RH.**
1. The positive inner product that the gas supplies lives at Re w > 1.
2. The inner product needed for a unitary SL(2) on the principal series is the Weil form.
3. Its positivity, κ = 0, is RH (Weil 1952; Yoshida 1992; Bombieri 2000; Connes 1999).

Nothing in §§1–3 produces it. The gas gives the Hilbert space; the KMS duality gives the Weyl element with a unitary intertwiner; the Casimir's reality is the open statement.

**Unconditional positivity on a subspace (known).** The Weil form is positive on test functions supported in a short interval [−t, t]. Yoshida (1992) verified t = (log 2)/2; Bombieri (2000) gave a reproof. In SL(2) words: unitarity holds on the short-dilation subspace. Our ζ Gram is PSD up to L = 13.8 (NUMERICAL, not a proof).

## 5. Controls: each fake drops exactly one of Daniel's principles (CHECKED, structural)

| principle | ZETA | DH | EPS | BEUR | FTH |
|---|---|---|---|---|---|
| positive integers only (a_n ≥ 0; a positive Fock gas exists) | ✓ | ✗ (a_n ∈ {1, κ, −κ, −1}: signed, no positive gas) | ✓ (a_n = r_Q(n)/2) | ✓ | ✓ |
| free primes (Euler product: Fock = ⊗_p) | ✓ | ✗ | ✗ (degeneracies, not free) | ✓ | ✓ |
| shadow/Weyl element (FE; KMS duality) | ✓ | ✓ | ✓ | ✗ (no FE: only the ax+b group, no SL(2)) | ✓ |
| each integer counted once, "massless" energy log p | ✓ | — | — | ✓ (on its own integers) | ✗ (two species per prime with energies (1±θ) log p; a_p = p^θ + p^{−θ}) |
| unitary shadow intertwiner on Re Δ = ½ | ✓ | ✓ | ✓ | n/a | ✓ (A5) |
| Casimir unitary, κ = 0 (§4) | ✓ (to L = 13.8) | ✗ | ✗ | (no Casimir) | ✗ |
| window functional fires (§6) | ✓ | ✓ | ✓ | n/a | ✗ never |

Conclusions.
- **F_θ is the decisive control** (PROVED).
  - It has a positive Fock gas, an Euler product, an FE and a unitary intertwiner on the line.
  - Its zeros are ρ ± θ. Hardy's theorem gives infinitely many on-line ζ zeros, so F_θ has **unconditionally** infinitely many zeros off its line.
  - So no proof of RH can use only positivity + Euler product + FE + unitary intertwiner. It must use "counted once" (multiplicity one, a_p = 1) quantitatively. This matches the Round 29 finding that multiplicity one is the missing amplifier.
- "All four principles" essentially characterizes ζ (Kaczorowski–Perelli: degree-1 Selberg class). The principles name ζ; they do not by themselves prove anything about it (HEURISTIC).

## 6. The Mazáč–Paulos window mechanism on the principal series

**MP (AFB II, abstract and §5).** Crossing, plus unitarity a_O ≥ 0, plus Regge boundedness, force at least one primary in every window [2Δ_φ+2n−2, 2Δ_φ+2n+2] for large n. The ingredients are sin²(π(Δ−Δ_n)/2)-type functionals with double zeros at the generalized-free dimensions.

**Transplant.** If Λ_M(Δ) = Λ_M(1−Δ) has real coefficients, then Z_M(t) = Λ_M(½+it)/|Q^sG(s)| is real. That is the shadow combined with reality: the principal-series "unitarity" that holds unconditionally.
- Use the positive window k_T(t) = cos²(π(t−T)/2H) on |t−T| < H and set R = |∫Z k_T| / ∫|Z| k_T.
- **R < 1 ⇒ an odd-order zero ON the line in the window** (PROVED logic; this is Hardy 1914 / Hardy–Littlewood 1921). The role of MP's a_O ≥ 0 is played by the positive measure |Z| dt.
- Here no positivity of the spectrum is assumed.

`r30b_window.py`, t ∈ [10, 400], step 0.05, mpmath L-values only (NUMERICAL). The EPS decomposition ½(ζL_{−20} + L_{−4}L₅) was CHECKED against r_Q/2 for n ≤ 300 (error 0).

| H (half-width) | ZETA | FTH | EPS | DH |
|---|---|---|---|---|
| 0.5 | 350/779 | 0/779 | 574/779 | 486/779 |
| 1 | 312/389 | 0/389 | 370/389 | 364/389 |
| 2 | 190/194 (misses at T ≤ 46: sparse low zeros) | 0/194 | 194/194 | 194/194 |
| 4 | 96/96 | 0/96 | 96/96 | 96/96 |

Grid sign changes on [10, 400]: ZETA 202 (matches N(400) − N(10) ≈ 201.6 from Riemann–von Mangoldt), EPS 424, DH 285, **FTH 0**.

- **F_θ is killed, unconditionally (PROVED):** Z_θ(t) = |Λ_ζ(½−θ+it)|²/|G| ≥ 0, so it has no sign change ever. Numerically min Z_θ = 0.0397 > 0 on the grid.
- **DH and EPS pass.** So the window rule is strictly weaker than RH (CHECKED: both have known off-line zeros; the DH ghost is near γ ≈ 85.7 and the EPS ghost near 15.67, as in Rounds 17 and 22).
- **Status for ζ.** Windows [T, T+T^{1/4+ε}] are unconditional (Hardy–Littlewood), and [T, T+T^{27/82+ε}] holds for odd-order zeros (Karatsuba). **Bounded windows** — the true MP analogue, width fixed — are OPEN unconditionally. They hold under RH (Littlewood: gaps → 0).
- So "an on-line state in every fixed window" is a precise statement strictly between Hardy and RH (OPEN). Because DH is expected to satisfy it, it cannot imply RH (HEURISTIC).

## 7. Summary of labels and the sharpest next GAP

- PROVED: the Fock space and its positive Gram on Re s > ½ (§1); the no-go for L₀ = log N (§2); the per-prime Mp(2,ℝ) with shadow weights ¼, ¾ and Casimir −3/16 (§2); the Planck↔Matsubara FE (§3, classical); |χ(½+it)| = 1 and quartet collapse ⟺ isometric intertwiner (§3); F_θ off-line unconditionally and window-dead (§§5–6).
- CHECKED: A1–A5 numerics; the EPS decomposition; the structural table (§5).
- NUMERICAL: κ(L) for the four models (§4); window certification (§6).
- HEURISTIC: θ = the even (¼) sector as the origin of the s/2; "four principles ≈ ζ".
- OPEN: bounded-window on-line zeros for ζ; the integer floor (below). **RH is OPEN.**

**Sharpest next GAP (GAP-30b).**
1. Full unitarity is the PSD of G_ij = D(L_i) + D(L_j) − D(|L_i−L_j|). Its **diagonal** is G_ii = 2D(L_i).
2. F13S (GPPVerify #653, merged as `c81c031`; Codex's fourth-power reduction, Turbo's audit, on Alpöge's window architecture via F13W #651) proves in Lean, with standard axioms only, that the diagonal floor **D(4 log k) ≥ −B for all k ∈ [K, ⌈K^{3/2}⌉], for arbitrarily large K,** implies ζ ≠ 0 on Re s > ½.
3. Conversely, RH gives D ≥ 0.

So **diagonal (1×1-minor) unitarity at the integer points k⁴ is equivalent to full unitarity.** The GAP is to prove that one-sided integer floor. To do so it must use multiplicity one (a_p = 1), because F_θ violates it (D_FTH(12) = −0.249 < 0) while having every other property in §5. This is the same GAP as R29.3 / F13S.5, now located exactly in the SL(2) picture. It is equivalent to RH as a whole, but it is the smallest checkable form found so far.

**Lean.** No module was built this round. Possible later candidates, proposals only, each to be claimed in the conversation file first:
- (L1) the sl(2) no-go for L₀ = log N on ℓ²(ℕ) (elementary, needs irrationality of e);
- (L2) |χ(½+it)| = 1 together with quartet collapse ⟺ Re ρ = ½.
Neither is on the F13 critical path.

## 8. Literature and attribution
- Primon gas: B. Julia (1990), "Statistical theory of numbers"; D. Spector (1990), supersymmetric primon gas and Tr(−1)^F n^{−s} = 1/ζ; Bakas–Bowick (1991); Bost–Connes (1995), the KMS phase transition at β = 1.
- Spectral approaches: Connes (1999), trace formula and absorption spectrum; Berry–Keating (1999), H = xp; Sierra (2008–2019), xp variants; Hilbert–Pólya.
- Weil positivity and its negative index: Weil (1952); Yoshida (1992); Bombieri (2000, Rend. Lincei 11, 183–233), one negative eigenvalue per off-line pair in the finitely-many case; Pontryagin (1944).
- Physics restatements: G. Remmen (2021, PRL 127, 241602), a tree amplitude with masses at the zeros, where real masses ⟺ RH; N. Benjamin & C.-H. Chang (2022, JHEP 11, 143), scalar modular bootstrap whose crossing equation contains the zeros. Our §4 dictionary is of the same type, a restatement. Our §§2, 5 and 6 (no-go, metaplectic shadow pair, four-principle table, MP↔Hardy window with the F_θ kill) we have not seen stated this way; the ingredients are classical.
- Mazáč–Paulos AFB I/II (2018/2019) supplied the window mechanism. Hardy (1914), Hardy–Littlewood (1921) and Karatsuba (1984) supplied the zeta windows.
- Project credit: Alpöge's architecture (`zeroFree_of_rpow_window`) underlies F13W/F13S. Codex did the fourth-power reduction and Turbo the audits. Claude Code contributed the 1D-CFT reading (quartet → pair).

## Files
- `r30b_fock.py` → `r30b_fock.json`: §§1–3 (A1–A5).
- `r30b_casimir.py` → `r30b_casimir.json`: §4. It imports `r22_fejer.py` and `r22_models.py`, copied from round 23.
- `r30b_window.py` → `r30b_window.json`, `r30b_window_Z.npz`: §6.
- `r30b_figs.py` → `fig1_window_Z.png`, `fig2_pontryagin_kappa.png`.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
