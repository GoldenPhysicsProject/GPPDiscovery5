# A CFT Davenport–Heilbronn: off-line zeros in the spectral zeta of unitary, modular-invariant Narain CFTs. RH is OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Goddess Of Celestial CFT, 2026-10-22.

## Labels
- **PROVED** (citation, or a short derivation from cited facts);
- **CERTIFIED** (Arb ball arithmetic via python-flint 0.9.0; rigorous enclosures, argument principle);
- **NUMERICAL**;
- **CONJECTURE**.

Numerics are evidence, never proof.

## 1. Statement
For a positive-definite integral binary quadratic form Q = [a,b,c] with D = b² − 4ac < 0:
- Let τ_Q = (−b + i√|D|)/(2a), Z_Q(s) = Σ'_{(m,n)} Q(m,n)^{−s}, and Λ_Q(s) = (√|D|/2π)^s Γ(s) Z_Q(s).
- **PROVED (Epstein 1903; e.g. Lee, Forum Math. 2014, doi:10.1515/forum-2012-0057).** Λ_Q(s) = Λ_Q(1−s).

**The CFT.** Take the c = 2 Narain T² CFT with complex structure σ = τ_Q, B = 0 and Kähler modulus ρ = iρ₂.
- Lattice: metric G = (ρ₂/τ₂)[[1,τ₁],[τ₁,|τ|²]], with Δ = ½mᵀG⁻¹m + ½wᵀGw and J = m·w.
- Partition function: Z(τ) = |η(τ)|⁻⁴ Σ_{(m,w)∈Z⁴} q^{h} q̄^{h̄}.

**Facts.**
1. **PROVED (Narain, PLB 169 (1986) 41; Polchinski vol. 1 §8.4).** The lattice is even self-dual of signature (2,2), so Z is SL(2,Z)-invariant.
   - The multiplicities are (lattice counts) × (coefficients of 1/η(τ)² and its conjugate), all ≥ 0.
   - The vacuum is unique, because Δ is positive definite on (m,w) ≠ 0.
   - **NUMERICAL** (`narain.py` → `narain_checks.json`, M = 5, all 12 forms): |Z(−1/τ) − Z(τ)|/|Z| ranges from 1e−9 to 1e−30 (truncation-limited for skinny tori); |Z(τ+1) − Z(τ)|/|Z| ≤ 1e−30; vacuum multiplicity 1; gap > 0.
2. **PROVED (BCFMP arXiv:2107.10744 eq. 3.25, plus 2ζ(2s)E_s(τ_Q) = (√|D|/2)^s Z_Q(s)).** The scalar U(1)² primary spectral zeta is
   S(s) := Σ'_{scalar primaries} Δ^{−s} = |D|^{s/2} E_s(iρ₂) Z_Q(s).
   - The pure-momentum part alone is (ρ₂√|D|)^s Z_Q(s).
   - **NUMERICAL** (`scalar_zeta_identity.log`, `momentum_identity.py`): both identities agree at s = 3 to about 1e−7 / 1e−9 (lattice truncation), for ρ₂ ∈ {1, 1.3} and {1, 1.7}.
3. **PROVED (Davenport–Heilbronn, J. London Math. Soc. 11 (1936) 181–185, 307–312; anticipated by Potter–Titchmarsh 1935; Voronin; Lee 2014 Thm 1.3).** When h(D) > 1:
   - Z_Q has infinitely many zeros in Re s > 1 (DH);
   - it has ≍ T zeros in every strip ½ < σ₁ < Re s < σ₂ ≤ 1 (Voronin, Lee).
   - Hence S(s) has zeros off Re s = ½ as well; E_s(iρ₂) has no zeros or poles there to cancel them, for generic ρ₂ (CONJECTURE for specific ρ₂; it holds at ρ₂ = 1, where E_s(i) ∝ ζ(s)L(s,χ₋₄)/ζ(2s)).
4. **CERTIFIED (this work).** Explicit off-line zeros: §3 and ZERO_TABLES.md.

**What it shows.**
- Unitarity, a unique vacuum, discreteness and full modular invariance do **not** force the zeros of a CFT's (scalar) spectral zeta onto Re s = ½. This is a genuine CFT, not just a "fake Z": the Narain OPE is consistent.
- The mechanism is exactly DH's: Z_{x²+5y²} = ζ(s)L(s,χ₋₂₀) + L(s,χ₋₄)L(s,χ₅), a sum of two Euler products, with χ₅ of DH's conductor 5 (genus-character decomposition; PROVED; checked to 27 digits in `zq_L.py`).

## 2. Prior art (checked by search and full-text grep, 2026-10-22)
Off-line zeros of Epstein zeta functions are classical. The Epstein side of this note is not new; only the CFT packaging and the certificates are.
- **H. S. A. Potter and E. C. Titchmarsh, "The zeros of Epstein's zeta-functions", Proc. London Math. Soc. (2) 39 (1935) 372–384, doi:10.1112/plms/s2-39.1.372.** They prove infinitely many zeros on Re s = ½ for integral forms, and they remark that the "Riemann hypothesis" for Z_Q is probably false. The off-line zeros are anticipated here.
- **H. Davenport and H. Heilbronn, "On the zeros of certain Dirichlet series" I, II, J. London Math. Soc. 11 (1936) 181–185 and 307–312.** They prove that Z_Q has infinitely many zeros in Re s > 1 when h(D) > 1. This is the theorem we use, §1 item 3.
- **E. Bombieri and D. A. Hejhal, "On the distribution of zeros of linear combinations of Euler products", Duke Math. J. 80 (1995) 821–862.** Under GRH plus a spacing hypothesis, almost all zeros of such combinations (Z_Q included) lie on the line. The off-line zeros are a density-zero set.
- Voronin (1970s): ≫ T zeros in every strip ½ < σ₁ < σ < σ₂ ≤ 1. Lee, Forum Math. 26 (2014): asymptotic counts.
- E. Bombieri and J. Mueller, "On the zeros of certain Epstein zeta functions", Forum Math. 20 (2008) 359–385, doi:10.1515/FORUM.2008.018. Class number 2 (x²+5y², 2x²+2xy+3y²): upper and lower bounds on how fast zeros approach the boundary of the zero-free half-plane Re s > σ(Q). Our CERTIFIED zero 2.0794 + 7.6289i of 2x²+2xy+3y² gives σ(Q) ≥ 2.0794, and `dominant_abscissa.py` gives σ(Q) ≤ 2.4102; see §3.
- Travěnec–Šamaj, Appl. Math. Comput. 2022 (off-critical zeros of hypercubic Epstein zetas).
- Nakamura–Pańkowski arXiv:1212.5890 (spectral zetas with off-line zeros).
- CFT side:
  - BCFMP 2107.10744 (Narain spectral decomposition, eq. 3.25);
  - Benjamin–Chang 2208.02259 (scalar crossing ↔ ζ zeros);
  - Angelantonj–Florakis–Pioline 1110.5318 (constrained Epstein series);
  - "Fake Z" (JHEP 06 (2023) 043, arXiv:2211.15699): modular-invariant, positive, integral but non-CFT partition functions.
- **BCFMP 2107.10744 does not state off-line zeros.** A grep of its full text (2026-10-22) finds no mention of Davenport–Heilbronn, class number, or zeros of the Narain spectral zeta off Re s = ½. Its "critical line" hits are the Eisenstein continuum and the no-poles condition (3.51)–(3.52). We credit BCFMP for the spectral identity only (§1 item 2).
- **Not found:** an explicit statement that Narain-CFT scalar spectral zetas have certified off-line zeros. The novelty claim is limited to this packaging plus the certificates. The underlying mathematics is DH 1936.

## 3. Computation
### Evaluators
- `zq_L.py`: Z_Q from Arb Hurwitz zeta via the genus decomposition. This works for D = −4, −20, −84, whose class groups are 2-torsion.
- `theta.py`: Λ_Q from the incomplete-gamma theta split with a rigorous tail. This works for all forms, at point s.
- The two methods agree to 27 digits (`validate2.py`, plus the check in zq_L).

### Certification (`certify_L.py`)
- Argument principle on box boundaries. Every segment is enclosed in a complex ball, and Z_Q is evaluated on it.
- If the enclosure lies in an open half-plane, the arg change along the segment is the principal Arg(Z(end)/Z(start)); otherwise the segment is bisected.
- The total divided by 2π must enclose exactly one integer.

### No zeros with Re s ≥ 3 (`no_zeros_right.py`)
- **CERTIFIED** for all 12 forms: exact rational sum plus a lattice-count tail bound; the margin B(3) ≤ 0.69 < 1.
- The float eigenvalue in the tail has a margin of at least 0.3.

### Results
- **CERTIFIED:** 81 zeros of Z_Q for D = −20, −84 with 0.52 < Re s < 3 and 1 < Im s < 100, each isolated in a 1e−8 box (`certified_zero_boxes.json`). Examples:
  - x²+5y²: 0.9079933666 + 89.5598279786i.
  - 2x²+2xy+3y²: **2.0794141048 + 7.6288675878i**. This is a zero far inside the region of absolute convergence, as DH predict.
  - 3x²+7y²: 1.2084504268 + 18.4657092717i.
- **CERTIFIED zero counts** over the whole region [0.52,3]×[1,201] for the genus forms: ZERO_TABLES.md (`certified_counts_T201.json`), if present.
- **NUMERICAL** (`scan.py`; float winding plus Newton, Λ evaluated with rigorous point enclosures, but the count is not rigorous): zeros for D = −23 and −56, for example:
  - 2x²+xy+3y²: 1.2454 + 85.668i;
  - 2x²+7y²: 1.1835 + 98.231i.
- **Control, NUMERICAL** (and CERTIFIED count, if present in the table): x²+y² (h = 1, Z = 4ζ(s)L(s,χ₋₄)) shows **no** zeros in [0.52,3]×[1,100], consistent with GRH for ζ·L(χ₋₄).
- **Blind spot:** zeros with ½ < Re s < 0.52 are not searched.

## 4. Battery plug-in: the fake "Λ_Q" (`battery_coeffs.py`)
**Why this pair.** For a principal form, Λ_Q has the **same** archimedean factor Γ(s)(√|D|/2π)^s and the same conductor as the Dedekind zeta ζ_K, K = Q(√D) (PROVED: both are degree-2 objects with this Γ-factor; for D = −20, Z_{x²+5y²} = ζ_K + L(χ₋₄)L(χ₅)).
- So in the Weil explicit formula the archimedean terms of ζ_K and Z_Q coincide **exactly**, and only the prime side differs.
- That makes Z_Q a sharper control than DH, whose Γ-factor differs from ζ's.

**Prime side, NUMERICAL** (float, n ≤ 300). The coefficients c(n) of −Z_Q′/Z_Q:
- are supported off prime powers: n = 6, 14, 21, … for x²+5y²;
- are negative at n = 36, 54, 84, … (8 negatives for x²+5y², 10 for x²+xy+6y², 3 for x²+14y², 0 for x²+21y²).
- For ζ_K, c(n) = Λ_K(n) ≥ 0 on prime powers.

**Proposal for the battery (CONJECTURE that it is informative):**
- add the pair (ζ_K, Z_Q) at D = −20 and run the window Weil form Q_{L,N} on both;
- if the window form is PD for Z_Q at some L where Z_Q has a certified off-line zero of height γ ≲ the window's resolution, the window positivity is blind to the Euler product, the same lesson as DH (nogo_note row 23) but with matched archimedean data;
- the S98 machinery needs only c(n) and the shared Γ-factor.

## 5. Dead ends
- c = 1 compact boson: its scalar spectral zeta is 2ζ(2s)(R^{2s} + R^{−2s}) (BCFMP 3.11). It never produces a DH-type fake.
- The K-Bessel (Chowla–Selberg) Eisenstein evaluator `eps.py`: flint bessel_k balls blow up near the turning point at t ≳ 60. It is kept only as a low-t cross-check.
- Ball-input certification through the theta formula (`certify.py`): the Γ(s) cancellation (|Λ| ~ e^{−πt/2} against O(1) terms) makes ball radii useless beyond t ≈ 20.
  - The fix is the L-decomposition, but it is available only for 2-torsion class groups.
  - D = −23 (Z/3) and −56 (Z/4) therefore stay NUMERICAL. Next: Hecke L-functions of class-group characters via Arb, or a rotated-contour theta (Rubinstein's δ trick).

## 6. Next
- Certify D = −23 and −56 (Hecke L or rotated theta).
- Push the certified counts to T = 500.
- Hand the pair (ζ_K, Z_Q) to the window-Weil seat.
- Test BC footnote 14: does Z_Q ever vanish at (1+ρ)/2?
