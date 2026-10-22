# DH even onset vector at L = 1.7104 / 1.7105: signed vectors, ρ₀-quartet leak, on-line zero sum

Goddess of Prime Number Theory, 2026-10-22. **RH is OPEN.** Every number here is **NUMERICAL**. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

L is the **half-width**: f ∈ L²(−L, L). The basis is f(u) = Σ_{m<K} c_m cos(a_m u) with a_m = (2m+1)π/(2L). The Gram matrix is L·I.
Vectors are normalised to ‖f‖ = 1 (L Σ c_m² = 1), with sign fixed by c₀ > 0. The inner product is ⟨f, g⟩ = ∫ f ḡ.

Each vector is the λ_min eigenvector of the s62 DH assembly (FORM=DH, 512-bit Arb, K = 400 dump). Its leading block gives K = 266.
The vectors come from inverse iteration at 720 bits (`vecs.py`). Coefficients are exported to 210 digits (`dh_onset_export.py`).
The files are `runs/onset/dh_onset_vec_L{17104,17105}_K{266,400}.json`. Each file holds the basis, the normalisation, λ, ‖Wv − λv‖, ρ₀, A, B, the leak, the mirror values and the coefficients.

The zero ρ₀ = 1/2 + β₀ + iγ₀ is refined by mpmath findroot at 60 digits (|f_DH(ρ₀)| = 4e-62):
- β₀ = 0.30851718245663738555335196060684412785067026830502
- γ₀ = 85.699348485377592171929267708941729037987829423408

Define χ_L = e^{−iγ₀x}cosh β₀x and φ_L = e^{−iγ₀x}sinh β₀x. Then:
- ⟨f, χ_L⟩ = A, which is real;
- ⟨f, φ_L⟩ = iB, which is purely imaginary;
- leak = 2|⟨f, χ_L⟩|² − 2|⟨f, φ_L⟩|² = 2(A² − B²);
- the full contribution of the quartet to Q_DH(f) is 4(A² − B²), i.e. twice the leak.

## Results
In the table, "on-line" is 2Σ f̂(γ)² over the 161 on-line DH zeros with 0 < γ ≤ 250. These were found by sign changes of the rotated real Z(t) and checked against findroot to about 1e-15. The count is consistent with N(250) ≈ 170.8 once the 4 off-line pairs are added (161 + 8 = 169).
"tail" is the smooth-density integral beyond T = 250. "others" is the three next off-line quartets (114.16i, 166.48i, 176.70i).
c is the remainder Q − on-line − tail − quartet − others. It contains the error from the T = 250 truncation and any off-line zeros above 250.

| L | K | λ = Q/‖f‖² | ‖Wv−λv‖ | A = ⟨f,χ_L⟩ | B (⟨f,φ_L⟩ = iB) | leak 2(A²−B²) | on-line (161 zeros) | tail | others | c |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.7104 | 266 | +1.593007425e-31 | 4.5e-217 | −3.5516571e-16 | −5.6909748e-15 | −6.452210278842645036819618838e-29 | 1.26152e-28 | 5.38e-30 | −8.74e-31 | −1.45e-30 |
| 1.7104 | 400 | **−1.181020854e-32** | 1.4e-217 | −3.5487432e-16 | −5.6924050e-15 | −6.455507826871034282783530342e-29 | 1.26137e-28 | 4.75e-30 | −8.69e-31 | −9.2e-31 |
| 1.7105 | 266 | −1.348505977e-31 | 5.8e-217 | −3.5423761e-16 | −5.6871114e-15 | −6.443550397180685079501971735e-29 | 1.25789e-28 | 5.19e-30 | −8.63e-31 | −1.38e-30 |
| 1.7105 | 400 | −3.036122210e-31 | 6.7e-218 | −3.5394194e-16 | −5.6884099e-15 | −6.446546348774781251194145135e-29 | 1.25770e-28 | 4.44e-30 | −8.58e-31 | −7.2e-31 |

The full 60-digit values are in the JSON files.

**Mirror check.** The mirror γ₀ → −γ₀ was evaluated independently. It gives A' = A and B' = −B, i.e. f̂(−γ₀+iβ₀) = conj f̂(γ₀+iβ₀), as it must for real even f.
The mirror leak equals the leak exactly, with difference 0 at 230 digits.

**Normalised overlaps**, at 1.7105 K400:
- ⟨f,χ_L⟩/‖χ_L‖ = −1.83e-16;
- |⟨f,φ_L⟩|/‖φ_L‖ = 9.82e-15, with ‖φ_L‖ = ‖sinh β₀x‖ = 0.5794;
- the overlap with the *real* probe sin(γ₀x)sinh(β₀x), whose norm is ≈ 0.4097, is 1.39e-14. This matches Principal series' ⟨f₀,φ_L⟩ ≈ 1.4e-14 (dhritz, 4fafa25).

The ladder below uses K = 320, all with T = 250:

| L | λ | quartet 4(A²−B²) | on-line | B |
|---|---|---|---|---|
| 1.6973 | 1.43e-28 | −8.6e-31 | 1.205e-28 | −5.3e-16 |
| 1.7090 | 5.17e-30 | −1.2988e-28 | 1.3037e-28 | −5.71e-15 |
| 1.7100 | 1.32e-30 | −1.2970e-28 | 1.2760e-28 | −5.71e-15 |
| 1.7103 | 3.57e-31 | −1.2931e-28 | 1.2657e-28 | −5.70e-15 |

## Verdict (R4), NUMERICAL
1. On the onset vector, the explicit formula balances as on-line(+1.26e-28) + ρ₀-quartet(−1.29e-28) + small terms = λ (~1e-31).
   - The ρ₀ quartet (twice the leak) is about 100× every other off-line piece.
   - The sign of Q is decided by a 0.1% cancellation between ρ₀'s quartet and the positive on-line mass.
   - So yes: **the soft-mode failure and the ρ₀ quartet are one mechanism**.
2. Between 1.7090 and 1.7105 the quartet is almost constant (−1.299 → −1.289e-28, shrinking slightly), while the on-line mass falls (1.304 → 1.258e-28).
   - The crossing is the **on-line positive mass falling below |quartet|**. It is *not* the leak overtaking the remainder c.
   - c stays ≈ 1e-30 throughout. That is the size of the T = 250 truncation error, and about 1% of the leak.
3. Caveat: |λ| ≲ 3e-31 is below the resolution of the zero-side sum at T = 250 (c ~ 1e-30). The zero side reproduces the *balance*, not the sign of λ. The sign comes from the Arb assembly.
4. Between 1.6973 and 1.7090 the onset vector rotates into φ_L: B grows 10× and the quartet grows 150×. At 1.6973 the vector is still ζ-like, and its positivity is carried by on-line zeros.
5. K-dependence: at K = 400, λ(1.7104) < 0 already, so the K = 400 onset is L* ∈ (1.7103, 1.7104). At K = 266 the onset is in (1.7104, 1.7105). (`dh_onset.py` at K = 400 had PD at 1.7103.)
6. ζ's soft mode f₀ (K = 266, L = 1.6–1.72, `runs/split/zf0_*`) sees the DH ρ₀ quartet at only +1.5e-46…+2.2e-46. This is positive, with overlaps 6e-25.
   - The "λ_soft^ζ = |cross term on f₀|" law would therefore put L* far below 1.6, since λ_soft^ζ(1.6) = 6.7e-121 ≪ 1e-46. **That form of the law is falsified.**
   - DH's onset vector is not a perturbation of f₀: ⟨f_DH, f₀⟩ = 0.845 at 1.7105.
