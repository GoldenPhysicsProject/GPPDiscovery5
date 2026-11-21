# Alpha search: why is 1/α = 137.035999177?

Daniel asked (t584u): "Try to find it", meaning find α in the numbers. Later he asked for a check of his old GGC theory too.
— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). 2026-11-21.

**Verdict: nothing derives α.** The best golden-ratio formula is a real but modest numerical
coincidence (NUMERICAL). It is not a derivation. RH-style rules apply: no claim is closed.

## 0. Target and what "the value" means

| source | 1/α | 
|---|---|
| CODATA 2022 (recommended) | 137.035999177(21) |
| Fan+ 2023, electron g−2 + QED (PRL 130, 071801) | 137.035999166(15) |
| Morel+ 2020, Rb recoil (Nature 588, 61) | 137.035999206(11) |
| Parker+ 2018, Cs recoil (Science 360, 191) | 137.035999046(27) |

- The Rb and Cs values disagree at 5.5σ. So "within 1σ of CODATA" is a fragile standard. We test at 2.1e-8 (1σ CODATA) and at 1.6e-7 (the Cs–Rb spread).
- **Running.** 137.036 is the Thomson-limit value at momentum transfer q → 0, the far infrared end. At the Z mass, 1/α̂(M_Z) = 127.951(9) (MS-bar; PDG 2023 EW review, Δα_had^(5) = 0.02768(7)).
  - Vacuum polarization screens charge, so 1/α decreases with energy. A "bare" (UV) 1/α is therefore smaller than 137.036, never larger.
  - Any scheme that dresses a bare 137.5 down to 137.036 has the sign of QED screening backwards.
  - There is no "colder" value than α(0). Lab temperature enters only through powers of kT/m_e ≈ 5×10⁻⁸ at 300 K, which is negligible.
  - Consequence for the golden angle: 360/φ² = 137.508 cannot be the zero-temperature or bare α.

## 1. Daniel's own material

**GPPVerify `GppVerify/NumberTheory/DecodingReality.lean` (main c81c031): no α theorem.**
- It proves `Nat.gcd 137 360 = 1` (`cp_phase_fraction_simplified`). The CP phase δ_CKM = 137π/360 is the placeholder `open_cp_phase_from_jacobi_sum : True := trivial`. That is a stub, not a proof.
- `SquarefreeDensityZeta.lean` proves ∏_p(1−p⁻²) = 6/π² and states explicitly that it makes no claim about α.
- In ONON v5.2.1.3 (tex), `conj:cp-alpha` is labelled an *empirical observation*: δ ≈ 137π/360, 0.04%. The same document lists "α_em⁻¹ = 137.036" under "Not yet derived (open problems)".

**GGC (Golden Geometric Coherence), Daniel's Drive, GPP folder, Nov–Dec 2025.**
- Contents: E8 projected onto dual 4-manifolds with golden-ratio weighting. ℏ = 2π/φ, a leakage parameter λ = φ⁻⁴, and closed forms for SM constants.
- **GGC III, Appendix B:** α⁻¹ = (360/φ²)(1 − (π/3)φ⁻¹² − (π/8)φ⁻¹⁶ − (π/668)φ⁻²⁰). Here 360 = 120 charge axes × 3 polarizations, and 668 = 4 × 167, with 167 = 240 − 72 − 1.
- The paper claims the value 137.03599908413 and a relative error of 2.0×10⁻¹².

Our recomputation (`candidates.py`, 50 digits) shows the formula as printed gives **137.0360516351**, not the claimed value:

| truncation | value | off by |
|---|---|---|
| bare 360/φ² | 137.5077640500 | |
| + C1 | 137.0605616195 | |
| + C2 | 137.0360943863 | |
| + C3 (as printed) | **137.0360516351** | +5.25×10⁻⁵, i.e. 3.8×10⁻⁷ relative, 2500σ |

- The printed number is the CODATA 2018 central value with digits appended. It was not computed from the formula. **ERRATUM, GGC III B2.**
- To hit CODATA 2022, the last denominator would have to be 299.947, not 668. Using 300 lands at +0.8σ, but that is a post-hoc fit with no derivation behind it. In the GGC grammar (n_i ≤ 400), the chance rate for a 1σ hit is 1.2% per target (`control_targets.py`).
- **GGC II:** α⁻¹ = 4π/φ² evaluates to 4.80, not 137.036. **ERRATUM**, probably meant 360/φ².
- **Other GGC III closed forms** (`ggc_other_constants.py`):

| claim | GGC says | actual value | observed |
|---|---|---|---|
| v = M_P φ⁻⁸⁰ | 246 GeV | 233.2 GeV | 246.22 GeV |
| m_H = vφ⁻⁴ | 125.1 GeV | 35.9 GeV | 125.2 GeV |
| η_B = φ⁻⁸/(2π) | 2.53×10⁻¹⁰ | 3.4×10⁻³ | ≈6.1×10⁻¹⁰ |

  All three are **ERRATA**. They stay on record.

## 2. Literature (verified sources)

| attempt | 1/α predicted | digits | status |
|---|---|---|---|
| Eddington (1929 → 1930s) | 136, then 137 | 2–3 | rejected: integer; changed to fit the data (Kragh, arXiv 1510.04046) |
| Wyler 1969, volumes of O(5,2)/O(5)×O(2) domains | 137.0360824 | 6 | Robertson PRL 27, 1545: radius set to 1 arbitrarily. Pease 1977: wrong Poisson-kernel coefficient. Now +4000σ off |
| 4π³+π²+π (folklore) | 137.0363038 | 5–6 | numerology |
| Atiyah 2018 ("Todd function") | no reproducible number | – | the Todd function is provably constant from its own axioms (open-mapping / identity-theorem arguments); physics objection that α runs (Carroll, Hossenfelder) |
| Singh 2022 (EPJ Plus 137:664), exceptional Jordan algebra | α = 0.00729713 (≈137.040) before a correction tuned at 144.53 GeV | 4–5 | correction scale fitted |
| GUT running | α_GUT and M_GUT are inputs; sin²θ_W is predicted given α, α_s | 0 | moves the question to α_GUT |
| String moduli | set by vacuum choice (landscape) | 0 | no unique prediction |
| Anthropic | windows such as 1/α ≈ 145 ± 50 (stellar Fe end-point) and 145 ± 9 (radiogenic heating) (Sandora, arXiv 1604.03151); the GUT window 1/180 < α < 1/85 (Barrow–Tipler) | order of magnitude | explains the ballpark, not the digits |
| Heyrovská: 360/φ² − 2/φ³ | 137.0356281 | 5.6 | numerology |
| 3-term golden: 360/φ² − 2/φ³ + (3φ)⁻⁵ (attributed to Sherbon; Pellis SSRN 4160769) | **137.0359991648** | **10.05** | see §3 |

## 3. Numerical search with look-elsewhere control

Scripts: `candidates.py`, `search.py`, `control_targets.py`, `control_A_fine.py`, `search_golden4.py`, `grammar_of_grammars.py`. JSON outputs sit next to them.

**Method.**
- Enumerate every expression in a grammar.
- Keep distinct values only. Golden identities such as φ⁻¹ = φ⁻² + φ⁻³ collapse many rewrites into one.
- Count hits at the real target. Then count hits at 2,000–200,000 decoy targets drawn uniformly near 137.036. The decoy rate is the honest chance rate.

| grammar | size | 1σ hits at real T | chance rate per target |
|---|---|---|---|
| A: 360φ⁻² ± aφ⁻ⁿ ± (bφ)⁻ᵐ or bφ⁻ᵐ (a,b,n,m ≤ 9) | 4.96×10⁴ | 1 (the 3-term golden formula) | 2.6–3.4×10⁻⁴ |
| B: GGC grammar, π/n_i with n_i ≤ 400 | 6.4×10⁷ | 1 (n = 3, 8, 300, not 668) | 1.2×10⁻² |
| C: r·xᵃ·yᵇ over 13 framework constants (π, φ, e, √2, √3, √5, √7, ζ(3), ζ(5), Catalan = L(2,χ₋₄), L(2,χ₋₃), ln 2, γ), r = p/q, p ≤ 400, q ≤ 12 | 2.4×10⁷ | **0** | 1×10⁻³ |
| D: 137 ± r·xᵃ·yᵇ, p, q ≤ 12 | 1.4×10⁶ | **0** | 4×10⁻² |
| 216 sibling grammars shaped like A (leads 360/φ², 137, 4π³+π²+π; base constants φ, π, e, √2, √3, √5) | 216 × ~5×10⁴ | 3 total | 0.32 expected |

Within the 216 sibling grammars:
- **360/φ² lead (72 grammars):** 1 hit observed, 0.0063 expected, P ≈ 0.6% (about 2.5σ).
- **4π³+π²+π lead:** 2 hits, 0.25 expected, P ≈ 2.7%. The hits are `−9π⁻⁹ − (6√2)⁻⁶` and `−e⁻⁸ + (6√5)⁻⁴`.
- **Adding a 4th golden term:** 21 distinct 1σ formulas appear. All of them are the 3-term hit plus negligible tails.

**Top candidates** (digits = −log₁₀ of relative error vs CODATA 2022):

| # | formula | value | digits | pull vs C22 / Fan / Rb / Cs | look-elsewhere |
|---|---|---|---|---|---|
| 1 | 360/φ² − 2/φ³ + (3φ)⁻⁵ = (176410 − 88447φ)/243 | 137.0359991648 | 10.05 | −0.6 / −0.1 / −3.8 / +4.4 | 3×10⁻⁴ in its grammar; ≈0.6% across 72 golden siblings; not significant once grammar choice and post-hoc selection are counted |
| 2 | GGC with 300 for 668 | 137.0359991937 | 9.91 | +0.8 / +1.9 / −1.1 / +5.5 | P ≈ 1.2% in its grammar; post-hoc |
| 3 | (4π³+π²+π) − e⁻⁸ + (6√5)⁻⁴ | within 1σ | ~10 | – | 2.7% per lead; shows how cheap such hits are |
| 4 | Wyler | 137.0360824 | 6.2 | +4000σ | historical |
| 5 | GGC as printed | 137.0360516 | 6.4 | +2500σ | erratum |

**Why #1 is still not a derivation:**
- **(a) Units.** 360 is degrees, a Babylonian convention. In radians, the golden angle is 2π/φ², so the formula hides a factor of 180/π.
- **(b) Field.** The value lies in ℚ(√5), with height about 1.8×10⁵ and denominator 243. Roughly 200 elements of ℚ(√5) of equal or smaller height fall within 1σ. Only the *short spelling* is rare.
- **(c) No running and no mechanism.** α(0) is fixed by QED vacuum polarization from e, μ, τ, and hadrons. A pure-number formula would need all of these masses to conspire with it.
- **(d) It is not a prediction.** It sat at +3.8σ from CODATA 2018 (137.035999084). It agrees only with the g−2 route and fails Rb at 3.8σ and Cs at 4.4σ. The Rb/Cs tension decides it: a future resolution toward Rb would kill it.
- **(e) Post-hoc choices.** Both the lead (chosen because it is near 137) and the term (3φ)⁻⁵ came from search. Two other unrelated families also produce 1σ hits.

## 4. Labels

| item | label |
|---|---|
| DecodingReality.lean: no α content; the CP-phase item is a `True` stub | CHECKED |
| GGC III B2 evaluates to 137.0360516, not 137.03599908413 | CHECKED (ERRATUM) |
| GGC II 4π/φ² = 4.80 | CHECKED (ERRATUM) |
| GGC v, m_H, η_B formulas do not give the stated numbers | CHECKED (ERRATA) |
| 360/φ² − 2/φ³ + (3φ)⁻⁵ matches CODATA 2022 at 0.6σ; chance rate ≈ 3×10⁻⁴ per grammar, ≈0.6% across 72 golden siblings | NUMERICAL |
| "α comes from the golden angle" | HEURISTIC (no mechanism; units and running objections) |
| A structural derivation of α(0) | OPEN; no candidate found |
| Daniel-framework arithmetic constants (ζ, L-values, Gauss-sum moduli √p, Casimir rationals) produce no 1σ hit in 2.4×10⁷ monomials | NUMERICAL (null) |

## 5. "On the Origin of Spacetime": every version (Daniel t587u)

Script: `origin_of_spacetime_versions.py`; output in `origin_of_spacetime_versions.json`. Daniel's drafts are not copied here. The formulas were transcribed from the Drive PDFs and from `AppendixG.tex` inside `Origin of Spacetime-Newest.zip`. The drafts write "log" without a base, so both ln and log₁₀ were evaluated, along with the "1/α vs α" alternative readings.

**Versions found** (Drive, GPP folder; dates from the title page or file time):

| date (2025) | file(s) | what it says about α |
|---|---|---|
| May 14 | On The Origin Of Spacetime V2.docx | qualitative only: "constants from projection". No formula |
| May 20–22 | Information_and_Holographic_Spacetime (3)–(9), titled "On the Origin of Spacetime(s)" | table entry α(θ) ∼ 1/κ_U(1)(θ), with the observed 1/137.035999 written beside it. No computable formula |
| May 24 | OnTheOriginOfSpacetime.pdf, OnTheOriginOfSpacetimeGG.pdf | G.4: α ∼ 1/(θϕ) at θ = φ, said to give "the observed value within projection uncertainty" |
| May 25 | MPNP (7), (8), _Final_Copy; Information_and_Holographic_Spacetime (13), (16), (24)–(27) | G.2: α(θ) = φ⁻⁶ · θ^1.472 / (1 + 1.78 log(1/θ)), "at θ = φ yields α⁻¹ ≈ 137.036". Also gives sin²θ_W and α_s formulas |
| May 27 | MPNP (11) | formula boxes are empty in the PDF: "Using : Matches the observed value to high precision" |
| May 27 → Jun 19 | MPNP (12), (16), (19)–(22), (28); _complete; MPNP__; Revised; v3; Newest.zip (Jun 11) | G.3: α⁻¹ = Σ_{k≥2} κ_α / (φᵏ[1 + (θ − φᵏ)²]) at θ = φ⁹, with "fitted κ_α = 19.56" giving α⁻¹ ≈ 137.036. Identical in every later draft |
| Jun 7–11 | MPNP (23)–(27), ~260 KB | actually the Yang–Mills mass-gap paper; no α |
| Jun 15 | On_the_dimensionality_of_spacetime.pdf | no α |
| Nov–Dec 2025 | GGC III / II (see §1) | 360/φ² series (gives 137.0360516); 4π/φ² (gives 4.80) |

**Formula by formula** (50 digits; CODATA 2022 = 137.035999177(21), Rb = …206(11), Cs = …046(27)):

| version | formula as printed | printed | computed 1/α | off by |
|---|---|---|---|---|
| May 24 | α ∼ 1/(θφ), θ = φ | "observed" | 2.618 (φ²) | −134.4 |
| May 24, alt. (ϕ = golden angle 137.5°) | θ · 360/φ² | – | 222.49 | +85.5 |
| May 25, ln | 1/[φ⁻⁶ θ^1.472/(1 + 1.78 ln(1/θ))], θ = φ | 137.036 | 1.2676 | −135.8 |
| May 25, log₁₀ | same with log₁₀ | 137.036 | 5.5495 | −131.5 |
| May 25, read as 1/α directly | – | – | 0.789 (ln) / 0.180 (log₁₀) | −136 |
| May 27 – Jun 19 | 19.56 · Σ_{k≥2} 1/(φᵏ[1 + (φ⁹ − φᵏ)²]) | 137.036 | **0.26182** | −136.8 |
| same, read as α | 1/(19.56 · S) | – | 3.8194 | −133.2 |
| May 27 (11) | (empty) | "high precision" | n/a | – |

- Every result is off by 85–137 in absolute terms, so the pulls against CODATA, Rb, and Cs are all 10⁹–10¹⁰σ.
- S(φ⁹) = 0.0133856279, so the κ_α that would actually give CODATA is **10237.547343**, not 19.56.
  - Rounding κ_α = 19.56 to its printed precision moves the result only within [0.26176, 0.26189].
  - The value 19.56 is called "fitted", but it is not a fit to 137.036 under any reading I could construct.

**Other constants in the May 25 draft:**

| quantity | printed | computed (ln reading) | computed (log₁₀ reading) |
|---|---|---|---|
| sin²θ_W at θ = φ | 0.231 | 0.995 | 0.94 |
| α_s at θ = φ | "≈ 1" | 0.235 | 0.336 |

From May 27 on, the drafts drop these formulas. sin²θ_W is quoted from observation (0.231), and α_s comes from the standard one-loop QCD formula. G "derived" from M_P = √(ħc/G) is circular.

**Look-elsewhere.** No version lands near 137.036, so there is no hit to discount. Two structural points still matter:
- Both computable versions contain fitted continuous parameters (1.472 and 1.78; κ_α).
- With one free continuous parameter, the formula reaches *any* target (κ = 10237.5 reaches every one of CODATA, Rb, and Cs). A correct match would therefore carry zero evidential weight (P = 1). Compare grammar A in §3, where no parameter is fitted and the chance rate is 3×10⁻⁴.

**Verdict.**
- No version of "On the Origin of Spacetime" gets α right.
- In every version where a number is printed, it is the observed 137.036 written in. The formula beside it does not produce it: computed values range from 0.18 to 222.
- Across all of Daniel's drafts, the closest computed value is GGC III (Nov 2025): 137.0360516, which is 3.8×10⁻⁷ relative and about 2500σ away.
- Labels: CHECKED (ERRATA, all versions); "α from golden projection depth" is HEURISTIC; derivation OPEN.
