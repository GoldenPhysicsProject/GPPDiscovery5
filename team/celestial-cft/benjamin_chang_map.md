# Benjamin–Chang (arXiv:2208.02259) mapped onto GPP objects. RH is OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Goddess Of Celestial CFT, 2026-10-22.

## Labels
- **PROVED** means it is in a cited paper, or it is a two-line consequence of cited equations.
- **CERTIFIED** means Arb.
- **NUMERICAL** means computational evidence only.
- **CONJECTURE** means a heuristic or an open question.

Equation numbers below are from BC (arXiv:2208.02259, the version served by arxiv.org/abs on 2026-10-22) and from BCFMP (arXiv:2107.10744).

## 0. Dictionary

| BC object | GPP object |
|---|---|
| scalar U(1)^c primaries 𝒮, constrained Epstein zeta 𝓔^c_s = Σ_{Δ∈𝒮}(2Δ)^{−s} (3.4) | the "spectral side" (dynamical, theory-dependent) |
| Λ(s) = π^{−s}Γ(s)ζ(2s), Λ(s) = Λ(½ − s) (2.8–2.9); E_s's constant term Λ(1−s)/Λ(s) (2.7) | ξ, entering kinematically through the SL(2,Z) scattering phase |
| poles at s = (1+z_k)/2 (App. A) | the zeros ρ = z_k of ζ (Re = 3/4 ⇔ Re ρ = ½) |
| tower coefficients δ_{k,c}(μ) (3.13); for c=2, (C.21): δ_k ∝ π^{z_k/2} Λ((1+z_k)/2)² E_{(1+z_k)/2}(ρ) E_{(1+z_k)/2}(σ) / (Γ(z_k/2) ζ′(z_k)) | the zero-side weights of an explicit formula |
| b(n) = Σ_{k\|n} kμ(k), the coefficients of ζ(2s)/ζ(2s−1) (3.17–3.18) | the arithmetic (Möbius / Euler-product) input |
| functionals F^φ[t^s] = ζ(s)M_φ(s) (3.29–3.32), from Bondarenko–Radchenko–Seip 2005.02996 | test functions whose Mellin transform carries a factor ζ(s), i.e. vanishes at every zero |
| RH ⇔ envelope T^{c/2−3/4} of the oscillating term in (5.1) | ModularResonanceDepth.lean: depth ¼ |

## 1. Which crossing statement corresponds to window Weil positivity?
- **PROVED (structure).** BC's scalar crossing equation (3.22) has the same form as a Guinand–Weil explicit formula. It reads: spectral sum over 𝒮 = Cardy terms + Σ_k Re(δ_k y^{…−z_k/2}) + dual spectral sum with kernel b(n)·U(…). The zero tower comes from the poles of 1/Λ(1−s) (App. A).
- **PROVED (from BC §3.3).** The positivity BC actually use is the opposite of Weil positivity:
  - BC remove the zero tower, either with F^φ (Mellin ∝ ζ(s), which kills every t^{z_k}) or with the k ↔ 1/k antisymmetric form (B.12).
  - Their positivity is then unitarity, d(Δ) ≥ 0, on the spectral side.
  - Weil positivity instead asks that the zero side be a positive quadratic form: W(f∗f̃) = Σ_ρ |f̂(γ_ρ)|² ≥ 0.
  - BC state explicitly that the δ_k are "not sign-definite for any k" (§3.3).
- **Verdict (PROVED as a reading of BC).** No crossing statement in BC is equivalent to window Weil positivity. The BC/BRS functional family is the "zero-blind" dual of our window problem (cft_hooks H1 item 3). It is not a positivity certificate for zeros.
- **Nearest analogue (PROVED algebra, CONJECTURE as a route).** At c=2 with ρ = σ, the numerator of δ_k is a perfect square, [Λ((1+z)/2) E_{(1+z)/2}(σ)]². That parallels the |f̂(γ)|² of Weil, but the weight 1/(Γ(z_k/2) ζ′(z_k)) is sign-wild. The analogue of window Weil positivity would therefore be:
  - "Re Σ_k w_k |a_k|² ≥ 0 for all moduli σ",
  - where a_k = Λ((1+z_k)/2) E_{(1+z_k)/2}(σ) and w_k = π^{z_k/2}/(Γ(z_k/2) ζ′(z_k)).
  - I have not checked whether this follows from, or implies, RH. It is open, and it is a candidate H-item.

## 2. Does the spectral-density framing distinguish ζ from DH and from the Narain fakes?
- **PROVED (BC App. A + BCFMP 3.25).** ζ enters BC **kinematically**. The tower locations (1+z_k)/2 are those of ζ for every U(1)^c CFT, fixed by SL(2,Z) through Λ(1−s)/Λ(s). The theory enters only through δ_k(μ).
- **PROVED (BCFMP eq. 3.25 + 2ζ(2s)E_s(τ_Q) = (√|D|/2)^s Z_Q(s)).** The Narain fake enters **dynamically**.
  - At σ = τ_Q, B = 0, ρ = iρ₂, the scalar spectral zeta is S(s) = Σ_{𝒮}Δ^{−s} = |D|^{s/2} E_s(ρ) Z_Q(s).
  - NUMERICAL check at s=3, to about 1e−7: `epstein/scalar_zeta_identity.log`.
  - Its zeros include the zeros of Z_Q. For D=−20 and −84 these are **CERTIFIED off the line** (`epstein/certified_zero_boxes.json`: 81 boxes for t ≤ 100, up to Re s ≈ 2.08).
- **PROVED (consequence).** Zeros of 𝓔^c_s are invisible to the BC high-T asymptotics (5.1), which see only poles. So:
  - a genuine unitary, modular-invariant CFT can have a spectral zeta with off-line zeros, and BC's RH criterion is unaffected;
  - the Narain fakes are **not** counterexamples to the BC reformulation;
  - they *are* counterexamples to any "spectral-zeta RH from unitarity + modular invariance" heuristic.
- **CONJECTURE (heuristic).**
  - DH itself cannot appear kinematically at any level. For Γ₀(N) the Eisenstein scattering matrix is built from Dirichlet L-ratios (Huxley 1984, from memory, unverified), never from non-Euler combinations like DH.
  - The Epstein fakes are the CFT avatars of DH: Z_{x²+5y²} = ζ·L(χ_{−20}) + L(χ_{−4})·L(χ_5), a sum of two Euler products, with χ_5 of conductor 5 just as in DH.
- **Consequence for the battery (PROVED structure, NUMERICAL data).** BC's framing distinguishes ζ from DH and the Narain fakes by *where* they enter (kinematic poles vs dynamical zeros), **not by any positivity**. So it does not supply the missing "Euler product ⇒ W ≥ 0" input.
- **Footnote-14 subtlety (PROVED algebra).** For σ = τ_Q, δ_k ∝ E_{(1+z_k)/2}(τ_Q) ∝ Z_Q((1+z_k)/2)/ζ(1+z_k). A ζ-tower is suppressed iff Z_Q vanishes at the exact point (1+z_k)/2, which has Re = ¾ under RH. Our Z_Q zeros near Re ≈ 0.75 make this a measurable coincidence test; CONJECTURE: it never happens.

## 3. Prior art for this note
- Benjamin–Chang 2208.02259;
- BCFMP 2107.10744 (eq. 3.25, the c=2 overlap 2Λ(s)E_s(ρ)E_s(σ));
- Angelantonj–Florakis–Pioline 1110.5318 and Obers–Pioline hep-th/9903113 (constrained Epstein series);
- Bondarenko–Radchenko–Seip 2005.02996 (the functionals).
- I found no CFT paper that points out off-line zeros of Narain scalar spectral zetas. This is a search, not a proof of absence (§ prior art in `epstein/EPSTEIN_DH_CFT.md`).
