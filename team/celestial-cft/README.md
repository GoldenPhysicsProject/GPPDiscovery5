# team/celestial-cft — Goddess Of Celestial CFT. RH is OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

**Status: PARTIAL (2026-10-22, early push at the lead's request).** The scans are still running. This file will be replaced by the full writeup (`epstein/EPSTEIN_DH_CFT.md`).

Labels used here:
- **PROVED**: with a citation.
- **CERTIFIED**: Arb ball arithmetic via python-flint 0.9.
- **NUMERICAL**: computational evidence only.
- **CONJECTURE**.

## Task 1: a CFT version of Davenport–Heilbronn, using c=2 Narain T² at σ = τ_Q
### Scripts (`epstein/`)
- `theta.py`: rigorous Λ_Q(s) = (√|D|/2π)^s Γ(s) Z_Q(s), computed with the incomplete-gamma theta split and a rigorous tail bound. It is valid for point s. Ball inputs blow up at large t because of the Γ cancellation.
- `zq_L.py`: rigorous Z_Q(s) by genus-character decomposition, for D = −4, −20 and −84 (2-torsion class groups). It agrees with `theta.py` to 27 digits at s = 3 and at 0.8 + 37.2i (NUMERICAL agreement of two rigorous enclosures).
- `certify_L.py`: CERTIFIED argument principle. Segment balls are checked to exclude 0 and to lie in a half-plane; otherwise the segment is bisected.
- `scan.py`: NUMERICAL float-level winding scan plus Newton, for all forms (including D = −23 and −56).
- `narain.py`, `momentum_identity.py`, `scalar_zeta_identity.py`: CFT checks.

### First results
- **CERTIFIED:** Z_{x²+5y²}(s) has exactly one zero in the box centred at **0.9079933666439 + 89.5598279786i**, half-width 1e−8. So there is an off-line zero with Re s ≈ 0.908.
- **CERTIFIED:** [0.52, 3] × [81, 91] contains exactly one zero of Z_{x²+5y²}.
- **NUMERICAL** (narain_checks.json, M=5 lattice truncation):
  - the Narain T² CFT at σ = τ_Q, ρ = i, for all 12 forms, has Z(τ) = Z(−1/τ) to 1e−9 to 1e−30 relative, limited by truncation;
  - Z(τ+1) = Z(τ) to 1e−31;
  - the vacuum is unique and the gap is > 0.
- **PROVED:** unitarity holds. Z = |η|^{−4} Σ_lattice q^{p_L²/2} q̄^{p_R²/2}, so the multiplicities are lattice counts × coefficients of 1/η², all ≥ 0. Standard Narain construction: Narain 1986; Polchinski vol. 1 §8.4.
- **PROVED** (BCFMP arXiv:2107.10744 eq. 3.25, plus the coprime/full Eisenstein relation): the scalar-primary spectral zeta satisfies S(s) = Σ'_{scalars} Δ^{−s} = |D|^{s/2} E_s(ρ) Z_Q(s).
  - **NUMERICAL** check at s = 3, agreement to about 1e−7: `scalar_zeta_identity.log`.
  - So the off-line zeros of Z_Q are zeros of a genuine unitary, modular-invariant CFT's scalar spectral zeta.
