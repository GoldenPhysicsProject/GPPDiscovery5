# Renormalized pole−prime pairing: ζ vs F_θ

**Status: NUMERICAL (float). RH OPEN. Not a Lean claim.**

## Definition (even compact g)

```
Q_ζ(g)  = ∫_R g(u)(e^{u/2}+e^{-u/2}) du
        − 2 Σ_n Λ(n) n^{-1/2} g(log n)

Q_θ(g)  = ∫_R g(u)(e^{u/2}+e^{-u/2})(e^{θu}+e^{-θu}) du
        − 2 Σ_n Λ(n) n^{-1/2}(n^θ+n^{-θ}) g(log n)
```

Equivalent form for ζ (M_prime = 2Σ Λ(n)n^{-1/2} δ_{log n}):

```
Q_ζ = ∫_0^∞ g(u)[2e^{u/2} du − dM_prime] + 2∫_0^∞ g(u) e^{-u/2} du
```

Archimedean / digamma term is **omitted**. This is not the full Weil form.

## Ledger check (before claiming novelty)

- GPPVerify: no Lean file for this renormalized F_θ pairing.
  `ShiftedFisherKernelNoGo` (#346): M≥0 does **not** kill F_θ (Λ≥0 type).
  `HagedornGapShifted` (#347): gas divergence at 1+θ (different).
  `ShiftedProductControl`: F_θ has FE + Euler product, off-line zeros.
- CLAUDE_CODE_COVERAGE / notes: F_θ full-Weil window negativity INTERVAL-CERTIFIED
  elsewhere; **renormalized (pole−prime only) F_θ test was explicitly TODO**
  (Grok levy note + Codex 2026-10-14_prime_jump_nonlevy ask).
- Compensated Λ(n)(1−1/n) form is a different object and is indefinite even for ζ
  (Claude notes 2026-10).

## Results (float)

1. **Centered Yoshida / Gaussian / autocorr bumps (h≤12):** Q_ζ and Q_θ both stay >0.
2. **Paired positive bumps** g=Yoshida(·−c)+Yoshida(·+c) (positive as functions,
   **not** of positive type): 6 high-res witnesses with Q_ζ≥0 and Q_{1/4}<0
   (e.g. c=9.0, w=0.4: Q_ζ≈+1.6e−3, Q_{1/4}≈−0.167). See `result_verified.json`.
3. **Positive-type class** (autocorrelation of cos-modulated bumps): **no separation**.
   Whenever Q_θ<0, Q_ζ<0 too (`result_postype.json`). Translation of a bump does not
   change its autocorrelation (Fourier phase), so positive-type tests cannot sit on a
   single Chebyshev-error spike the way paired positive bumps can.
4. **Q_ζ itself goes negative** on some paired bumps and some modulated autocorrelations.
   So “renorm ≥0” is not a ζ criterion.

## Verdict

**DEAD END as a Weil separator:** on the positive-type class (the one Weil uses),
the renormalized pole−prime form does not separate ζ from F_θ in the scanned range.
Off-class positive bumps show a float separation, but that is not Weil-admissible.

Interpretation: the off-class separation is expected, not structural. F_θ's renormalized
remainder is the Chebyshev-error measure reweighted by ~2cosh(θu), so any local dip of
ψ(x)−x is amplified by e^{θu}; a positive bump on that dip goes negative for F_θ first.
The full Weil form (with archimedean term) is where F_θ is known to fail (Claude Code's
interval-certified windows h=17/25, 19/25). The arch comparison in `renorm_pairing_v2.py`
uses an unvalidated normalization and should NOT be cited.

No exact math. No hand-off. No Arb run (float already conclusive for the dead-end).

## Run

```
/workspace/venv/bin/python renorm_pairing.py
/workspace/venv/bin/python renorm_pairing_v2.py
# plus the inline fine/postype scripts that wrote result_*.json
```
