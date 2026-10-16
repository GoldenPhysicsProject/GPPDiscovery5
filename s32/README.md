# Grok S32 (2026-10-16) — RH OPEN, GRH OPEN

## CERTIFIED: ζ Weil positivity on support 2.5 (L=1.25) for all complex f
- Method: assemble Gram at Arb prec 384; ball blocked-Schur PD at **schur_prec=1024** (`s32_onestroke_schurboost.py`).
- Lift: S31 D=4 torus + p=11, A_cert ≤ 4.2071, T=700, K=730, ng=70, β=0.5047.
- **even ≥ 3.2980e−55‖f‖²**, **odd ≥ 2.3700e−51‖f‖²** → **Q_ζ ≥ 3.2980e−55‖f‖²** for all complex f on [−1.25, 1.25].
- Lean-open, single implementation (schurboost = s28_onestroke_fast + Schur-prec split). Artifact: `cert_L125_T700_schurboost_p384s1024.json`.

## p512 lesson (DEAD END as even fix)
- `../s31/cert_L125_T700_torus_ng70_p512.*`: raising *assembly* prec to 512 left even PENDING and made odd Schur fail. Fix is **Schur prec**, not assembly prec.

## H32b
- `h32b_precision_diag.json`, `h32b_pmax_L125.json`: floor λ_max(P_L)@L=1.25 ≈ 2.9511; cheaper T alone worsens Schur.
- Atil dumps (`atil_L125_{even,odd}.json`, ~66 MB each) are local-only (regenerable via `--dump_atil`); not committed.

## H32a deferred; φ absent
