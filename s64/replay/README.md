# Support-2.8 replay package (Grok S64, 2026-10-17)

RH OPEN. This is a finite-window certificate only: Q_ζ(f) ≥ 7.4196e−78 ‖f‖₂² for every complex f ∈ L²[−7/5, 7/5] (support 2.8 = log 16.44). It is an external Arb + exact-integer computation. It is **not Lean** and has **not been independently replayed**, and it is not to be called a record until it is.

Convention: S14/Zhu normalisation, which is the same as Mori VIII (1.1) (checked: we reproduce his Schur constant R_w = 4.0915 exactly, see `../lift_cmp`). f ∈ L²[−L, L], L = 7/5, L² norm.

## What is being certified (S35 conditional theorem, one parity at a time)
- **H1 (lift).** On frequencies t ≥ T, the lifted prime symbol satisfies S_lift(t) ≤ A_cert. Here S_lift = Σ_{n<e^{2L}} 2Λ(n)n^{−1/2}cos(t log n) + h(t), and h is a sum of "invisible" characters cos(k·log p·t) with |k·log p| ≥ 2L + ε. Each such character integrates to 0 against |f̂|² for f ∈ L²[−L, L].
  - even: `torus_L1.40_buf_D9D8_cert.json`, A_cert = 5.189 (S58).
  - odd: `torus_L1.40_buf_D10D8_cert.json`, A_cert = 5.0536 (S61).
  - Each is certified by two B&B implementations (float `s34_bb_D4` and libm-free interval `s35_interval_bb`).
- **H2 (assembly).** Ã is the 384-bit midpoint Gram matrix (one-stroke Zhu-type, spherical-Bessel basis, `s53_onestroke_par.py`) of Q restricted to frequencies below T, minus the tail. All quadrature errors are collected into rE.
- **H3 (PD).** λ_min(Ã) ≥ λ₀ > 0. Then Q ≥ min(λ₀ − rE − δ_str, lamD − ε_B)·‖f‖².

| parity | T | K | β | Ã sha256 | λ_min(Ã) ≥ | final |
|---|---|---|---|---|---|---|
| even (S61) | 1310 | 1497 | 0.150142 | 94238400e0b9…3e1a69 | 7.4197e−78 | 7.4196e−78 |
| odd (S63) | 1277 | 1459 | 0.260009 | db6365a4f8d0…524071 | 1.24599e−73 | 1.2459e−73 |

The full sha256s are in `MANIFEST.sha256`. Certificate JSONs: `cert_L140_EVEN_S61_FINAL.json`, `cert_L140_ODD_S63_FINAL.json`, `cert_L140_SUPPORT28_S63.json` (summary, including simplicity λ₂(even) ≥ 8.6759e−70).

## Level 1: PD replay from the stored dumps (about 25 min per parity on 1 core, peak about 3 GB)
The dumps are on the Grok box at `/workspace/GPPDiscovery5/s58/run140/atil_L140_T1310_K1497_even.json` (279 MB) and `/workspace/GPPDiscovery5/s58/run140/atil_L140_T1277_K1459_odd.json` (265 MB). They are not in git: they are too big. Get them from Daniel's box, or regenerate them with Level 2.

```
./replay_pd.sh even /path/atil_L140_T1310_K1497_even.json
./replay_pd.sh odd  /path/atil_L140_T1277_K1459_odd.json
```
Writes `replay_out_<parity>/REPLAY_RESULT.json`. PASS requires all of the following:
- the sha256 matches;
- `s36_cholres` (Arb, 640 bits) and `s37_ldlt_exact` (pure integers, no Arb) both give λ_min > 0, agreeing with the certificate to 4 significant digits;
- the negative control at 1.02·λ_est FAILS;
- the s55 inertia shows 0 negative and 0 zero pivots.

**The best replay is a third, independent PD implementation written by the replayer.** The math to check is just "Ã − λ₀I is positive definite" for the decimal matrix in the dump, with each entry's error at most 1 unit in its last printed digit. Any rigorous method works (exact rational LDLᵀ, interval Cholesky, etc.).

## Level 2: re-assemble Ã (about 2 h per parity with 3 workers, peak about 9 GB)
From GPPDiscovery5 at the commit of this package:
```
cd s58/run140
# even (S61): S58 lift
SPHJ_BOOST=2.0 T=1310 K=1497 PARS=0 W=3 TORUS=torus_L1.40_buf_D9D8_cert.json ./run.sh
# odd (S63): S61 lift
SPHJ_BOOST=2.0 T=1277 K=1459 PARS=1 W=3 TORUS=torus_L1.40_buf_D10D8_cert.json ./run.sh
```
Then compare the sha256 (expected to be deterministic at prec 384; this has NOT been tested) and run Level 1. If the sha differs, compare entrywise (expect agreement to at least 1e−100). Check that `rE=inf` does not appear in the log. SPHJ_BOOST ≥ 2.0 is required: 1.75 gave a void run in S60.

## Level 3: re-certify the lifts (about 10 min to a few hours)
- Even lift: `s58/run140/lift/lift.sh`, or just the two B&B checks on `bb_p235/torus_L1.4_D9_cert.json` and `bb_p71113/torus_L1.4_D8_cert.json` with `s34_bb_D4.py` and `s35_interval_bb.py`, then `s35_check_merge.py`.
- Odd lift: `s61/lift/cert_r3/cert.sh`.
- Independent check (easy to write): every character in the torus JSON has |k·log p| ≥ 2L + ε. For the sup bound, sample sup_t S_lift(t) on a fine grid; it must stay below A_cert. A sample is necessary, not sufficient.

## Known caveats
- H2 relies on Zhu's Theorem 1.1 framework (arXiv:2608.24827) for the one-stroke tail inequality, as validated in S27.
- The lift is our own method. S64 found it is a *frequency-side* bound on the same quantity λ_max(C_p|L²[−L,L]) that Mori VIII and Liu B.4 bound on the position side with a weighted Schur test. At L = 1.4 it is numerically *weaker* (see `../lift_cmp`).

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
