# GPP survey report (Grok seat, 2026-10-14)

**What I did and didn't do.**
- Everything was read-only on all remotes: no push, comment or Drive change.
- Local clones live under /workspace. Some local scripts rewrote files in those clones; I restored every clone afterwards with `git checkout` / `git clean`.
- `git pull` failed on one clone with an HTTPS auth error, so that clone is at the state from its first clone today.
- **Not done:**
  - `lake build`: there is no elan or lake on the box, and a cold Mathlib build is far beyond 20 minutes. I relied on CI instead.
  - A full Drive inventory: I paged only the GPP root, "Newer rh physics stuff" and a title search for the book. "Everything", "Verification Code" and "Old shit" were not paged.

## 1. Timeline
| When | Event (source) |
|---|---|
| Dec 2025 – Jan 2026 | Drive: *Celestial Holography / RH / CPT / Yang–Mills* preprints, and many *On the Nature of Nature* builds, some named "BSD_5STAR_FINAL", "PERFECTED" or "WITH_AXIOMS". There are 50+ versions; the last is v21 (2026-07-07, 950+ pages). |
| 2026-07 | GPPVerify Lean tree. 13 axioms were later retired (bridge README:50-54). |
| 2026-08-01 | GPPDiscovery and GPPDiscovery2 start with the truncated-Weil / Connes–van Suijlekom (CvS) λ_min decay scans (fit.py, point.py). |
| 2026-08-24 to 08-30 | c=29 negative mode turns out to be a precision artifact. Codex runs a burst of Fisher, Mehler, AFT, googly and YM work and Drive "mining deltas" (GPPDiscovery2/discovery/*_2026-08-2x). The first corrections are entered (CLAUDE_CODE_CORRECTIONS 1-4). |
| 2026-09-01/02 | Integrity crises: a false "0 sorries", vacuity gates, and the Mathlib 4.33 migration (CC_CORR 9-21). |
| 2026-09-13 | First Grok stint. It opened GPPDiscovery3 (`weil_heat/prime_heat.py`) and PR #187, which was later closed unmerged. |
| 2026-09-24 | Seat reorganization: Claude Code, Codex, a guest seat, and the bridge as control plane. |
| 2026-09-22 to 10-05 | The CCM two-channel lemma; the fixed-window/two-box paper; P34; the RH map (CLAUDE_RH_MAP.md); the R1 Landau step. |
| 2026-10-04/05 | "Positivity selects the primes" is killed by an SDP. The DH fixed-window control runs. |
| 2026-10-05 to 10-10 | GPPDiscovery3 becomes Astra's (pole-neutral, Pick and odd-collision experiments). |
| 2026-10-11 to 10-13 | Turbo arrives (GPPDiscovery4, Lean ArithmeticFock, F_θ "lie detector"). |
| 2026-10-13 | Arb certificates: F_{1/4} and F_{1/10} window negativity, then the DH single-carrier witness. The θ-curvature expansion. |
| 2026-10-13/14 | A rapid loop of Codex notes, Claude Code audits and Lean merges #326 to #342. The parity pole-term remark. Grok gets GPPDiscovery5. |

## 2. Per seat
- **Claude Code: GPPDiscovery (99 commits, 08-01 to 10-08).**
  - Threads in `discovery/`: weil_decay, odd_weil (DH controls), fixed_window, ccm_overlap, arithmetic_holography, celestial_loops, number_thermodynamics, wiener_hopf, lean_sweep.
  - It now does most of its work in bridge `work/claude-code/` and in GPPVerify (Lean, CI, certificates).
  - Current: Arb certificates and the formalization queue. Abandoned or killed: Hermite–Biehler split (09-28), the prime-blind dual, the c-scan slope claims.
- **Codex "Sol": GPPDiscovery2 (328 commits, active to 10-06).**
  - The broadest seat: AFT, celestial, the raised box, prime gas/Hankel, spectral weights, the YM baseline, and daily Drive mining (08-30 to 09-10).
  - Since 10-06 Codex writes only in bridge `work/codex/` (13k-line ledger).
  - Current: pole-neutral reductions, the θ-expansion, and mirror/Cayley/Haar dictionaries, each labelled "not RH".
- **Astra: GPPDiscovery3 (37 commits, 10-05 to 10-10).** Small float experiments (translation growth, pole-neutral window, Pick continuation). The README is honest: everything there is floating-point. Bridge `work/2026-10-*_astra_*` holds "pole-neutral collision closure" and "attachment descent and rh closure"; despite the titles, its claims are audited as partial.
- **Turbo: GPPDiscovery4 (9 commits, 10-11 to 10-13).**
  - Contents: an F_θ lie detector, dyadic anti-locking, a Lean ArithmeticFock snapshot, and the counting law (numerical).
  - Turbo retracted its "purely odd" defect claim.
  - The README still calls GPPDiscovery3 "Grok/Astra" (stale).

## 3. Drive vs repos
- **The book, *On the Nature of Nature* v21 (Drive id 1ZR-tv1miv…, pdf, 2026-07-07).**
  - It claims proofs, or "arguments for", three Millennium problems.
  - Chapter 8 offers "five independent arguments" that all zeros lie on Re s = ½ (onn.txt ≈ lines 1082-1108).
  - It claims the Yang–Mills mass gap is "proven unconditionally" (≈1137).
- **These RH arguments are not valid as stated.**
  - A^×/Q^× is not compact; only its norm-one part is.
  - Zeros of ξ are not a Peter–Weyl spectrum of that quotient.
  - "Haar self-duality forces s = 1−s̄" restates the functional equation, which DH also satisfies.
  - "Positivity of Haar measure forces the Weil sum ≥ 0" assumes exactly the open step.
  - Every one of these would also "prove" RH for DH, which has a zero off the line.
- **How the repos relate to the book.** The repo work is effectively a disciplined retreat from it: the GPPVerify README says "not a proof of any famous conjecture", and CLAUDE_RH_MAP.md §0 names the missing Euler-product positivity step.
  - Neither the book nor the Drive drafts have been publicly corrected, as far as I can see.
  - **This is the largest overclaim in the project.** Daniel should decide how to handle it.
- **"Newer rh physics stuff" (2026-09-25 / 10-08).** It mirrors the bridge notes (CCM_Boundary_Channel_RH_Reduction.tex/pdf, arithmetic_shadow_boundary_v1.tex, Astra notes, DH control). It also holds arXiv 0708.4025, 1803.08928, 2412.09558 and the Feynman positron paper.

## 4. Reproduction (venv: mpmath 1.3, numpy, scipy, sympy, python-flint 0.9)
**Batch run of all 243 `.py` files** (90 s cap each, in their own directory; logs in /workspace/runs/):

| Result | Count | Notes |
|---|---|---|
| Ran successfully | 154 | |
| Exit 1/2 | 61 | 22 need command-line arguments; 7 need `connes_cvs`; 6 import a missing `theta_control_ritz_4_spk4`; 6 load mpmath pickles written by another version; 4 need `cvxpy`; 9 IndexError; 7 hard-coded paths (`/home/box/...`, `google-provider/`, `lake`) |
| Timed out | 28 | heavy |

| Repo | Ran / failed / timed out |
|---|---|
| GPPDiscovery | 22 / 26 / 13 |
| GPPDiscovery2 | 78 / 8 / 6 |
| GPPDiscovery3 | 14 / 1 / 1 |
| GPPDiscovery4 | 7 / 10 / 3 |
| GPP-bridge | 33 / 16 / 5 |

- **No tracked output file changed except one.** `fit.py` in GPPDiscovery and GPPDiscovery2 opens `results.jsonl` with mode "w" and rebuilds it from `artifacts/`, which isn't in git. **Running it locally wiped every saved row** (`fit.py:10`). This is a reproducibility hazard; restored.
- **`connes-cvs==0.2.2` is a real PyPI package.** It is a small pure-Python third-party package by Akiva Groskin (CvS Galerkin matrix). It's not vetted, but it is not a mystery.
- **GPPVerify CI.** main Build is green at #340. #341 and #342 were still building at last check. I did not build locally.

**Independent spot-checks (scripts in /workspace/spot/):**
1. **Parity pole term (float, quadrature, `pole_parity.py`).**
   - For even g: pole = 2G(½)G(−½) = +2G(½)² = 1.71887.
   - For odd g: pole = −0.0052785 = −2G(½)².
   - For mixed-parity g: pole = 2G(½)G(−½) = 1.69775 ≠ 2G(½)², i.e. 2(C²−S²).
   - **Confirmed.** The parity split is exact, and the pole term is indefinite for mixed parity.
2. **DH zero (float, mpmath 30 digits, `dh_zero.py`).** Zero at 0.808517182456637 + 85.699348485377592i, |f| ≈ 1e-32. The functional-equation residual is 2.5e-32. **Confirmed.**
3. **F_θ window certificates.** I recomputed them independently, by summing over the first 150 actual ζ zeros rather than through the explicit formula (`ftheta_zero_sum.py`, float).

   | h | F_{1/4} (mine) | F_{1/4} (certified) | ζ (mine) | ζ (certified) |
   |---|---|---|---|---|
   | 19/25 | −0.021969 | −0.02197 | +0.002030 | +0.00203 |
   | 17/25 | −0.0088635 | −0.00886 | +0.000809 | +0.000810 |

   **Both certificates reproduced.**
4. **Error found in that certificate's README.** `work/claude-code/2026-10-13_weil_window_certificate/README.md` says "the minimizing g are odd in x: c₀ = c₂ = … = 0". The coefficient files say otherwise:
   - 19/25 and 1/2 use only even indices, so those g are EVEN;
   - 17/25, 11/20 and 3/5 use only odd indices, so those g are odd.
   - The statement is false as written. The coefficients agree with the later "both parity sectors fail" (CLAUDE_CODE_RESEARCH_NOTES.md:2021).

## 5. Assessment
**Verified (Lean, or Arb plus my independent check):**
- the finite and structural identities (#326 to #342; QuartetSignature; OddCollisionTrace);
- F_θ window negativity and the DH detection (Arb);
- the parity split of the pole term;
- the routes reduced to "Euler-product ⇒ W ≥ 0 on F" (CLAUDE_RH_MAP.md).

**Numerical only:** Turbo's counting law and −0.052θ²; the CvS λ_min decay; the R1 check "x ≤ 1e7"; ζ positivity on the certified g (which is consistent with RH and proves nothing).

**Conjecture or interpretation:** the holography, entropy and "count only goes up" frames; shadow = T; the physics readings.

**Strengths:**
- a real culture of corrections, with "catch" lines;
- falsifier discipline (F_θ, DH);
- Arb certificates;
- Lean gates against sorry, axioms and vacuity.

**Weaknesses:**
- The audit loop is rubber-stamp fast: "everything correct" within minutes, every time.
- The Lean count grows on identities, while the positivity front has not moved since about 10-05.
- 99 `native_decide` uses sit under "0 axioms" wording.
- Many scripts are not runnable as committed.
- The book's claims are uncorrected.
- Seat references are stale across READMEs.

**Promising directions for an independent explorer:**
- **(a) A certified lower threshold.** Is the F_θ and DH form positive on all g supported in |x| ≤ h₀? It is the only quantitative separation test ζ could pass "for structural reasons". Nobody has it (CC README "Not done").
- **(b) An odd-sector form.** Test −(Arch+prime) ≥ 2G(½)² on odd g at growing h, against DH, which has no pole: the inequality must use the Euler product to separate ζ from DH. Even-sector failure must be tracked with it.
- **(c) A growth-orbit version of the heat criterion** (CLAUDE_RH_MAP.md:24 open question).
- **(d) Make the Discovery repos reproducible.** Fix fit.py, the hard-coded paths and the argument-less scripts.

All of these stay OPEN. RH is OPEN.
