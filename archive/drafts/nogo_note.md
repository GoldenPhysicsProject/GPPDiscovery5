# Grok: attack on the no-go ledger and the classical RH equivalences (2026-10-17)

Seat: Grok (Loverbot). Asked by Daniel (Oct 17): "Run through our no-go list in Verify and see if you can see a way to defeat any of them, and also look up other equivalencies of other people in the past and see if you can do things they couldn't."

**RH is OPEN. GRH is OPEN.** Nothing here proves or implies RH. Numerics are evidence, not proof.

Labels: **PROVED** (complete argument, Lean or paper), **CERTIFIED** (Arb interval arithmetic, not Lean), **NUMERICAL**, **HEURISTIC**, **CONJECTURE**, **CLASSICAL** (known in the literature, not ours).

Sources read:
- GPPVerify `b2292fb`: `GppVerify/NoGo/*` (22 files plus `Ledger.lean`), all `RiemannHypothesis/*NoGo*.lean` (11 files), `StubRegistry.lean`, `RHCriteriaAtlas.lean`, `LiCriterion.lean`.
- GPP-bridge `b5f080e`: `CLAUDE_RH_MAP.md`, `CLAUDE_CORRECTIONS.md`, `GROK_RESEARCH_NOTES.md` (S1–S60), and grep over the `CODEX_*` and `CLAUDE_*` ledgers.
- Web: the AIM RH-equivalents list, Zhu arXiv:2608.24827 v2, Liu (Sept 14 2026), Mori's papers VI–VIII and his Lean log 5 proof (Zenodo), dBN 2026 bounds, and GORZ (PNAS 2019).
- Scripts: `GPPDiscovery5/nogo_attack/`. I did not touch the running hourly jobs (PID 83815 family and the s61 job).

## 0. A correction to our record claims (most important item)

The ledger check found published work that **nobody on the team had logged**: a grep for Liu, Mori and 17/16 over GROK_*, CLAUDE_*, CODEX_*, CONVERSATION.md and work/ came back empty.

| Who | Date | Window (full support = 2L) | Test functions | Status |
|---|---|---|---|---|
| Zhu, arXiv:2608.24827 v2 | Aug/Sept 2026 | 1.6 | all complex | Arb certificate; his support-2.38 certificate was withdrawn in v2 |
| V. Liu, "Certified Weil positivity beyond the unit window" | Sept 14 2026 | 2.0 (Q ≥ 2^−151) and 2·17/16 = 2.125 (Q ≥ 2^−49162) | all complex | exact integer Gram certificates, not Lean |
| R. Mori, paper VI | Oct 8 2026 | log 10 = 2.303 | all complex C² | Arb, external |
| R. Mori, paper VII | Oct 9 2026 | log 12 = 2.485 | all complex C² | Arb, external |
| R. Mori, paper VIII (zenodo 23281235) | **Oct 10 2026** | **log 15 = 2.70805**, Q ≥ 2^−2904‖f‖² | all complex C² | Arb; finite margins 2.40e−70 (even) and 3.95e−66 (odd); 832×832 matrix per parity |
| R. Mori, Lean 4 (zenodo 23034002) | 2026 | log 5 = 1.609 | all complex C² | **fully Lean-kernel-checked**, from the explicit formula to the last inequality |
| **Us** (S27–S57) | Oct 15–17 | full complex 2.0, 2.38, 2.5, 2.6, 2.7, **2.701**; **even only 2.76** | | Arb, one or two implementations |

Consequences:
1. **Our messages to Daniel that support 2.0, 2.38, 2.5, 2.6 and 2.7 were "bigger than anyone has published" were WRONG.** Liu already had 2.0 and 2.125, and Mori had 2.303, 2.485 and 2.708 before our posts. Our best full-complex window, 2.701, is just below Mori's 2.708. This goes into GROK_CORRECTIONS.
2. The parts that are still ours:
   - the **even-sector** certificate at 2.76 > log 15;
   - the support-2.8 run now in progress. If both parities certify, 2.8 would be the first full-complex window past log 15 *and* past log 16 = 2.7726, where 16 = 2⁴ enters the prime comb. This is CONDITIONAL; it is not certified yet.
   - the invisible-harmonic lift idea itself. Its novelty against Liu's weighted Schur row-sum bound (his App. B.4) and Mori's "new weight in the Schur bound for the prime operator" is **OPEN**. I have not read Mori's paper VIII body yet.
3. Cross-check (NUMERICAL consistency, not a proof): Mori's finite even margin at L = 1.35403 is 2.4e−70. Our certified even lower bounds at L = 1.35–1.3505 are 5.3e−70 and 1.2e−69. Both are valid lower bounds at slightly smaller L, so the two independent pipelines agree on the order of magnitude.
4. **Other literature the ledgers lack:**
   - Liu's **Obstruction Theorem**: the "tail-dropping localization" lower comparison has a negative high-frequency limit once L > ½log 8, and no fixed compact correction repairs it. This is a published no-go; I add it to §1.
   - dBN upper bounds in 2026: Platt–Trudgian 0.2 (published); Mosaic 0.1965, Gomila 0.1787854 and Gordon 0.157657 (unreviewed, Gordon with a conditional Lean companion).
   - Alpöge–Furman arXiv:2608.13637 ("more than two thirds of zeros simple and on the line"), cited by Liu; I have not checked it myself.

## 1. The no-go ledger: hypotheses, and whether each can be dropped, weakened or sidestepped

Verdict key:
- **WINDOW-DEFEATED**: our window certificates violate the no-go's conclusion on every finite window, but not uniformly in L.
- **SIDESTEP-KNOWN**: the hint names a known structure that avoids the hypothesis.
- **HARD**: the only way round it is RH-strength.
- **INFO**: a detector or lesson, not an obstruction.
- **N/A**: not an RH obstruction.

| # | No-go (Lean name) | Exact hypotheses | Weaken/drop? | Verdict |
|---|---|---|---|---|
| 1 | FreeJacobiNoEigenvalue | constant-coefficient (free) Jacobi operator | Allow perturbed coefficients. But any real zero sequence is the spectrum of *some* Jacobi matrix, so dropping "free" is vacuous unless the coefficients come from arithmetic. The natural candidate is the Grommer/Hankel Jacobi matrix of ξ's moments; all its b_n > 0 ⇔ RH (S51, CLASSICAL). | HARD |
| 2 | GramDiagonalTautology | only the diagonal window energy is used | The hint asks for a *signed upper* off-diagonal bound. That is exactly what Zhu's comb envelope and our invisible-harmonic lift provide: an upper bound S(t) ≤ A_L on the prime comb. | **WINDOW-DEFEATED** (all certificates up to 2.701/2.76). The lift ceiling ≳ e^L (S28, with Codex) and Zhu's doubly-exponential frequency threshold block uniformity. |
| 3 | EulerDegreeBlind | uses only the degree defect D_p | add the splitting witnesses or the archimedean completion | SIDESTEP-KNOWN (all window work uses the completion) |
| 4 | EulerPoleTruncation | finite Euler truncation plus a bounded correction, Pick nodes approaching the pole | On compact support the explicit formula sees **exactly** the finitely many n < e^{2L}, so there is no truncation error. The pole and primes are renormalised jointly. | **WINDOW-DEFEATED** (structural: this is why window certificates exist at all) |
| 5 | LaguerreCutoff | cutoff N ≤ e^{n/48} for Laguerre order n | As 4: a window has the exact cutoff e^{2L}. HEURISTIC dictionary: Cayley–Laguerre order n ↔ window L ≈ n/96. It says nothing about sufficiency. | WINDOW-DEFEATED (same mechanism) |
| 6–11 | TwoNodePick, QuartetThreeNode, QuartetSignature, QuartetFourNode, FiniteNoMasking(Multi), FiniteLoewnerNoMasking | **finitely** many on-line pieces | The hypothesis "finite" cannot usefully be dropped: with infinitely many, masking = Weil positivity = RH. | INFO: detection is local, masking needs the infinite sum. Matches S46/S47: detection needs a window that resolves one local zero spacing, L ≈ ½log(γ/2π). |
| 12 | CayleyGap44 | (detector) | — | INFO |
| 13 | GaussianHeatOS | Gaussian subordination of one prime length with a positive prime measure | Use the completed K(t); the Gamma part must restore OS positivity | HARD. S52 (PROVED): finite-order CM is unconditional above τH² ≥ 52.16k + 32.1, so RH-blind. Planted CM depth k* ~ γ²log γ/δ. |
| 14 | TwirlNotCollapse | book correction | — | N/A |
| 15 | ConstantVertexSchur | constant vertex | use a vertex built from ζ's own data | SIDESTEP-KNOWN: our pole Schur complement / DC lift *is* ζ's own vertex |
| 16 | FiniteNormResolvent | finite-norm Ω, accretive A | rigged Ω | SIDESTEP-KNOWN. S60 (PROVED): ‖e_s‖² = ξ(Re s), the obstruction sits at σ = 1, and the resulting kernel positivity is RH-blind |
| 17 | ConfluentStieltjesGram | (detector) | — | INFO |
| 18 | HaarSquareNotReflection | Haar positive type only | needs *reflected* (OS) positivity | HARD |
| 19 | **CompactCommutingMixing** | the operator **commutes exactly** with a weakly mixing scale flow | Break scale invariance in a controlled way: the **Connes–Moscovici prolate operator** (PNAS 2022) and the CCMvS semilocal program do exactly that. | **SIDESTEP-KNOWN, and it points to Opening C** (§4) |
| 20 | HermiteFourierIndefinite | Fourier unitarity on the K-trivial sector | positivity from a half-space (OS) structure | HARD |
| 21 | **ShiftedProductControl (F_θ)** | FE + degree-2 Euler product with **non-unitary** Satake p^{±θ} | Add local unitarity / Ramanujan. HagedornGapShifted (Lean) and SatakeTwoPointNoGo (Lean, our S9/S10) both exclude F_θ. Our support-2.6 certificate also excludes it (S45 CERTIFIED). | **DEFEATED by an extra hypothesis we already have in Lean.** The cost: DH is not excluded by unitarity, see 29. |
| 22 | Sharp-window energy (Littlewood) | sharp cutoff | smooth cutoff | SIDESTEP-KNOWN |
| 23 | DH (classical) | FE + growth, **no Euler product** | must use the Euler product | HARD. §2 has the new numbers: window positivity does **not** see DH up to support 2.7 or beyond. |
| 24 | One-sided Cayley–Laguerre (Pringsheim) | — | the one-sided bound is Landau-type ⇔ RH | HARD |
| 25 | First-prime-trace | — | restatement | HARD |
| 26 | BostConnesFiniteNoGo | **static** range projections, finite N | use the dynamical (KMS / time-evolution) channels | OPEN, but that is Connes' adelic program, which is RH-strength |
| 27 | FiniteModelSpaceReflectionNoGo | quotient-level reflection only | annihilate the boundary defect. That is the Nyman–Burnol distance d_N → 0 (§3, row NB) | HARD |
| 28 | FisherZeroLogConcavityNoGo | generic positive density, any log-concavity | Restrict to Pólya-type kernels. Rodgers–Tao (Λ ≥ 0) says ξ has **no slack**, so no "strong positivity with margin" can work | HARD (a meta no-go, see R15) |
| 29 | PositiveSquareOddSectorNoGo | grading + D² ≥ 0 only | needs an arithmetic no-ghost theorem | HARD. NUMERICAL: in every ζ window certificate the **odd** sector is the slack one; even binds, consistent with "even simple ground state" (CCMvS) |
| 30 | PrimeInnerCascadeNoGo | raw product of prime inner channels at a_p = p^{−1/2} | Quotient the common zero at i/2 (the k = 0 mode is the pole). This is the same one-channel pole renormalisation as 15 and 33 | SIDESTEP-PARTIAL: the k ≠ 0 accumulation remains |
| 31 | PrimeWeightNoGo | **diagonal** Hilbert reweighting | non-diagonal completion | WINDOW-DEFEATED in spirit: the lifts are non-diagonal in prime-character space. NUMERICAL S59: cross-group characters buy only ≈1% |
| 32 | ReciprocalPositivity / SatakeTwoPoint / ShiftedFisherKernel | FE + Λ ≥ 0 (+ non-unitary) | add (C): tempered local factors with one pole pair | DEFEATED for F_θ by (C). DH still needs (B) |
| 33 | **ShiftedPoleRankOneNoGo** | **two** negative polar directions; one positive rank-one correction | Note the converse: with **one** negative direction a rank-one positive correction *can* heal. | **Structural crack, §2:** ζ has polar index 1, so the pole Schur complement works for ζ and provably cannot work for F_θ |
| 34 | Liu's Obstruction Theorem (external, new to us) | tail-dropping localization plus compact repair, L > ½log 8 | Do not localize. Zhu's one-stroke, our lift/Schur and Mori's Galerkin all avoid it | N/A for our pipeline |
| 35 | Zhu frequency threshold plus our lift ceiling (S28 with Codex) | **one-stroke** envelope certificates | Non-one-stroke (direct Galerkin plus complement, as Liu and Mori do). Under RH, λ_min(L) ≤ exp(−L e^L) (Zhu), so *every* finite-window certificate loses its margin super-exponentially | **HARD wall for any finite-window numerics as an RH route** |

**Bottom line of §1 (HEURISTIC synthesis, but every row cites a Lean or certified item).** No no-go falls in a way that touches RH. Every "defeat" is either window-wise (rows 2, 4, 5, 31) or an extra hypothesis (rows 15, 19, 21, 32, 33). Read together, the ledger says a proof must use all five of these at once:
- (i) exactly one polar channel (rows 21, 32, 33);
- (ii) the Euler product with Λ(n) ≥ 0 on prime powers (rows 23, 32);
- (iii) local unitarity (rows 21, 32);
- (iv) infinitely many primes and zeros (rows 4–11);
- (v) no slack (rows 28, 35; Rodgers–Tao).

The current falsifiers each violate a different item: F_θ violates (i) and (iii), DH violates (ii), H_{−t} violates (v), and the Haar witness violates all of them. **That is exactly the Selberg-class-with-Ramanujan setting, where nobody knows a counterexample.** It is a map, not a crack.

## 2. New computation: the "polar index" and whether window positivity sees DH

**Polar index (PROVED, elementary).** The pole term of the Weil form is P(f) = 2Re F f(i/2)·conj F f(−i/2). Its Lean/Liu form is 2(|⟨cosh(u/2),f⟩|² − |⟨sinh(u/2),f⟩|²), which has signature (1,1), with the negative direction in the odd sector.
- ζ: one pole pair {0,1}, so **index 1**.
- F_θ: poles ±θ, 1±θ, so **index 2** (ShiftedPoleRankOne, Lean).
- DH: entire (a combination of L(s,χ), L(s,χ̄) mod 5), so **index 0**. There is no polar negative direction at all.

So the ranking is F_θ 2 > ζ 1 > DH 0. Our rank-one pole Schur complement removes exactly ζ's one channel, cannot fix F_θ (Lean), and has nothing to do for DH.

**HEURISTIC consequence.** On small windows DH should be *more* positive than ζ: no polar negative direction, plus a smoother sech archimedean difference. Any separation of ζ from DH must therefore come from the prime side, through b(n) being signed and non-multiplicative.

**NUMERICAL test** (`GPPDiscovery5/nogo_attack/`, S47 exact Arb entries; the Ritz eigen-solve is in mpmath, so this is NUMERICAL):

| M | L (support 2L) | DH even | DH odd | ζ even | ζ odd | source |
|---|---|---|---|---|---|---|
| 88 | 1.3 (2.6) | 9.45e−19 | 8.60e−19 | 2.39e−49 | 9.12e−47 | S47 |
| 128 | 1.35 (2.7) | 3.26e−21 | 2.99e−21 | 5.77e−60 | 2.23e−57 | this note |
| 160 | 1.35 (2.7) | 1.03e−22 | 9.78e−23 | 7.55e−65 | 2.33e−61 | this note |
| 128 | 1.5 (3.0) | 2.78e−21 | 2.73e−21 | 1.63e−68 | 8.70e−67 | this note |

(Ritz λ_min on Legendre·(1−u²)⁴ of degree < M, mp dps 160–200, Arb entry radii ≤ 3e−210. Ritz values are upper bounds for the true sector infimum. The M = 160 runs at L = 1.5 and 1.75 were stopped to free CPU.)

Reading:
- DH is PD at every tested (L, M), with λ_min **20–40 orders of magnitude above ζ's** on the same space. So the generic polynomial window basis does **not** separate ζ from DH through support 3.0 (L = 1.5).
- DH's λ_min keeps falling with M (×~30 per +32 in M at L = 1.35, from 3.3e−21 at M = 128 to 1.0e−22 at M = 160), so its true infimum is OPEN.
- The local-spacing rule (S47) predicts DH's first off-line zero (γ ≈ 85.7) becomes resolvable at L ≈ ½log(85.7/2π) ≈ 1.31. It has **not** produced a negative direction at L = 1.35 (M ≤ 160) or L = 1.5 (M = 128). DH λ_min is almost L-independent while ζ plunges (×1e−8 from L = 1.35 to 1.5), so DH’s floor here is set by basis resolution, not by its off-line zero.
- So: CERTIFIED ζ positivity at support 2.7 (ours) or 2.708 (Mori) is, as far as this basis can tell, a statement that **DH also satisfies** (NUMERICAL). Window positivity at this scale is not yet evidence that the Euler product is being used. That is the falsifier rule (CLAUDE_RH_MAP §4) applied to our own headline.
- Claude Code's carrier-adapted DH detection (half-width 2.5 float, certificates at 3 and 3.5; `2026-10-13_dh_modulated_window`) gives the current upper bound on DH's separation window.

## 3. Historical equivalences: where each is stuck, and whether our tools add anything

| Row | Equivalence (who) | Where the known approach is stuck | Can our tools do something they couldn't? | Verdict |
|---|---|---|---|---|
| W | **Weil positivity** (Weil 1952; Bombieri 2000; Yoshida 1992; Connes–Consani 2021; Zhu, Liu, Mori 2026) | Positivity is known only on finite windows; the cost is doubly exponential (Zhu), and under RH λ_min ≤ exp(−Le^L). | Our certificates are at the frontier: 2.701 full, 2.76 even, 2.8 running. In Lean: log 5 (Mori, fully kernel-checked); GPPVerify has only rung 0 (support < log 2), conditional on arch ≥ 0 (WeilSupportLadder). On paper: log 2 (Yoshida 1992, Connes–Consani 2021). Nothing uniform in L. | Frontier work, not an RH route (row 35) |
| CCM | **Connes–Consani–Moscovici–van Suijlekom** semilocal / prolate program | Needs, per window, a *simple even* ground state (Zhu certified this at 1.6 only); then uniformity. | We can certify evenness and simplicity at larger L. Our data already shows even binding at all L ≤ 1.4 (NUMERICAL). It needs **two-sided** bounds (Opening C). | **Opening C** |
| Li | **Li / Keiper** λ_n ≥ 0 (Li 1997; Bombieri–Lagarias 1999) | λ_n = Weil form on non-compact test functions; positivity for n ≤ N follows from RH verified to height ~√N; asymptotics need RH. | GPPVerify has λ₁ > 0 (LiCriterion.lean). The window machinery does not apply (non-compact support). | DEAD END for us |
| NB | **Nyman–Beurling / Báez-Duarte** density, d_N → 0 | Burnol: d_N² ≥ (Σ m_ρ²/|ρ|²+o(1))/log N, so convergence is at best logarithmic and numerics cannot decide. | Our Gram/Schur tools give certified finite-N values, which say nothing about the limit. FiniteModelSpaceReflectionNoGo = the Nyman–Burnol ghost. | DEAD END |
| Rob | **Robin** σ(n) < e^γ n log log n, n > 5040; **Lagarias** σ(n) ≤ H_n + e^{H_n}log H_n; **Nicolas** (totient) | One-sided arithmetic inequality; known for large classes of n and numerically far out; full ⇔ RH. | Same shape as our S48 result: a one-sided bound suffices by Landau, so G ≥ −CX^ε ⇔ RH. We have no tool that proves one-sided bounds for all x. | DEAD END for us (noted as a Landau-type row) |
| Mer | **Mertens-type** M(x) = O(x^{1/2+ε}); **Redheffer** det R_n = M(n); **Franel–Landau** Farey discrepancy; **Riesz** Σ(−x)^k/((k−1)!ζ(2k)) = O(x^{1/4+ε}); **Hardy–Littlewood** analogue | All are the growth of 1/ζ, equivalent to zero-free Re s > ½; there is no positivity structure. | None. Positivity and Gram tools have nothing to act on. | DEAD END |
| BSY | **Balazard–Saias–Yor**: ∫ log|ζ(½+it)|/(¼+t²) dt = 0; **Volchkov** double-integral identity | Each is a single nonnegative number that equals 0 iff RH; it cannot be certified *equal* to 0. | Arb could bound it ≤ ε, which only excludes off-line zeros with log|ρ/(1−ρ)| > ε/2π, i.e. tiny low ones that are already excluded by Platt–Trudgian. | DEAD END |
| Sp | **Speiser** (ζ′ ≠ 0 in 0 < Re < ½); Levinson–Conrey proportion (now > 2/3 per Alpöge–Furman 2026, cited not checked) | Mollifier length barrier, so 100% is out of reach. | None (different technology). | DEAD END |
| J | **Pólya / Jensen**: all J^{d,n} hyperbolic. **GORZ 2019**: every d for n ≫ d, and all n for d ≤ 8. CNV 1986 (d=2), Dimitrov–Lucas (d=3). | Uniformity in d at small shift n. | S50 (PROVED): a cone + Gauss–Lucas proof of d=2 for all n, plus Laguerre bands. **New estimate (HEURISTIC):** feeding the Platt–Trudgian cone (all u-zeros in \|arg u\| ≤ 1/H, H = 3·10¹²) into the S50 band law gives hyperbolicity only up to order ≈ ¼log(H/2π) ≈ 6–7 < GORZ's 8. | **DEAD END** (logged) |
| dBN | **de Bruijn–Newman** Λ = 0 (Rodgers–Tao 2020: Λ ≥ 0) | Upper bounds improve like 1/log H (0.2 published; 0.158 unreviewed). Λ ≤ 0 is RH. | Our Lean heat layer and S52 CM depth are RH-blind at finite order. **Lesson:** Rodgers–Tao is a *meta-falsifier*. Any mechanism with slack also holds for H_{−t}, which has non-real zeros (S60: kernel positivity holds for H_{−0.2}). | Add H_{−t} to the standard falsifier battery |
| HP | **Hilbert–Pólya / Berry–Keating** xp, Bender–Brody–Müller | No self-adjoint operator with the right domain and spectrum; xp needs boundary conditions that encode ζ. | Our no-go rows 1, 16, 19, 26 and 29 kill every naive candidate. The only survivor consistent with the whole ledger is a scale-breaking (prolate-type) operator, i.e. Opening C. Daniel's Fock/Gibbs "QFT for numbers" (PrimeFockPartition etc.) is the Julia/Spector/Bost–Connes primon gas, CLASSICAL; it lives at Re s > 1 (Hagedorn point s = 1, HagedornGapShifted). | No new crack; Fock layer is a Re s > 1 object |
| Gr | **Grommer / Hankel** (1914), Laguerre–Pólya (Codex 09-28) | Finite N certifies only finite height; N_max ~ log H (S51). | Done (S51 CERTIFIED H6). | DEAD END (logged S51) |
| Su | **Suzuki** screw function −g_ζ ≥ 0 (JLMS 2023); Mittermeier certified to 10^10 | One-sided, Landau-type. | Codex 6215141315 already identified our S(t) with it. | Known, logged |

## 4. Ranked openings (most promising first)

**A. DH separation window: find the first support where certified ζ positivity *fails* for DH.**
- *Crack:* ledger row 23 says a proof must use the Euler product. §2 shows our certificates have not yet demonstrated that they do.
- *Test:* run the actual DH Weil form (S47 exact assembly) on a **carrier-adapted** basis (Legendre + cos/sin(γx) with γ ≈ 85.7, DH's first off-line zero) at L = 1.35–2.5. Pin the smallest L_DH with a CERTIFIED negative direction (one Arb vector suffices). Then push the ζ certificate to that L.
- *Payoff:* the first CERTIFIED theorem of the form "Weil positivity on [−L,L] holds for ζ and fails for DH". That would show window positivity uses the Euler product. It is a statement about *what the proof must contain*, not RH. Realistic only if L_DH ≲ 1.5, and the generic basis already fails at 1.5, so it needs the carrier basis.
- *Owner:* Astra (carrier basis) + Grok (ζ side).

**B. Reconcile with Liu and Mori, and cross-certify.**
- *Crack:* none (this is hygiene), but it decides what is new.
- *Test:* (1) replay Mori's public log-15 certificate data with our independent code, as we did for Zhu in S27; (2) read Mori VIII's Schur weight and Liu's App. B.4 against our invisible-harmonic lift; (3) finish support 2.8 (both parities). It would be the first full-complex window past log 15 and log 16.
- *Owner:* Turbo (independent replay), Grok.

**C. CCMvS spectral hypothesis (simple, even ground state) beyond Zhu's 1.6.**
- *Crack:* CompactCommutingMixing (row 19) kills every scale-invariant boundary operator; the prolate operator is the one known sidestep, and its program needs exactly this per-window hypothesis.
- *Test:*
  - **Evenness:** CERTIFIED lower bound on the odd λ_min above an Arb-evaluated Ritz *upper* bound on the even λ_min. Today: certified odd lower bound ≥ 1.5659e−65 at L = 1.35 (S53). The rigorous L²-normalised even Ritz upper bound at L = 1.35 is RITZ_EVEN (M = 160, `nogo_attack/ritz_l2.py`; validated at L = 0.8 against Zhu's bracket). Evenness is certified iff odd_lower > even_upper: RITZ_VERDICT. (The plain-coefficient-norm values in the §2 table are NOT L²-normalised and are used only for signs.)
  - **Simplicity:** certify Q_even + C|v⟩⟨v| ≥ μ with μ above the even Ritz upper bound (min–max ⇒ λ₂ ≥ μ).
  - Start at L = 0.9–1.1, where margins are ~1e−20 to 1e−30.
- *Payoff:* certified inputs for the CCMvS real-zero approximants on larger windows. Not RH.

**D. Lean pipeline via Mori's kernel-checked log-5 proof.**
- *Crack:* this removes our standing caveat "Arb, not Lean" (and StubRegistry `WeilPositivityAll` stays a Prop).
- *Test:* Codex compares `moriryota/weil-positivity-log5-lean` (license first) with GPPVerify `WeilSupportLadder` / `WeilInvisibleFrequencySupport` / `WeilResidualCoercivity`. Is its explicit-formula-to-finite-inequality chain reusable, and can an ~880×880 integer Gram residual check scale in-kernel?
- *Owner:* Codex.

**E. Falsifier battery = {F_θ, DH, H_{−t}} for every new mechanism.**
- *Crack:* rows 21/23/28 plus Rodgers–Tao.
- *Test:* a standing rule: any proposed positivity mechanism must state which of the three it fails and why. S60 already showed coherent-kernel positivity passes all three, hence RH-blind.
- *Owner:* all seats.

## 5. Dead ends logged (this session)
- "Our support 2.0–2.7 certificates are the largest published": **FALSE** (Liu 2.125; Mori 2.303/2.485/2.708, all earlier).
- Using finite-window Weil certificates as RH evidence or as a DH filter at support ≤ 2.7 (NUMERICAL: DH is PD there).
- The Jensen/GORZ cone route, which reaches order ≈ 6–7 < 8 even with Platt–Trudgian (HEURISTIC estimate).
- Li, Nyman–Beurling/Báez-Duarte, Robin/Lagarias/Nicolas, Mertens/Redheffer/Franel–Landau/Riesz/Hardy–Littlewood, BSY/Volchkov, Speiser/Levinson: no foothold for positivity/Gram/Schur tools.
- Dropping "finite" in the no-masking theorems (that is RH itself).
- The Fock/Gibbs "QFT for numbers" layer as a critical-strip object (it is the primon gas, Re s > 1).

## 6. Hand-offs
- **Codex (Lean):**
  - (D) above;
  - a one-line lemma, the converse of ShiftedPoleRankOne: a 1-dim negative sector *is* healed by a rank-one positive correction with a large enough weight (Schur);
  - log Liu's Obstruction Theorem in `NoGo/Ledger.lean` as an *external, not formalized* row.
- **Astra:** Opening A (DH carrier basis at L = 1.35–2.5; first certified negative vector).
- **Turbo:** Opening B(1), an independent replay of Mori's log-15 data.
- **Sol:** none.

φ: absent. It appears only in DH's κ = √(φ+2) − φ, from the input choice of conductor 5. RH OPEN. GRH OPEN.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
