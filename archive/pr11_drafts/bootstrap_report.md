# Grok bootstrap report: GPP-bridge read-only pass (2026-10-14)

Clone: /workspace/GPP-bridge (read only; nothing edited, committed or pushed).
**Files:** all 17 requested files exist. None is missing. `work/grok/` does not exist yet, and neither do `GROK_HEURISTICS.md` or `GROK_CORRECTIONS.md`.
**Depth:** CODEX_RESEARCH_NOTES.md is 13,063 lines. I read its section index for 2026-09-30 to 10-13 (lines 8515–13061) and the 10-13 parity note in full, not every line. CODEX_HEURISTICS.md (995 lines) and rules/GPPVERIFY.md were skimmed only.

## 1. Current front
- **The frame (CLAUDE_RH_MAP.md:10, :66).** Every route is a statement of the form "W ≥ 0 (or tempered) on some family F". The missing piece is always the same: show the Euler product forces W ≥ 0 on F. Davenport–Heilbronn (DH) is the falsifier for every route.
- **Most active routes:**
  - R1, fixed window / two box (:30). Missing: an arithmetic upper bound on S_ℓ(t).
  - R2/R17, the semilocal odd Weil matrix and the CCM finite core (:31, :45). Missing: the closure step, ⟨k,Tk⟩·𝔐 → 0, with the η channel kept open.
  - R18, Koebe depth moments (:46).
- **Live work (10-13 to 10-14):**
  - Arb-certified compact-window negativity for F_{1/4} and F_{1/10}, with ζ positive on the same windows (CLAUDE_CODE_RESEARCH_NOTES.md:2006-2018).
  - The θ-expansion Q_θ = 2Q₀(F) + θ²Q₀(t²F) + O(θ⁴). Q₀(t²F) < 0 is certified on the witnesses (:2024-2026; work/codex/2026-10-13_zero_margin_theta_curvature_parity_response.md:30-36).
  - A DH single-carrier witness, g = u(1−u²)^k cos(857x/10). Arb gives −0.64137 at h=3, k=3 and −1.11918 at h=7/2, k=4 (CLAUDE_CODE_RESEARCH_NOTES.md:2034-2040).
  - Turbo's counting law: bundle = D − 2N(Ω) − 2 (CONVERSATION.md, Turbo 10-13 counting entry). This is numerical only.
  - The new parity pole-term observation (CLAUDE_CODE_RESEARCH_NOTES.md:2116-2117; CONVERSATION.md:737).
- **Formalization in the last two days:** a run of Lean merges, #326 to #340, all standard axioms only (CONVERSATION.md:~700-735). These are structural or finite identities. None is a positivity mechanism.

## 2. Dead routes (CLAUDE_CODE_CORRECTIONS.md, CLAUDE_CORRECTIONS.md)
- **Local prime-kernel positivity (K_p−1) as global progress.** The prime term enters the form with the opposite sign (CC_CORR:30-47).
- **Pointwise multiplier sign from convolution positivity.** −W_p takes both signs (CC_CORR:71-87).
- **Fixed-Mellin Γ-prime contraction, or a finite-prime Koszul homotopy.** The kernel is indefinite and small divisors break it. Recorded as "do not revive" (CC_CORR:89-104).
- **"All three dead routes share one invariant."** This holds for 2 of the 3. Hardy reflection died of Fredholm index −1 ≠ +1, while the GCD generator B and co-Poisson compression died of N(T) = e^T counting rigidity (CC_CORR:603-620).
- **"Compensation breaks ζ positivity."** The compensated form is not a Weil form (CC_CORR:708-710).
- **"Pointwise Widder tests can't detect F_θ."** Wrong: they detect it at Hankel size ≥ 4–5 (CC_CORR:714-726).
- **R_N → 0 at fixed λ.** False: the values fall to a minimum and then grow (CLAUDE_CORRECTIONS.md:5).
- **"UV decoupling controlled".** The UV term is exactly the unresolved-zero tail (:6).
- **The boundary Wronskian identity.** It just restates βᵀξ = 0 (:7).
- **"Positivity selects the primes" / prime-blind archimedean positivity (E).** An SDP finds feasible b with margin 6e-3 far from Λ (:24).
- **"Wider windows expose DH sooner."** False: short windows win (:26).
- **"The open step is bounded ⇒ tempered".** The fixed-window paper already proves that step. The open step is the arithmetic bound (:28).
- **"One negative square per off-axis quartet."** It is two, signature (2,2) (:30; Lean QuartetSignature #283).
- **Turbo's "F_θ defect purely odd".** Retracted by Turbo (CONVERSATION.md, Turbo 10-13 "Parity correction ACCEPTED").
- **Also dead (CLAUDE_RH_MAP.md:81).** A one-sided bound with real λ and a negative coefficient is not enough, because −e^s ≤ 0.

## 3. Problems found
1. **The parity pole term: the algebra is right, the conclusion overreaches.**
   - *Check.* Take real g, G(s) = ∫g(x)e^{sx}dx and f = g * g̃. Then f's Mellin values at the two poles each equal G(½)G(−½). Write C = ∫g cosh(x/2) (even part only) and S = ∫g sinh(x/2) (odd part only). The pole term is 2G(½)G(−½) = 2(C² − S²). That gives +2G(½)² for even g and −2G(½)² for odd g. The even and odd parts do not mix. This is consistent with #264 SheetParitySplit, "sinh pairs only with the odd part" (CLAUDE_CODE_COVERAGE.md:309). The pole enters with a + sign on the zero side, where Σ_ρ = pole − arch − prime.
   - *Hidden assumptions.* The statement assumes real g and needs G(½) defined. For complex g the pole term is G(½)·conj(G(−½))-type, not a square. The statement also changes for F_θ, whose poles sit at {±θ, 1±θ} (CLAUDE_CODE_RESEARCH_NOTES.md:2118), so the parity split does not carry over to the falsifier as written.
   - *The overclaim.* "The odd sector is the only place where Arch+prime must dominate" (CONVERSATION.md:737; CLAUDE_CODE_RESEARCH_NOTES.md:2117) does not follow. In the even sector you still need −(arch + prime) ≥ −2C², which is a real constraint. Claude Code's own certificate says F_{1/4} fails in BOTH sectors, the even one first at h ≈ 0.50, and "not pole-driven" (CLAUDE_CODE_RESEARCH_NOTES.md:2021). Codex says the same (parity note:44, "refutes ... odd-only defect"). Line 2117 also calls the split "consistent with Turbo's 'F_θ defect purely odd'", which Turbo has retracted.
   - *Status.* PROVED identity, misread consequence. It should be relabelled: the odd sector is where the pole works *against* positivity. It is not the only sector where the arithmetic must win.
2. **The earlier Grok stint conflicts with the current seat map.** GROK_RESEARCH_NOTES.md:34-35 says Grok created GPPDiscovery3 and `GPPVerify:grok/workbench`. GROK.md:3, :15 now say GPPDiscovery3 is Astra's and off-limits.
   - Grok's old script `weil_heat/prime_heat.py` still sits in GPPDiscovery3.
   - The branch was renamed `guest/workbench` (rules/GPPVERIFY.md:725). It still holds `ArithmeticWeilHeat.lean`, which is not on GPPVerify main (404).
   - PR #187 is CLOSED and unmerged. GROK_RESEARCH_NOTES.md:56-68 still describes it as pending CI, with a possible `vonMangoldt_nonneg` name fix.
   - Open items left from that stint: that file's Lean content (prefactor positivity, even collapse, primeSide ≥ 0, heatTrace ≥ 0 ↔ arch dominance) was never merged or audited. The notes stop on 2026-09-13; there is no closing entry.
3. **Stale or contradictory pointers.**
   - CONVERSATION.md (10-14 first note) says private work goes to `work/guest/` and `GUEST_RESEARCH_NOTES.md`. It was corrected in the next note to `work/grok/` and GROK_*.
   - The same note describes the Discovery5 layout as research/experiments/notes/scripts. Daniel specified discovery/<thread>/, research/, experiments/, scripts/.
   - rules/PROJECT.md is headed "Claude operating instructions" (:1). It points at `/home/user/gpp-bridge` (:100) and tells agents to query `gpp_vault` at session start (:59-63). Grok has no vault access and must not read credentials.
   - CLAUDE_CODE_HEURISTICS.md:1 is titled "Claude heuristics", left over from the seat rename.
   - CLAUDE_CODE_COVERAGE.md:3 says "Updated 2026-10-03", and its counts (:8-16) are stale. The table has rows through 10-14 (:362-374).
   - The README says the bridge is a full ledger with Supabase as the mirror (:43-46). PROJECT.md:191 says the bridge "is not the sole research ledger". That is a mild inconsistency.
4. **Status labels.**
   - R14 relies on `True := trivial` placeholder stubs (CLAUDE_RH_MAP.md:43). These are honestly flagged.
   - R1's "zeta holds for all integers x ≤ 1e7" (:30) is floating-point evidence. It sits in a row whose status is "paper".
   - I found no THEOREM label attached to a purely numerical claim in the files I read.

## 4. Where GROK.md differs from Daniel's instructions
- **Mirror.** GROK.md:10 asks for a periodic Supabase `guest.research_notes` checkpoint. Daniel says Grok has no access yet, and to checkpoint only if he grants it.
- **Lab content.** GROK.md:7 bans "process files" from the lab. Daniel wants a README, requirements.txt and a GitHub Actions runner there. Those count as reproducibility files, not process files, so I see no real conflict, but it's worth stating.
- **Channel.** GROK.md:11 says to post to CONVERSATION.md "only when another seat can act". Daniel's message says to post drafts to him.
- **Missing from GROK.md:** the falsifier rule, the status labels, the hand-off format, the credit line and the "full book title" rule. GROK.md's boundaries do match Daniel's list: no pushes to the other seats' repos or branches (:15), and no edits to the CLAUDE_*, CLAUDE_CODE_*, CODEX_* or TURBO_* files (:16).
- **Lean.** GROK.md:8 and Daniel agree: Grok does not push to GPPVerify unless Daniel asks.

## Instructions found inside the data
None addressed to Grok beyond the seat rules above.

## 5. PR #11 (GPP-bridge chat channel)
1,175 issue comments, all posted from the GoldenPhysicsProject account and prefixed with the sending seat. There are 0 review comments.
- **Addressed to Grok:** none. No comment mentions Grok or the guest seat.
- **Recent traffic (10-13 18:48 to 10-14 03:14 UTC).** Codex posts a note built on Daniel's physical pictures, and Claude Code audits it within minutes ("everything checks") and merges the exact parts into GPPVerify. Topics in order:
  - prime-Gibbs/Kakutani all-prime detector, #332–#334;
  - Bost–Connes and the conductor-5 Gauss sum, #335;
  - prime fermion / theta mirror, #336 (with the supertrace correction);
  - zero↔prime hierarchy;
  - Turbo's θ-branch sewing test (NUMERICAL), where F_θ's defect is "sewn-only" and Claude Code adds a caveat;
  - the conserved-(q,t) entangled pair;
  - the mirror Gram lemma, overlap 2hδ/sinh(2hδ), #338;
  - the Cayley circle, #339, and the celestial sphere, #340;
  - the arithmetic-seed Euler horodisk, #341;
  - quartet/cosh Cayley geometry, #342 (CI running).
  - Every note is labelled "not RH".
- **Source of the parity claim:** comment 6159594763 (Claude Code, 10-14 02:33), point 9 onward.
  - Convention: ĝ(½+z) = G(z) = ∫g e^{zt}, with g real and f = g * g̃.
  - Pole term: f̂(0) + f̂(1) = 2ĝ(1)ĝ(0) = 2G(½)G(−½). This gives +2G(½)² for even g and −2G(½)² for odd g.
  - **The derivation is internally consistent.** It agrees with my C² − S² check in §3.1.
  - **Two problems in the same comment:**
    1. It cites Turbo's "purely odd F_θ defect" as support, but Turbo had already retracted that claim.
    2. "The even sector has the pole term on its side" plays down that the even sector still needs −(Arch + prime) ≥ −2G(½)². F_{1/4} fails there first (CLAUDE_CODE_RESEARCH_NOTES.md:2021).
  - The comment itself does say F_θ's pole term does not split by parity, and that "odd-sector domination is not easier than RH". That is correct and is the honest reading.
- **Companion comment 6159500176:** the pole term is a cross-face pairing 2Re[ĝ(1)·conj ĝ(0)], which is indefinite. DH has no pole at all, so a pole-based mechanism cannot explain DH's failure.
- **Pattern worth noting, as a skeptic:** audits come back "all correct" within 1–3 minutes, for every note. What the audits verify are identities. None supplies a positivity mechanism, so the Lean count grows while the RH front itself does not move.

## 6. GPPVerify survey (shallow clone of main, read-only)
- **Setup.**
  - Toolchain: leanprover/lean4:v4.33.1, built with lakefile.toml. The tree has 1,017 `.lean` files and 432 of them are under `GppVerify/RiemannHypothesis`.
  - Top-level folders: Analysis, CelestialHolography, Cosmology, GeneralRelativity, NoGo, NumberTheory, QFT, QuantumGravity, QuantumInformation, RiemannHypothesis, StandardModel, StringTheory, YangMills, Thread{HT,S,WeilParity,WeilSemibound}, Upstream.
  - README.md:1-6 says plainly that it "is not a proof of any famous conjecture".
- **CI.**
  - main Build: success on #340 (run 38503290294) and #339.
  - PR #341 CayleyEulerHorodisk and PR #342 LiQuartetMoment: Build in progress, Blueprint green.
- **Gate grep.**
  - `sorry`/`admit`: only in comments and docstrings. No live hits.
  - `axiom` declarations: 0.
  - **`native_decide`: 99 uses in 18 files**, for example ShadowEulerIdentity.lean:160-168, YangMills/MassGap.lean (7) and StandardModel/UnitarySheetRelationalMatter.lean (8). It relies on the `Lean.ofReduceBool` trust axiom. The README and notes claim "0 axioms" / "standard axioms only" (bridge README:112). That is true only if `ofReduceBool` is not counted. `scripts/check_axioms.lean` never mentions it. **Flag this to Claude Code.**
  - `True := trivial` placeholders: 112. These are known, named open stubs (CLAUDE_RH_MAP.md:43).
- **RH line.**
  - `def RiemannHypothesis` (NumberTheory/ZetaProperties.lean:131) is used only as a hypothesis. No theorem proves it.
  - The explicit-formula pieces are separate modules (ArchimedeanFloor, ArchimedeanTailMass, WeilSupportLadder, ...). I found no single module proving the full Weil explicit formula for ζ.
  - Pole and parity modules: SheetParitySplit.lean:97 `pole_pairing_odd_only` (the sinh pairing sees only the odd part), PoleNeutralCongruence (#299), PoleNeutralFilter, PoleNeutralMultiplier, RadicalPoleThreshold.
  - **The ±2G(½)² parity pole statement itself is not in Lean.** It exists only in PR #11 comment 6159594763.
  - **The DH zero detection is not in Lean.** It is an Arb/Python certificate in bridge `work/claude-code/2026-10-13_dh_modulated_window/`. In Lean, DH appears only in docstrings, plus NoGo/Ledger.lean:54, which still says "needs DH zero computation". That line is stale relative to the Arb certificate.
  - Recent claims I checked against the tree are all present: MirrorGramWindow, CayleyCircleHarmonic, CelestialCayleySphere, QuartetSignature, OddCollisionTrace, PrimeScatteringPhase, ConductorFiveMellin.
- **Cross-check against CLAUDE_CODE_COVERAGE.md.**
  - I found no claimed-done row whose module is missing.
  - The counts header (:8-16, "Updated 2026-10-03") is stale.
  - The old Grok `ArithmeticWeilHeat.lean` lives only on branch `guest/workbench` (closed PR #187), not on main. It has no coverage row, so it is formalized-but-unmerged and unaudited.
