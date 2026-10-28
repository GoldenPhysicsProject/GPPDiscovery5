# ARP §S — Spectrum & Measure: the Weil distribution as an indefinite (Kreĭn/Pontryagin) measure

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Spectrum & Measure seat, 2026-10-28. Companion to `derivation/ARP.md` (origin/main 36aa6be), assignment row §S:
"Bochner for the mean field (§G.3) and the signature of the fluctuation."

**RH is OPEN. Nothing in this file proves or implies RH. No step below claims progress on RH.**
No numerics were run for this file. All statements are PROVED (citation or a short proof written out here),
PROVED-LEAN (GPPVerify name), CONJECTURE, or GAP, as in ARP §0.

---

## S.0 Conventions (fixed once; every constant below uses them)

- Window: f ∈ C_c^∞([−L, L]) (complex-valued allowed), f̃(u) := conj f(−u), g := f ∗ f̃, so supp g ⊂ [−2L, 2L].
- Fourier: f̂(z) := ∫ f(u) e^{izu} du (ARP §G, Suzuki 2206.03682 §3, Bombieri 2000 via Suzuki (5.15)).
  Then ĝ(z) = f̂(z)·conj f̂(z̄); on ℝ, ĝ = |f̂|² ≥ 0. Plancherel: ∫|f|² du = (1/2π)∫|f̂|² dt.
- Zeros: ρ = ½ + iγ_ρ over nontrivial zeros, with multiplicity. Q_L(f) := Σ_ρ ĝ(γ_ρ) (symmetric limit).
- **Alignment with main's μ_ρ (G.4d, 7328e39).** Main writes ρ = ½ + β + iγ (γ real) and μ_ρ := harmonic (Poisson-kernel) measure
  of ρ relative to the critical line, dμ_ρ(t) = (|β|/π) dt/((t − γ)² + β²), a probability measure, μ_ρ = δ_γ if β = 0. In §S
  notation γ_ρ = γ − iβ (complex), and the off-line term in Q_L is the exact point value m_ρ ĝ(γ_ρ) (explicit formula, PROVED),
  not ∫ĝ dμ_ρ. Exact split: Q_L(f) = Σ_ρ m_ρ ∫ĝ dμ_ρ + Σ_{β≠0} m_ρ [ĝ(γ_ρ) − ∫ĝ dμ_ρ]. The first sum is ≥ 0 (ĝ = |f̂|² on ℝ, μ_ρ ≥ 0);
  the second is minus main's "off-line residual". μ_ρ = μ_{1−ρ̄} (same Poisson kernel), so the residual pairs as {ρ, 1−ρ̄} like
  the signature-(1,1) pair form of S(b)3.6, and κ_Q(L) ≤ (negative index of the residual) since the first sum is PSD. We use only this algebraic split. The
  limit identity lim_Λ Δ_Λ = Σ_ρ m_ρ ∫ĝ dμ_ρ is PROVED for function fields only (Connes 1999 §VIII Lemma 3) and is CONJECTURE for ℚ
  (G.4d); §S does not use it.
- Explicit formula (Weil 1952; in this normalization Suzuki eq. (5.15), derived there from Bombieri 2000 p. 186):
  Q_L(f) = ĝ(i/2) + ĝ(−i/2) − Σ_{n≥2} Λ(n) n^{−1/2} [g(log n) + g(−log n)] + W_∞(g),
  W_∞(g) = −(log π) g(0) + (1/2π) ∫ Re ψ(¼ + it/2) ĝ(t) dt.
- Weil criterion (PROVED: Weil 1952; Yoshida 1992; Bombieri 2000; stated in this convention in Suzuki §3.2 eq. (3.4)):
  RH ⇔ Q_L(f) ≥ 0 for all f and all L. PROVED-LEAN (finite, free-coefficient form):
  `GppWeilCriterion.rh_iff_weil_pairedForm_nonneg` (GPPVerify `RiemannHypothesis/WeilPositivityCriterion.lean`).
- **Fluctuation measure: E_Λ, not E₁.** On (0, ∞):
  E_Λ := Σ_{n≥2} Λ(n) n^{−1/2} δ_{log n} − e^{u/2} du  (half-density of ψ(x) − x).
  Reason (ARP G.3d, repaired 36aa6be, PROVED, classical): ψ(x) = Σ_{k≥1} θ(x^{1/k}), so the prime-power layers
  μ_k = Σ_p (log p) p^{−k/2} δ_{k log p} sum to μ_Λ = Σ Λ(n) n^{−1/2} δ_{log n}, the half-density of ψ, whose PNT mean is the
  single density e^{u/2} du. There is no separate k-th layer mean to compensate (in particular no k = 2 rank-one).
  The old "term-by-term" cancellation via θ = Σ μ(k) ψ(x^{1/k}) is false (k = 4: μ = 0; k = 6: μ = +1); we do not use it.
- Decomposition: Q_L = MF_L + Fl_L with
  MF_L(f) := [ĝ(i/2) + ĝ(−i/2)] − ∫_0^∞ e^{u/2}[g(u) + g(−u)] du,
  Fl_L(f) := W_∞(g) − ∫_{(0,∞)} [g(u) + g(−u)] dE_Λ(u).
  (Exact; the pole is spent once, in MF.)

---

## S(a) The mean field: Bochner measure, with constants (positive part)

**S(a)1 (PROVED, elementary; = Astra's Ward identity G.3a at v = 0).**
MF_L(f) = ∫_ℝ g(u) e^{−|u|/2} du, for every f (no evenness needed).
*Proof.* ĝ(i/2) + ĝ(−i/2) = ∫ g(u)(e^{u/2} + e^{−u/2}) du. Split ℝ at 0: the u > 0 and u < 0 halves give
∫_0^∞ [g(u) + g(−u)](e^{u/2} + e^{−u/2}) du; subtracting ∫_0^∞ e^{u/2}[g(u) + g(−u)] du leaves ∫_0^∞ [g(u) + g(−u)] e^{−u/2} du. ∎
Matches G.3a (Astra): pole − PNT mean = ∫ G₀(w) e^{−|v−w|/2} dw, with coefficient exactly 1.
R1 uniform in L (identity for every L) | R2 uses pole + PNT mean only, no FE/Euler product jointly | R3 no prime-weight sign | R4 contains the pole state (below) | R5 see table S.4.

**S(a)2 (PROVED, elementary; Cauchy/Poisson kernel).** With the S.0 Fourier convention,
e^{−|x|/2} = ∫_ℝ e^{−itx} dμ_MF(t),  **dμ_MF(t) = dt / (2π (t² + ¼))**,
a positive finite measure of **total mass 1** (= κ(0), since (1/2π)∫dt/(t²+¼) = (1/2π)·2π = 1).
Equivalently ∫ e^{−|x|/2} e^{itx} dx = 1/(t² + ¼) = 4/(1 + 4t²). Hence (Bochner, exact constants)

  **MF_L(f) = ∫ ĝ dμ_MF = (1/2π) ∫_ℝ |f̂(t)|² / (t² + ¼) dt.**

- Constant/convention statement (answers Golden Eye's flag (1)): density **1/(2π(t²+¼))** with respect to dt, for
  f̂(t) = ∫ f e^{itu} du; κ(u,v) = e^{−|u−v|/2} with κ(0) = 1. In the unitary convention F f(t) = (2π)^{−1/2}∫ f e^{itu} du the
  same form is ∫ |Ff|² /(t² + ¼) dt, density 1/(t²+¼) (ARP G.3a's wording). **No ½-vs-¼ ambiguity exists in the MF itself.**
- Bounds (PROVED, from the formula): 0 < MF_L(f) ≤ 4 ‖f‖₂² for f ≠ 0 (density > 0 and f̂ entire ≢ 0; sup of 1/(t²+¼) is 4).
  So MF_L is a strictly positive, bounded, positive-definite form on every window; as an operator on L²[−L, L] it has the
  continuous positive-definite kernel e^{−|u−v|/2}, hence is trace class (Mercer).
- OS square (ARP G.3a, PROVED): ∫∫_{u,v≥0} f(u) f(v) e^{−(u+v)/2} = |∫ f e^{−u/2}|²: rank one, vector = pole state e^{−u/2}.
R1 every L, constants L-free | R2 none needed (no RH content) | R3 yes | R4 the pole state and the soft mode are in the MF-positive space (MF > 0 on all f ≠ 0) | R5 S.4.

**S(a)3 Normalization ½ vs ¼ (Golden Eye flag) — reconciled for the Toeplitz part; moot after G.3d.**
The flag concerns the k = 2 layer (ARP G.3d's "−½ ĝ(0)" vs Astra's even bulk "−¼·11ᵀ"), not the MF constant.
- (PROVED, elementary.) k = 2 mean density ½ du on (0, ∞) acting on g(u) + g(−u) gives −½∫_ℝ g = −½ ĝ(0) = −½|f̂(0)|².
  With Astra's bumps b(u) = sin²(πu) 1_{[0,1]}(u), ∫b = ½, f = Σ c_i b(· − a_i): |f̂(0)|² = ¼(Σ c_i)², so the term is
  −⅛ (Σ c_i)² = −⅛·11ᵀ, which is exactly Astra's **Toeplitz** component (bridge `work/2026-10-27_astra_prime_pnt_os_ward_identity.md`).
- Astra's −¼·11ᵀ is the sum Toeplitz (−⅛) + Hankel (−⅛) in the even sector of her doubled (OS) bookkeeping. The Hankel −⅛ is
  taken from Astra's note, **not re-derived here**.
- Status: ½ (full-line, per |f̂(0)|²) and ¼ (even sector, Toeplitz+Hankel, per 11ᵀ with ∫b = ½) are the same term in two
  bookkeepings for the Toeplitz half: PROVED. **GAP (normalization, minor): independent replay of Astra's Hankel −⅛.**
  Either way the term has no role after G.3d (E_Λ carries no k = 2 remainder).

---

## S(b) The fluctuation as a signed / Kreĭn measure; negative index via zero counting

**S(b)1 Spectral side of Fl (PROVED, from S.0 + S(a)2).**
Fl_L(f) = Σ_ρ ĝ(γ_ρ) − ∫ ĝ dμ_MF. I.e. the Fourier-side "measure" of the fluctuation is
  **ν_Fl = Σ_ρ δ_{γ_ρ} − μ_MF**  (acting on ĝ, ĝ entire of exponential type 2L).
- Under RH: ν_Fl is a real signed measure on ℝ with Jordan parts ν₊ = Σ δ_γ (discrete) and ν₋ = μ_MF (absolutely continuous,
  mass 1). Note ν₋ is **not** the mean zero density (that is (1/2π) log(|t|/2π) dt, which W_∞'s Re ψ(¼+it/2) carries).
- Off RH: an off-line pair γ, γ̄ (i.e. ρ and 1 − ρ̄) contributes ĝ(γ) + ĝ(γ̄) = 2 Re[f̂(γ)·conj f̂(γ̄)], a Hermitian form of
  signature (1, 1): ½(|A+B|² − |A−B|²) with A = f̂(γ), B = f̂(γ̄). This is the mechanism formalized (finite, free coefficients) in
  `GppWeilCriterion.rh_iff_weil_pairedForm_nonneg` (pairedForm Σ conj c(1−ρ̄)·c(ρ)), PROVED-LEAN.
- **Absorption sign dictionary.** In Connes' trace-formula convention the zeros enter with a minus sign: the critical zeros are an
  *absorption* spectrum, off-line zeros appear as resonances (Connes 1999, Selecta 5, abstract), and RH ⇔ Σ_v W_v ≤ 0 on the
  test functions with ĝ(±i/2) = 0 (Connes–Consani 2021, Selecta 27:77, introduction). In ARP's convention the same identity reads
  Q_L = Σ_ρ ĝ(γ_ρ) ≥ 0. For Fl: −Fl_L = μ_MF(ĝ) − Σ_ρ ĝ(γ_ρ) = "continuum minus absorption lines". Same identity; only the
  overall sign convention differs, and we state everything in ARP's sign.
R1 identity for each L | R2 this is just the explicit formula | R3 — | R4 — | R5 fake-blind identity (any form with a pole and a zero set).

**S(b)2 Window indices (definitions).** κ_Q(L) := negative index (max dim of a negative subspace) of Q_L on C_c^∞([−L, L]);
κ_Fl(L) likewise for Fl_L; κ_Q(∞) := negative index of Q on ∪_L C_c^∞([−L, L]); N_{>½} := #{distinct zeros with Re ρ > ½}
(= number of off-line pairs {ρ, 1−ρ̄} counted once; two per quadruple).

**S(b)3 What is known at window L (monotone, finite, eventually global).**
1. *Monotone* (PROVED, elementary): L ≤ L′ ⇒ κ_Q(L) ≤ κ_Q(L′), κ_Fl(L) ≤ κ_Fl(L′) (test-space inclusion).
2. *Small windows* (PROVED): κ_Q(L) = 0 for L small — Yoshida 1992 Lemma 2 (positivity on K(a), a small; as reported in
   Suzuki 2206.03682 §4.2); Bombieri 2000 reproves it; for L < (log 2)/2 no prime enters (PROVED-LEAN
   `GppWeilLadder.primeSide_eq_zero_of_support_lt_log_two`) and Connes–Consani 2021 give the operator-theoretic reason (Thm 1).
3. *Finite at each L* (PROVED): Yoshida 1992 Lemma 3 (as stated in Suzuki §4.2): for a₀, μ > 0 there is N with
   Q ≥ μ‖·‖² on the finite-codimension subspace K_N(a), 0 < a ≤ a₀ (Suzuki's transcription of the norm power not re-checked; the
   inequality is read as homogeneous of degree 2). Since C(a) ⊂ K(a), C(a) ∩ K_N(a) has codim ≤ 2N+1 in C(a). Hence κ_Q(L) ≤ 2N+1 < ∞. Taking μ > 4 ≥ ‖MF_L‖ (S(a)2)
   also gives Fl_L ≥ (μ − 4)‖·‖² > 0 on K_N, so **κ_Fl(L) < ∞** for each L.
4. *Ordering* (PROVED, elementary): Q_L = Fl_L + MF_L with MF_L ≥ 0 ⇒ every Q-negative subspace is Fl-negative ⇒ κ_Q(L) ≤ κ_Fl(L).
5. *Eventually global* (PROVED, elementary): any finite-dim negative subspace of Q lies in some window, so
   κ_Q(∞) = sup_L κ_Q(L) = lim_L κ_Q(L); if κ_Q(∞) < ∞ then κ_Q(L) = κ_Q(∞) for all L ≥ some L*. L* is not effective.
6. *Upper bound by zero counting* (PROVED, elementary; written out here): κ_Q(∞) ≤ N_{>½}.
   Proof: Q(f) = Σ_{γ real} |f̂(γ)|² + Σ_{pairs} m_ρ·2Re[f̂(γ)conj f̂(γ̄)] (absolutely convergent for f ∈ C_c^∞: f̂ decays rapidly in
   the strip |Im z| ≤ ½). The first sum is ≥ 0; each pair term is a form of signature ≤ (1,1) (multiplicity m_ρ does not raise its
   rank). A form = (PSD) + (sum of k forms of negative index ≤ 1) has negative index ≤ k. ∎ (Holds verbatim for fakes: fake-blind.)
7. *Lower bound / exact count* (PROVED in Bombieri's setting): **Bombieri 2000, verbatim (abstract):** "if the Riemann Hypothesis is
   false but only with finitely many non-trivial zeros off the critical line we show that the number of negative eigenvalues is
   precisely one-half of the number of zeros failing to satisfy the Riemann Hypothesis, provided the truncation is big enough."
   One-half of all off-line zeros = N_{>½} (counting as in Bombieri; his multiplicity convention not re-checked — full text
   inaccessible this session). Also (Bombieri, intro): 0 is not an eigenvalue of the finite approximations, and the finite-dim
   results carry over to the infinite-dim case for negative eigenvalues bounded away from 0.
   Unconditionally: N_{>½} ≥ 1 ⇒ κ_Q(∞) ≥ 1 (Weil). Transition points: if RH fails, the first L where κ_Q jumps is a
   degeneracy point (Yoshida 1992 Thm 2 / Suzuki 2206.03682 Thm 1.4 and Thm 6.1 for the screw-function form).
   **GAP (verification/transfer):** identify Bombieri's "truncation" with the support window [−L, L] and conclude
   κ_Q(L) = N_{>½} for L ≥ L*. Opened in GAP2_bombieri_window.md (this directory): the window statement is proved there under
   "finitely many off-line zeros", with N_{>½} counting distinct zeros; only the identification of Bombieri's own truncation stays GAP.
8. *Kreĭn–Langer side* (PROVED, literature): Kreĭn–Langer 1977 (Math. Nachr. 77, 187–236), Satz 5.9 — bijection between screw
   functions g ∈ 𝒢_∞, g(0)=0, and Nevanlinna Q with Q(iy)/y → 0 (as cited by Suzuki, Thm 1.2 there: RH ⇔ −Ψ is a screw function).
   The indefinite (𝒢_κ / generalized Nevanlinna N_κ) count "κ = #poles of i ξ′/ξ(½−iz) in ℂ₊" for ζ is the S97 claim
   (`KREIN_LANGER_S97.md`); **not re-verified in a source this session → GAP (citation)**; and per S97 §B a finite-window
   extension certificate does not constrain κ(∞) (no RH leverage).
- **Kreĭn/Pontryagin structure** (PROVED, standard): by 3, (C_c^∞([−L, L]), Q_L) modulo null vectors completes to a Pontryagin
  space Π_{κ_Q(L)} (Iohvidov–Kreĭn–Langer 1982, general theory; Yoshida: Q nondegenerate on C(a)). RH ⇔ κ_Q ≡ 0.
R1 items 1–5 are per window; only 5–7 touch L → ∞, and none is uniform with an effective L* | R2 item 7 uses ζ's zeros (both) | R3 no
weight signs | R4 — | R5 items 1–6 fake-blind (DH: infinitely many off-line zeros ⇒ κ_Q(∞) may be ∞ while each κ_Q(L) is finite).

**S(b)4 The fluctuation alone has negative index even under RH (PROVED, unconditional, non-effective; written out here).**
There is L₀ such that **Fl_L is not positive semidefinite for every L ≥ L₀.** Under RH moreover κ_Fl(L) → ∞.
*Proof (case split).* Fix φ ∈ C_c^∞([−1, 1]) real, ‖φ‖₂ = 1 (nodeless allowed), and f_L(u) := L^{−1/2} φ(u/L), supported in
[−L, L], ‖f_L‖₂ = 1, f̂_L(t) = L^{1/2} φ̂(Lt).
- MF: MF_L(f_L) = (1/2π)∫ |φ̂(s)|² / ((s/L)² + ¼) ds ↑ (4/2π)∫|φ̂|² = 4 (monotone convergence; Plancherel).
- If RH holds: Q_L(f_L) = Σ_γ L|φ̂(Lγ)|² ≤ C_φ² L^{−3} Σ_γ |γ|^{−4} → 0 (|φ̂(s)| ≤ C_φ (1+|s|)^{−2}, all γ real, |γ| > 14,
  Σ|γ|^{−2} < ∞). Hence Fl_L(f_L) = Q_L(f_L) − MF_L(f_L) → −4 < 0. Taking k L²-orthonormal such φ_j, the k×k Gram of Fl on
  span{f_{j,L}} tends to −MF-Gram → −4·I_k, so κ_Fl(L) ≥ k for L large.
- If RH fails: by Weil, Q_L takes a negative value for some L₁, hence for all L ≥ L₁ (item 1); Fl_L ≤ Q_L.  ∎
R1 statement is for all large L, but L₀ is non-effective | R2 not used (only: zero set discrete, Σ|γ|^{−4} < ∞) — **fake-blind, no
RH content** | R3 — | R4 the witnesses are dilated nodeless bumps: the soft mode (S97/S98) is exactly where Fl ≈ −MF | R5 S.4.

---

## S(c) Why Bochner on the mean field cannot extend to the fluctuation

**S(c)1 (PROVED, S(b)4).** Bochner's theorem characterizes positive-definite distributions as Fourier transforms of *positive*
measures. ν_Fl has the negative part μ_MF on all of ℝ even under RH (S(b)1), and Fl_L is not PSD for L ≥ L₀ (S(b)4). So
"Bochner for Fl" is false, unconditionally.

**S(c)2 (PROVED: Weil; Yoshida; Suzuki §3.2).** Bochner for the *sum* MF + Fl = Q is Weil's criterion: Q is positive definite
on all windows ⇔ Σ_ρ δ_{γ_ρ} is a positive measure on ℝ (all γ real) ⇔ RH. So any "extension of Bochner from MF to Fl" is either
false (Fl alone) or a restatement of RH (Fl + MF) — the §0.2 "Schur/Kreĭn/de Branges restatement" dead end.

**S(c)3 Off-line atoms are not Bochner-admissible (PROVED, elementary; PROVED-LEAN analog).** For Im γ ≠ 0, x ↦ e^{iγx} is not
positive definite (|k(x)| ≤ k(0) fails). Continuing the positive Cauchy/Abel character off its strip produces indefinite Gram
matrices: `GppAbelContinuedKernel.abelCont_gt_one`, `GppAbelContinuedKernel.contForm_indefinite` (eigenvalues 2/3, −4/3), and
`GppAbelContinuedKernel.no_invariant_mean` (GPPVerify `RiemannHypothesis/AbelContinuedKernelIndefinite.lean`, Astra source; its
own scope note: says nothing about ζ's zeros). Same lesson for μ_MF: its positivity is a property of the integral against
1/(t²+¼) on ℝ and does not transport to complex atoms.
*Remark (PROVED, elementary).* A positive "harmonic-measure" smearing of an off-line atom onto the line (cf. Connes 1999 Lemma 3,
function-field setting) does not reproduce ĝ(γ_ρ): Poisson reproduction needs bounded harmonic functions on the half-plane, and
ĝ (exponential type 2L) is unbounded there. So positivity of such smeared measures gives no bound on Q.

**S(c)4 No slack (PROVED, unconditional; written out here).** Let r(L) := inf_{f ≠ 0} Q_L(f)/MF_L(f) (well-defined: MF_L > 0).
r is nonincreasing in L, and **lim_{L→∞} r(L) ≤ 0**; RH ⇔ r(L) ≥ 0 for all L ⇔ lim r(L) = 0.
*Proof.* Monotone by inclusion. RH case: r(L) ≤ Q_L(f_L)/MF_L(f_L) → 0/4 (S(b)4). Non-RH case: r(L) < 0 for L ≥ L₁. ∎
Consequence (PROVED; scope narrowed after Golden Eye's audit, which found S(c)4 itself sound): **no MF-relative margin** — there is
no c > 0 with Fl_L ≥ −(1−c)MF_L for all L (equivalently no c > 0 with Q_L ≥ c·MF_L for all L). So any argument of the form "MF
dominates the fluctuation by a fixed factor" is ruled out, unconditionally.
*What is NOT proved.* This says nothing about Q_L ≥ c·N_L(f) for a different norm N_L that degenerates on the dilated nodeless bumps
(e.g. one weighting only frequencies near the zeros, such as N_L(f) = Σ_γ w(γ)|f̂(γ)|²-type expressions); S(c)4 does not exclude such
a margin.
**CONJECTURE (well supported; not proved):** the closing step for G.4d must be an exact identity (zero margin in every natural norm),
rather than an inequality with slack. Support: S(c)4 (MF-relative), S98 (fake-blind positive cone), §0.2 restatement dead ends.
Credit: Golden Eye audit of S(c)4 (2026-10-28).
The MF-relative saturating directions are the dilated nodeless bumps (R4 soft mode; consistent with S98's fake-blind positive cone and the project's tiny window λ_min,
which is NUMERICAL orientation only).
R1 this is precisely the R1 obstruction made quantitative | R2 not used — fake-blind | R3 — | R4 saturation on the soft mode | R5 S.4.

**S(c)5 Cokernel sign (SETTLED on main: derivation/ARP.md G.4d, commit 7328e39, after Golden Eye's three-test sign audit, SOUND on
all three).** The positive object is the **deficit on coker E = ker E\***, i.e. a trace over the complement of range(E)
(S_Λ − Q′, which is free), not ‖P^⊥Eφ‖². Range positivity enters Δ_Λ with a minus sign and gives only the trivial upper bound, so
"additive positivity transfers through E" is WRONG-SIGN and withdrawn (G.4d). The two proved instances agree with this sign:
Connes 1999 Thm 5 (positive characteristic; proof uses Δ_Λ = Tr((S_Λ − Q′_{Λ,0}) V(f)) of positive type; the 2 log Λ term cancels
exactly against Tr(S_Λ V(h)), Connes 1999 §VIII eq. (16) and eqs. (17)–(26), as cited in G.4d) and Connes–Consani 2021 Thm 1
(archimedean, W_∞(g∗g*) ≥ Tr(ϑ(g) S ϑ(g)*) for supp g ⊂ [2^{−1/2}, 2^{1/2}] with the vanishing conditions of their Thm 1; Sonin-space
trace). Per G.4d, Connes–Consani 2021 confirms the deficit sign only: it is the Λ = 1 case, where 2 log Λ = 0, so it cannot test the
cancellation. What stays open is unchanged: the number-field, all-cutoff statement lim_Λ Δ_Λ = Weil (δ_Λ → 0) is G.4d itself
(GAP, RH-equivalent), and the ℚ-case of the μ_ρ residual formula is CONJECTURE (G.4d).

---

## S.4 Fakes (R5) for every §S step

| step | F_θ (θ<½) | DH | H₋ₜ | DMV (Beurling) / P_X | G_θ′ (must pass) | Epstein x²+5y² |
|---|---|---|---|---|---|---|
| S(a)1–2 MF = Cauchy, mass 1 | MF ≥ 0 (two Cauchy kernels, G.3c) | no pole ⇒ no MF | same MF | same MF | passes (twist identity) | MF ≥ 0, same shape |
| S(b)3 items 1–6 (index bookkeeping) | holds | holds; κ(∞) can be ∞ | holds | holds | κ ≡ 0 iff ζ's is | holds |
| S(b)4 / S(c)4 (Fl not PSD; no slack) | holds | Fl = Q (no MF), Q indefinite | holds | holds (agrees with ζ on bounded windows, R1) | holds | holds |
| decisive break | G.1 temperedness / G.4 polarization | G.1/G.4 | G.2/G.4 | G.4a | never | G.1d/G.4b |

Reading: every §S statement is fake-blind. §S locates no RH content; it sharpens *how* the open step G.4d must look
(no MF-relative margin, saturated on the soft mode — PROVED; "must be an exact identity" — CONJECTURE, see S(c)4). Consistent with ARP G.4e.

---

## S.5 Prior-art check

- Weil 1952 (criterion); Yoshida 1992 (window positivity for small support; positivity on finite-codim K_N(a); degeneracy at the first
  failure) — all S(b)3 window facts are his or elementary.
- Bombieri 2000 (index = half the off-line count; minimum attained on the unit ball with support in [−t, t]) — S(b)3.7.
- Burnol, C. R. Acad. Sci. Paris Sér. I 331 (2000) 423–428, "Sur les formules explicites I : analyse invariante" (conductor operator;
  Weil positivity verified under a support condition) — consistent with S(b)3.2.
- Connes 1999; Meyer 2005; Connes–Consani 2021 — absorption sign, cokernel, Sonin-space positivity (S(b)1, S(c)5).
- Suzuki, arXiv 2206.03682 (screw function of ζ; Thms 1.2–1.5, 6.1; Kreĭn–Langer Satz 5.9) — the cleanest continuous-kernel version
  of the Yoshida–Bombieri window analysis; RH ⇔ Ψ ≥ 0 (Thm 1.7).
- The MF/Fl split with Cauchy density and S(b)4/S(c)4: the facts that the Weil form's infimum on the unit ball tends to 0 (under RH)
  and that MF is a Cauchy kernel are elementary consequences of the above. **No novelty is claimed** for any §S result.

## S.6 Dead-ends log (this file; cite §0.2, not re-derived)

1. "Bochner for Fl" — false (S(c)1). 2. "Bochner for MF + Fl" — Weil restatement (§0.2 row 2). 3. "MF dominates Fl with margin" —
   ruled out unconditionally for the MF-relative margin only (S(c)4; other norms not excluded). 4. Analytic continuation of the MF character to off-line exponents — indefinite (S(c)3,
   `GppAbelContinuedKernel`). 5. Window-index certificates → global — Kreĭn–Langer extension freedom (S97 §B) and R1 (DMV).
6. PNT error bounds for Fl — exponentially too weak at half density (ARP G.4d(i)). 7. Positive harmonic-measure smearing of off-line
   atoms — no bound on Q (S(c)3 remark). Also §0.2: log-time OS false; positive cone fake-blind (S98); Selberg; Connes uniformity;
   η_ζ = 0; prime weights (H100a/H96a); S97 continuous Kreĭn reflection.

## S.7 Discrepancies found in ARP.md (origin/main 36aa6be; not edited — this commit touches only this file)

- **D1.** G.3b line 165: "E₁ … This is the Chebyshev discrepancy ψ(x) − x". E₁ (primes only) is the θ-discrepancy θ(x) − x;
  the ψ-discrepancy is E_Λ. Suggest: "θ(x) − x; superseded by E_Λ (G.3d)".
- **D2.** G.4d line 210: said "E₁-fluctuation" at 3889125; **already corrected to E_Λ at 36aa6be.** §S uses E_Λ throughout.
- **D3.** G.5 (T1) still says the k = 2 rank-one "cancels against θ's −√x" — the withdrawn term-by-term wording; should cite
  ψ = Σ_k θ(x^{1/k}).
- **D4.** G.3d "(normalisation match to Astra's ¼ not checked)" — reconciled for the Toeplitz half in S(a)3.
- **D5.** §0.3 writes `GppWeilSupportLadder.primeSide_eq_*`; in GPPVerify (local clone 8a49a74) the namespace is `GppWeilLadder`.
- **D6.** G.4c cites "Connes 1999 Thm" for the S-local trace formula: that is Connes 1999 Theorem 4 (§VII); the local version is
  Theorem 3 (§V).

## S.8 Final GAP (identical to ARP G.4d; RH-equivalent; no progress claimed)

**GAP (G.4d).** Fl_L(f) ≥ −MF_L(f) for all f ∈ C_c^∞([−L, L]), uniformly in L, with Fl_L the **E_Λ**-fluctuation.
Equivalent to RH (Weil 1952; Bombieri 2000). §S adds only: (i) MF_L = ∫ĝ dt/(2π(t²+¼)) exactly; (ii) the inequality has no MF-relative
margin and is saturated on dilated nodeless bumps (S(c)4, PROVED; margins in norms degenerating on those bumps are not excluded;
"exact identity needed" is CONJECTURE); (iii) the window negative index is monotone, finite, eventually
global and ≤ N_{>½}.
Secondary open items: GAP (verification/transfer) κ_Q(L) = N_{>½} on support windows (S(b)3.7); GAP (citation) Kreĭn–Langer
𝒢_κ count for ζ (S(b)3.8); GAP (normalization, minor) Astra's Hankel −⅛ (S(a)3); cokernel sign SETTLED (deficit on coker E = ker E*; G.4d on main, 7328e39; S(c)5).

## S.9 References (verified this session unless marked)

- Bombieri, E., "Remarks on Weil's quadratic functional in the theory of prime numbers, I", Atti Accad. Naz. Lincei Rend. Lincei
  (9) Mat. Appl. 11 (2000), no. 3, 183–233. MR1841692. (Abstract and introduction verified; full text not accessed.)
- Burnol, J.-F., C. R. Acad. Sci. Paris Sér. I 331 (2000), no. 6, 423–428; arXiv math/0101068.
- Connes, A., "Trace formula in noncommutative geometry and the zeros of the Riemann zeta function", Selecta Math. (N.S.) 5 (1999)
  29–106; arXiv math/9811068.
- Connes, A., Consani, C., "Weil positivity and trace formula, the archimedean place", Selecta Math. (N.S.) 27 (2021), no. 4,
  Paper 77; arXiv 2006.13771.
- Iohvidov, I. S., Kreĭn, M. G., Langer, H., Introduction to the Spectral Theory of Operators in Spaces with an Indefinite Metric,
  Akademie-Verlag, Berlin, 1982 (not re-verified this session).
- Kreĭn, M. G., Langer, H., Math. Nachr. 77 (1977) 187–236 (Satz 5.9 as cited by Suzuki [14]).
- Meyer, R., Duke Math. J. 127 (2005), no. 3, 519–595; arXiv math/0311468.
- Suzuki, M., "Aspects of the screw function corresponding to the Riemann zeta-function", arXiv 2206.03682.
- Weil, A., "Sur les « formules explicites » de la théorie des nombres premiers", Comm. Sém. Math. Univ. Lund (1952).
- Yoshida, H., "On Hermitian forms attached to zeta functions", Adv. Stud. Pure Math. 21 (1992) 281–325 (Lemmas 2, 3, Thm 2 as
  reported by Suzuki §4.2, §5).
- Project: `derivation/ARP.md` (36aa6be); `team/spectrum-measure/KREIN_LANGER_S97.md`; GPP-bridge
  `work/2026-10-27_astra_prime_pnt_os_ward_identity.md`; GPPVerify `WeilPositivityCriterion.lean`, `WeilSupportLadder.lean`,
  `AbelContinuedKernelIndefinite.lean`.
