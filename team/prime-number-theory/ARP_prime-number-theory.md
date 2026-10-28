# §N Goddess of Prime Number Theory: the Weil form as a limit of finite Euler-product forms

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Seat: Goddess of Prime Number Theory, 2026-10-28 (S103), branch `team/prime-number-theory`, built on GPPDiscovery5 main 8079b30
(ARP §G incl. 36aa6be/66cbedf/9184613, 0244ce4/f097a52/7328e39; §S 9a10243, §P, §Q, §C merged; aligned in N.9). Read-only consistency with §P at `team/principal-series` 1c46f59.

**RH is OPEN. Nothing in this section proves or implies RH.** The results are structural. On the support ladder P_X *is* ζ, and past
the ladder every finite Euler product fails Weil positivity within O(1) of the ladder edge. So no finite set of primes closes G.4d.

**Labels** (ARP §0): PROVED (citation or elementary proof given here), PROVED-LEAN (GPPVerify name), CONJECTURE, GAP,
NUMERICAL (orientation or constant evaluation only, never a proof step), CERTIFIED (ball arithmetic, earlier runs).
`derivation/LEAN_INVENTORY.md` does not exist. ARP §0.3 ("Inventory") says it serves as that inventory, so Lean names are taken
from there and checked against the GPPVerify source (`GppVerify/RiemannHypothesis/WeilSupportLadder.lean`, `namespace GppWeilLadder`).
Every step has an R-line. R1 is uniform in L. R2 is Euler product and FE jointly. R3 is no reliance on prime-weight signs.
R4 is that the positive space contains the soft mode. R5 is where F_θ, DH, H₋ₜ, DMV, Epstein break, with G_θ′ as the must-pass control.

---

## N.0 Conventions (those of ARP §G and §S S.0, unchanged)

- Window: f ∈ C_c^∞([−L, L]) real (L is the half-width), g = f ∗ f̃, supp g ⊂ [−2L, 2L], f̂(t) = ∫ f e^{itu} du, ĝ = |f̂|² on ℝ.
  Write g_s(u) := g(u) + g(−u) on (0, ∞) and F(s) := ∫ f(u) e^{su} du, so f̂(∓i/2) = F(±½).
- ζ (§S S.0, Weil 1952; Suzuki 2206.03682 eq. (5.15)):
  Q_L(f) = ĝ(i/2) + ĝ(−i/2) − Σ_{n≥2} Λ(n) n^{−1/2} g_s(log n) + W_∞(g), W_∞(g) = −(log π) g(0) + (1/2π)∫ Re ψ(¼ + it/2) ĝ(t) dt.
- Split (§S S.0, ARP G.3a/G.3b/G.3d): Q_L = MF_L + Fl_L,
  MF_L(f) = ĝ(i/2) + ĝ(−i/2) − ∫_0^∞ e^{u/2} g_s(u) du = (1/2π)∫ |f̂|²/(t² + ¼) dt (§S S(a)1–2, = G.3a at v = 0),
  Fl_L(f) = W_∞(g) − ∫_{(0,∞)} g_s dE_Λ, E_Λ = μ_Λ − e^{u/2} du, μ_Λ = Σ_n Λ(n) n^{−1/2} δ_{log n}.

**Convention N.0 (definition of the P_X "system").** P_X(s) = ∏_{p≤X} (1 − p^{−s})^{−1} has no functional equation, no pole at
s = 1 and no Γ-factor, so it has no Weil explicit formula of its own. As in `finite_windows_dmv.md` §1 we **define**
  Q_L^{P_X}(f) := ĝ(i/2) + ĝ(−i/2) − ∫ g_s dμ_Λ^{P_X} + W_∞(g),  μ_Λ^{P_X} := Σ_{p≤X, k≥1} (log p) p^{−k/2} δ_{k log p},
i.e. **ζ's pole terms and ζ's archimedean functional W_∞ are kept (borrowed)**, and only the prime-power measure is truncated
(−P_X′/P_X(s) = Σ_{p≤X,k} (log p) p^{−ks}). In Weil/Connes language, Q^{P_X} is the pole terms plus Σ_{v ∈ S_X} W_v with
S_X = {∞} ∪ {p ≤ X}, the semi-local sum of the local distributions of Connes 1999 Thm 4 / ARP G.4c, without the cutoff term.
The only fully RH-relevant object is X = ∞. Q^{P_X} is a comparison form, not "the Weil form of P_X".

## N.1 (a) The mean field is X-independent

**N.1 (PROVED, by definition plus §S S(a)1).** Put E_Λ^{P_X} := μ_Λ^{P_X} − e^{u/2} du (the *same* subtracted mean) and
Fl_L^{P_X} := W_∞(g) − ∫ g_s dE_Λ^{P_X}. Then
  **Q_L^{P_X} = MF_L + Fl_L^{P_X}, with MF_L exactly ARP's/§S's MF_L, independent of X.**
*Proof.* Add and subtract ∫_0^∞ e^{u/2} g_s du in Convention N.0. MF_L involves only the pole terms and the PNT density e^{u/2}du,
and neither depends on X. ∎
*Remark (convention, stated honestly).* e^{u/2}du is **not** the PNT mean of μ_Λ^{P_X}: ψ_{P_X}(x) = Σ_{p≤X} (log p)⌊log x/log p⌋ =
O_X(log x). Keeping MF fixed therefore puts a non-decaying mean −e^{u/2}1_{u>log X} du into E_Λ^{P_X} (see N.4). This is ARP
G.3c's "P_X: same pole, same PNT mean, same MF", made precise. The MF is the same *by the choice of decomposition*.
R1 identity for every L | R2 not used | R3 not used | R4 MF contains the pole state e^{−u/2} (G.3a) for every X | R5 fake-blind like G.3c.

## N.2 (b) The support ladder, sharp

Let p⁺(X) be the least prime > X.

**N.2 Lemma (PROVED, elementary; the tail-vanishing step is PROVED-LEAN `GppWeilLadder.primeSide_eq_truncation`,
`GppWeilLadder.primeSide_term_eq_zero`, with support doubling `GppWeilLadder.convolution_hasSupportIn`).**
(i) For every X ≥ 2 and every L with **2L ≤ log p⁺(X)** (sharp form), Q_L^{P_X}(f) = Q_L^ζ(f) and Fl_L^{P_X}(f) = Fl_L^ζ(f) for all
f ∈ C_c^∞([−L, L]). The simple sufficient form is **L ≤ (log X)/2**.
(ii) The sharp form is optimal. If 2L > log p⁺(X), there is f ∈ C_c^∞([−L, L]) with Q_L^{P_X}(f) ≠ Q_L^ζ(f).
*Proof.* Q^ζ − Q^{P_X} = Fl^ζ − Fl^{P_X} = −D_X(f) with
  **D_X(f) := Σ_{p > X, k ≥ 1} (log p) p^{−k/2} g_s(k log p) = ∫ g_s dμ_miss,** μ_miss := μ_Λ − μ_Λ^{P_X}.
supp g ⊂ [−2L, 2L] (Lean: `convolution_hasSupportIn` with HasSupportIn at 2L), and g(±2L) = 0 because the integrand of
g(2L) = ∫ f(u) f(u − 2L) du is supported at the single point u = L. So only atoms with k log p < 2L contribute. In Lean this is
`primeSide_term_eq_zero` for log n > 2L, and the boundary atom log n = 2L contributes g(2L) = 0. A prime power p^k with p > X
satisfies p^k ≥ p ≥ p⁺(X), so k log p ≥ log p⁺(X) ≥ 2L and every term of D_X vanishes. This is (i). Since p⁺(X) > X, the simple
form follows. For (ii), put a := log p⁺(X) < 2L and f := φ(· − a/2) + φ(· + a/2), with φ ∈ C_c^∞, φ ≥ 0, φ ≢ 0, supported in
(−η, η), η < min(L − a/2, (a − log X)/4, a/4). Then g ≥ 0 is supported in {|u| < 2η} ∪ {|u ∓ a| < 2η}, and g(a) = ‖φ‖² > 0.
The only missing atom with |log n − a| < 2η is p⁺(X), and missing atoms satisfy log n ≥ a, so D_X(f) ≥ (log p⁺) p⁺^{−1/2} g(a) > 0. ∎
Boundary bookkeeping: prime powers p^k with p ≤ X are *all* in P_X (25 = 5², 27 = 3³ ∈ P_23). The criterion is therefore "every
prime power n < e^{2L} has its prime ≤ X", which is "every prime < e^{2L} is ≤ X", which is p⁺(X) ≥ e^{2L}. Example: X = 23,
p⁺ = 29, sharp ladder L ≤ ½ log 29 = 1.68365, simple ladder L ≤ ½ log 23 = 1.56775.

**N.2′ What the FE contributes (PROVED, elementary, given N.0).** P_X has no FE. Its form carries ζ's FE data only through the
borrowed pole + W_∞, which are identical for ζ and P_X. On the ladder the difference Q^ζ − Q^{P_X} = −D_X is identically 0. So
**the FE is invisible to the comparison for L ≤ ½ log p⁺(X)**, and off the ladder the whole difference is the missing-prime
functional D_X on (log p⁺(X), 2L). The FE does not enter D_X at all. It enters only through the *match* between the borrowed W_∞ +
pole and the complete prime side. That match is Weil's explicit formula for ζ and is absent for every finite X (ARP R2; DMV row of
G.4e).
**N.2″ No-slack on the ladder (PROVED, from N.2 plus §S S(c)4).** With r^{P_X}(L) := inf_f Q_L^{P_X}/MF_L (MF_L > 0, §S S(a)2),
r^{P_X}(L) = r(L) for L ≤ ½ log p⁺(X). So r^{P_X} is nonincreasing there, and RH ⇔ r^{P_X}(L) ≥ 0 on the full ladder for every X.
(⇒ by Weil. ⇐ because every L lies on the ladder of some X.) Past the ladder see N.3.4. The other direction fails: r^{P_X} → −∞.
R1 the ladder is per window. Uniformity would need X → ∞, i.e. ζ itself | R2 the FE is invisible on the ladder: this *is* the R2
statement for P_X | R3 not used | R4 trivial (the forms coincide) | R5 DMV/P_X agree with ζ here (the content of
`finite_windows_dmv.md`). F_θ, DH, G_θ′ and Epstein differ from ζ already at n = 2, so the ladder is not about them. H₋ₜ is not a
prime-side system.

## N.3 (c) The negative index past the ladder

### N.3.1 Sign bookkeeping: removing primes adds an *indefinite* form (PROVED, elementary)
Q^{P_X} = Q^ζ + D_X. μ_miss is a positive measure, but **f ↦ D_X(f) = ∫ g_s dμ_miss is not a positive form**. A single atom gives
f ↦ g(a) = ∫ f(u) f(u − a) du, which is indefinite for a > 0. With f = φ(· − a/2) + ε φ(· + a/2) and φ of support < a/2,
g(a) = ε‖φ‖², with either sign. So the heuristic "removing the subtracted term −ΣΛ n^{−1/2} g_s can only make the form more
positive" is **false as a statement about forms**. It is true only along test functions with g ≥ 0 at the missing atoms, and the
soft mode does not stay in that class.
R3 this is R3 on the prime side: the weight sign of the removed atoms carries no sign for the form | R5 identical for any truncation.

### N.3.2 Multiplier form of Q^{P_X} (PROVED; explicit formula convention §S S.0, Plancherel)
For finite X, μ_Λ^{P_X} is a finite measure (mass Σ_{p≤X} log p/(√p − 1)). For f ∈ C_c^∞([−L, L]) and every L,
  **Q_L^{P_X}(f) = 2F(½)F(−½) + (1/2π)∫_ℝ |f̂(t)|² m_X(t) dt,**
  m_X(t) := Ω(t) − 2 Re D_X(t), Ω(t) := Re ψ(¼ + it/2) − log π, D_X(t) := Σ_{p≤X} (log p) z_p/(1 − z_p), z_p = p^{−1/2−it}.
*Proof.* ĝ(±i/2) = f̂(i/2) f̂(−i/2) = F(−½)F(½) for real f. Also g(a) = (1/2π)∫ |f̂|² e^{−ita} dt, so ∫ g_s dμ_Λ^{P_X} =
(1/2π)∫ |f̂|² · 2Re Σ_{p≤X,k} (log p) p^{−k/2 − ikt} dt (a finite sum of absolutely convergent geometric series), and W_∞ is as in
N.0. ∎
Properties (PROVED, elementary): m_X is real-analytic and |2 Re D_X(t)| ≤ A_X := 2Σ_{p≤X} log p/(√p − 1). Ω(t) = log(|t|/2π) + O(1/t²)
(Stirling), so m_X(t) ≤ A_X + Ω⁺(t). Also m_X(0) = Ω(0) − A_X < 0, with Ω(0) = −γ − π/2 − 3 log 2 − log π = −5.3722.
A_X ~ 4√X by PNT and A_X = O(√X) by Chebyshev's bound (NUMERICAL values: A_2 = 3.347, A_23 = 20.576, A_{10^7} = 12659.2 vs
4√X = 12649.1; `sN/omega0_check.out`, `sN/odd_pole_bound.out`).
The pole term is a **hyperbolic rank-2 form**: 2F(½)F(−½) = ½(F(½) + F(−½))² − ½(F(½) − F(−½))², of signature (1, 1). For ζ the full
prime tail cancels its growing part (Ward identity G.3a). For P_X nothing cancels it, because the prime multiplier is bounded.
Equivalently (§P P.6), the Euler operator L_{S_X} is a finite product, bounded on the critical-line space α = ½.

### N.3.3 Proposition (PROVED, elementary, unconditional): every P_X has a negative direction within O(1) of the ladder
Let b ∈ C_c^∞([−w, w]) be even and real with b ≢ 0. Put B := ∫ b(u) e^{u/2} du > 0 (for b ≥ 0), ω_b := (1/2π)∫ |b̂|² Ω⁺ dt, and for
c > w let f_c(u) := b(u − c) − b(u + c) ∈ C_c^∞([−(c + w), c + w]). If
  **8 sinh²(c/2) B² > 4(A_X ‖b‖² + ω_b),**
then Q_L^{P_X}(f_c) < 0 for L = c + w. Hence ind_−(Q_L^{P_X}) ≥ 1 for all L ≥ L₁(X) := w + 2 arcsinh √((A_X‖b‖² + ω_b)/(2B²)),
and **L₁(X) = ½ log X + O(1)** (because A_X = O(√X)).
*Proof.* F(½) = (e^{c/2} − e^{−c/2}) B = −F(−½), so the pole term is −8 sinh²(c/2) B². f̂_c(t) = 2i sin(ct) b̂(t), so the
multiplier term is (1/2π)∫ 4 sin²(ct) |b̂|² m_X ≤ 4(1/2π)∫ |b̂|² (A_X + Ω⁺) = 4(A_X‖b‖² + ω_b) (Plancherel). Monotonicity in L is
test-space inclusion. ∎
- NUMERICAL evaluation of the constant (`sN/odd_pole_bound.py`, b = cos²(πu/2w), w = 1: B = 1.0164, ‖b‖² = 0.75,
  ω_b = 9.06e−5): L₁(23) = 4.461, L₁(10³) = 6.245, L₁(10⁵) = 8.522, L₁(10⁷) = 10.819. The excess L₁(X) − ½ log X is
  2.893, 2.791, 2.766, 2.760, tending to w + log(8‖b‖²/B²) = 2.760 with this b. The bound is crude: the certified failure of P_23 is
  at 1.7105 (N.3.5).
- Consistency check (firing on the ladder would contradict RH): the excess L₁(X) − ½ log X is ≥ 2.76 at every tested X
  (NUMERICAL). For every b, the asymptotic excess w + log(8‖b‖²/B²) ≥ w + log(4/sinh w) > log 8 > 0, by Cauchy–Schwarz
  B² ≤ ‖b‖²∫_{−w}^{w} e^{u} du = 2‖b‖² sinh w (PROVED, elementary).
- **Corollary (PROVED).** For every X, **G.4d fails for P_X**: Fl_L^{P_X}(f_c) < −MF_L(f_c) for L ≥ L₁(X), since
  Q^{P_X} = MF + Fl^{P_X} by N.1.
R1 explicit and uniform in X: one formula L₁(X) for all X | R2 uses neither the Euler product nor the FE, only that P_X's prime
multiplier is bounded (§P P.6: L_{S_X} bounded at α = ½) | R3 uses |weights| (A_X), not their signs | R4 the witness is the **odd**
partner of the pole state: the negative direction of the hyperbolic pole form, F(½) = −F(−½). It is not the even nodeless soft mode.
ζ's full prime mean turns pole − PNT mean into MF ≥ 0 (G.3a), and truncation undoes exactly that | R5 the hypothesis "bounded prime
multiplier plus ζ's pole" holds for P_X only. F_θ, G_θ′, DMV and Epstein have unbounded prime sides, DH has no pole (entire), and
H₋ₜ has no prime side, so the Proposition is silent on all of them. In particular **G_θ′ is not falsely failed**.

### N.3.4 Growth in L for fixed X, and the no-slack ratio (PROVED unless marked)
- (a) ind_−(Q_L^{P_X}) → ∞ as L → ∞, for each fixed X (PROVED). Write M(f) := (1/2π)∫ |f̂|² m_X. Choose δ > 0 with
  m_X ≤ m_X(0)/2 =: −c₀ < 0 on |t| ≤ δ. Take φ ∈ C_c^∞([−½, ½]), ‖φ‖ = 1, and k translates f_j = ℓ^{−1/2} φ((· − c_j)/ℓ) with disjoint
  supports in [−L, L] (kℓ ≤ 2L). For F = Σ a_j f_j, Cauchy–Schwarz gives |F̂|² ≤ k Σ a_j² |f̂_j|², so
  M(F) ≤ −c₀‖a‖² + k‖a‖² T(ℓ), T(ℓ) := (1/2π)∫_{|s|>δℓ} |φ̂(s)|² (c₀ + A_X + Ω⁺(s/ℓ)) ds → 0 as ℓ → ∞. So ind_−(M) ≥ k once
  kT(ℓ) < c₀ and kℓ ≤ 2L. The pole form has ind_+ = 1, so ind_−(Q^{P_X}) ≥ ind_−(M) − 1 ≥ k − 1. ∎
- (b) Rate (CONJECTURE, Szegő-type heuristic). ind_−(Q_L^{P_X}) ~ (L/π)·|{t : m_X(t) < 0}| as L → ∞, the truncated-convolution
  eigenvalue count. It is PROVED for indicator symbols by Landau–Widom, J. Math. Anal. Appl. 77 (1980) 469–481. A version for the
  symbol m_X (unbounded, log-growing) is not checked: theorem number unverified. |{m_X < 0}| is finite (Ω → ∞) and grows with X.
- (c) **No-slack past the ladder (PROVED).** MF_L ≤ 4‖·‖² (§S S(a)2) and ‖f_c‖² = 2‖b‖², so whenever Q^{P_X}(f_c) < 0,
  r^{P_X}(c + w) ≤ [−8 sinh²(c/2) B² + 4(A_X‖b‖² + ω_b)] / (8‖b‖²) → −∞, like −e^{L} B² e^{−w}/(4‖b‖²).
  Contrast with ζ (§S S(c)4): r(L) is nonincreasing with lim ≤ 0, and RH ⇔ lim r = 0. **ζ's ratio stays at the edge only with the
  complete prime tail. Truncating at any X sends it to −∞ exponentially, within O(1) of ½ log X.**

### N.3.5 Data, and the difference between plain truncation and glued Beurling systems
- **Plain truncation, CERTIFIED.** At L = 1.7105 (= 3421/2000), K = 200 (the px pipeline's cos basis), the window sees n < e^{3.421} = 30.6. With
  GPRIMES = {2, …, 23} (the `del29` run, `precision/runs/px/px_results.jsonl` lines 13–16, local run log; summarized in `px_positivity.md` §A), this is exactly Q^{P_23} (25 and 27 are
  included; 29 is the only missing atom). Ball-arithmetic LDLᵀ gives a certified negative pivot [−8.3402e−14 ± 4.94e−19] at index 5, so by Sylvester
  **ind_−(Q_{1.7105}^{P_23}) ≥ 1, CERTIFIED**, at 0.027 past the sharp ladder edge 1.68365.
- **Glued systems (different objects).** The `px_positivity.md` §A runs (CERTIFIED PD at L = 1.6035, CERTIFIED non-PSD at 1.6043,
  K = 160) use P = {p ≤ 23} ∪ {24.7, 28.3}. That is a Beurling system that **adds** a non-rational atom (adding −2(log q)q^{−1/2} g(log q))
  *and* removes 29. Its first deviation from ζ is at q₁ = 24.7 (L₀ = 1.60337), not at p⁺(23). Those data are evidence for the glued
  DMV-type fakes of `finite_windows_dmv.md`, not for plain truncation. Both kinds fail just past their first deviation, with
  opposite-sign perturbations, as N.3.1 predicts (indefinite either way).
- **CONJECTURE N.C1 (onset at the edge).** For every X, the onset L*(X) := sup{L : Q_{L′}^{P_X} ≥ 0 ∀ L′ ≤ L} satisfies
  L*(X) − ½ log p⁺(X) → 0. Evidence: the glued onset within 8·10⁻⁴ of its edge, the plain P_23 failure within 0.027 (one X only),
  and the near-null spectrum of Q^ζ at the edge (λ_min ~ 10⁻¹²¹ at K = 160). PROVED bounds: L*(X) ≤ L₁(X) = ½ log X + O(1)
  (N.3.3). Under RH, L*(X) ≥ ½ log p⁺(X) (N.2 + Weil). Unconditionally, L*(X) ≥ min(½ log p⁺(X), sup of windows where Q^ζ is
  proved PSD), with the latter ≥ (log 2)/2 (Yoshida 1992; ARP G.2b).
- **Answer to (c).** Yes. **ind_−(Q_L^{P_X}) ≥ 1 for some L, for every X: PROVED** (N.3.3, at L ≤ ½ log X + O(1); explicit
  L₁(23) ≤ 4.461), and CERTIFIED at L = 1.7105 for X = 23. Moreover ind_− → ∞ as L → ∞ (PROVED, rate CONJECTURE).
- **DMV (Theorem 1 of Diamond–Montgomery–Vorhauer, "Beurling primes with large oscillation", Math. Ann. 334 (2006) 1–36,
  doi:10.1007/s00208-005-0638-2; verified this session).** For 1/2 < θ < 1 and a > (4/e)(1 − θ) there are discrete Beurling primes
  with N(x) = κx + O(x^θ), ζ_B analytic on σ ≥ θ apart from a simple pole at 1, infinitely many zeros on σ = 1 − a/log t, and
  ψ_B(x) = x + O(x exp(−2√(a log x))) and ψ_B(x) − x = Ω_±(x exp(−2√(a log x))) (Thm 1 (i)–(v); as transcribed in
  `finite_windows_dmv.md` §2). Glued at X (`finite_windows_dmv.md` §2–3, PROVED), it agrees with ζ on the ladder. Unlike
  plain P_X it **keeps a full prime tail with PNT mean e^{u/2}du**, so its multiplier is unbounded and N.3.3 does **not** apply.
  **GAP (DMV past the ladder):** is Q_L^{P_DMV} non-PSD for some L? No explicit formula ties the borrowed W_∞ to DMV's zeros (no FE),
  so DMV's off-line zeros give no handle (`finite_windows_dmv.md` §4.1). CONJECTURE: yes, by the glued evidence above.
R1 the onset results are per X, and the lower bound needs ζ positivity (RH) | R2 P_X has the Euler product, no FE: N.3.3 shows a
bounded Euler operator cannot hold the pole | R3 none of N.3 uses weight signs | R4 the failing directions are the odd pole partner
(N.3.3) and, in the certified even-sector run, a negative LDLᵀ pivot in the cos basis (index 5) | R5 as in N.3.3. DMV is the open case.

## N.4 (d) The ψ-fluctuation for P_X (S103 rewrite)

**N.4a (PROVED, classical, cited as ARP G.3d with Golden Eye's repair, 36aa6be).** ψ(x) = Σ_{k≥1} θ(x^{1/k}). So the prime-power layers
μ_k = Σ_p (log p) p^{−k/2} δ_{k log p} sum to μ_Λ, the half-density of ψ, whose PNT mean is e^{u/2}du. Taking the fluctuation
against ψ leaves **no separate k = 1, 2 mean-field layers**: the pole is spent once, on the Ward identity G.3a (Astra, Grok-replayed
PROVED, MF covariance e^{−|u−v|/2}). The same identity holds verbatim for P_X, with ψ_{P_X}(x) = Σ_k θ_{P_X}(x^{1/k}) and
θ_{P_X}(x) = Σ_{p ≤ min(x,X)} log p. (We do not use the term-by-term θ = Σ μ(k)ψ(x^{1/k}) argument; it is false at k = 4 and 6.)

**N.4b (PROVED, elementary).** With μ_miss = Σ_{p>X,k} (log p) p^{−k/2} δ_{k log p} ≥ 0,
  **E_Λ^{P_X} = E_Λ − μ_miss = E_Λ − e^{u/2}1_{(log X, ∞)}(u) du − E_miss,**  E_miss := μ_miss − e^{u/2}1_{(log X,∞)} du,
  **Fl_L^{P_X} = Fl_L^ζ + ∫ g_s dμ_miss = Fl_L^ζ + ∫_{log X}^∞ e^{u/2} g_s du + ∫ g_s dE_miss.**
Reading. (i) On the sharp ladder every term past Fl^ζ vanishes (N.2). Note that ∫_{log X}^{2L} e^{u/2} g_s cancels against
∫ g_s dE_miss there, since μ_miss has no atoms below p⁺(X). (ii) Off the ladder the middle term is the **unspent part of the Ward
identity**. MF_L + ∫_{log X}^∞ e^{u/2} g_s = pole − ∫_0^{log X} e^{u/2} g_s, i.e. G.3a with its PNT mean cut at u = log X. The
uncancelled pole growth e^{c}B² on the odd pole partner is what N.3.3 exhibits. (iii) E_miss is the (signed) ψ-discrepancy of the
missing primes. Its size is ζ's PNT error on (X, e^{2L}), which is RH-strength to control (G.4d(i)).
R1 identity for all L, X | R2 not used | R3 μ_miss ≥ 0 does not make ∫ g_s dμ_miss ≥ 0 (N.3.1) | R4 the unspent mean acts on the pole
state's hyperbolic partner | R5 fake-blind (G.3d).

**N.4c Normalization of Astra's k = 2 even bulk −¼·11ᵀ against ARP's ½ (PROVED, elementary; explicit computation).**
Astra (`work/2026-10-27_astra_k2_dom_bulk_sewing.md`): b = sin²(πu)1_{[0,1]}, ∫b = ½, G₀ = b ∗ b̃, ∫G₀ = (∫b)² = ¼, bumps
b_j = b(· − a_j), test vector f_v = Σ v_j b_j.
- *Layer density.* The k = 2 layer mean in u = log n is ½·du (Σ_{p≤y} (log p)/p = log y + O(1), u = 2 log p). In ARP's convention its
  contribution to any window is −∫_0^∞ ½ g_s du = −½∫_ℝ g = −½ ĝ(0) = **−½(∫f)²**. That is ARP G.3d's coefficient.
- *Toeplitz block.* For f = f_v: (∫f_v)² = (½Σv_j)², so the bulk is −½·¼(Σv)² = vᵀ(**−⅛·11ᵀ**)v. This is exactly Astra's Q^{(2)}
  bulk −⅛ (b_i ∗ b̃_j = G₀(· − (a_i − a_j)) reproduces her Q^{(2)}_{ij} formula).
- *Hankel block and even sector.* Astra's C^{(2)}_{ij} = −∫ G₀(u − a_i − a_j − 1) dμ₂ is the Weil pairing of b_i with the reflected
  bump b_j(−·). Indeed b_i ∗ (b_j(−·))~ (u) = (b ∗ b)(u − a_i − a_j) = G₀(u − a_i − a_j − 1), using b(u) = b(1 − u); the g(−u) half
  vanishes since a_i + a_j + 1 > 1. Put F_v := f_v + f_v(−·) (even). The Weil form is invariant under f ↦ f(−·), so
  Q(F_v) = 2vᵀ(Q + C)v, i.e. **vᵀP⁺v = ½Q(F_v)**. With ∫F_v = Σv_j, the k = 2 bulk of ½Q(F_v) is ½·(−½)(Σv)² = vᵀ(**−¼·11ᵀ**)v. ✓
- *Result.* −¼ = 2 × (−½) × (∫b)², i.e. ARP's layer coefficient ½, times Astra's bump normalization (∫b)² = ¼, times the even-sector
  factor 2 from P⁺ = Q + C = ½ Q(F_v) with ∫F_v = 2∫f_v. **ARP's ½ and Astra's ¼ agree.** Golden Eye's guess ("even-sector factor")
  is right for the factor 2, and the remaining factor is (∫b)² = ¼. The odd sector F_v^− = f_v − f_v(−·) has ∫F_v^− = 0, which is why
  Astra's odd k = 2 bulk vanishes. (NUMERICAL: ∫b = 0.5, (∫b)² = 0.25, `sN/omega0_check.out`.)
- Per N.4a this rank-one is a ψ/θ splitting artefact and needs no compensation (T1 withdrawn).
R1 identity | R2 — | R3 — | R4 the vector is 1 (constant, nodeless): the k = 2 bulk sits on the soft direction, consistent with S97/S98 |
R5 F_θ: Astra's two-species bulk −2M(θ)cosh cosh (even), +2M(θ) sinh sinh (odd), with M(0) = ¼ consistent with the above.

## N.5 The ℚ^*-sum in the finite-S picture (consistency with §P P.6)

**N.5 (PROVED as a restatement of §P P.6 plus N.3.3).** §P shows that the only non-isometric step in 𝒮(𝔸) → L²(C_ℚ) is the global
ℚ^*-sum (Meyer 2005 Lemma 5.3, the Euler operator L_S). On the critical-line space it is the multiplier ∏_v L_v(½ + it), with no
operator-norm limit. The P_X form keeps only the places in S_X, so the ℚ^*-sum is replaced by the sum over the S_X-units
{±∏_{p≤X} p^{n_p}} acting semi-locally, and the Euler operator by the **finite product L_{S_X}**, bounded at α = ½ (norm
∏_{p≤X}(1 − p^{−1/2})^{−1}). Its log-derivative gives the bounded prime multiplier 2Re D_X of N.3.2. N.3.3 then shows: **a bounded
Euler operator cannot absorb the Poisson boundary/pole term**. Its hyperbolic negative direction survives once e^{L} ≳ A_X. So
positivity on window L needs the non-unitary ℚ^*-sum to be effectively unbounded at scale e^{L}. Concretely: the primes in
(e^{2L − 2C}, e^{2L}) cannot be dropped, with C = O(1) (C ≈ 2.9 at the tested X).
- Relation to Astra's support bound (ARP G.4d(ii), PROVED by Astra, replayed in §C C.6d; per main 8079b30 the necessity half of T2
  is conditional on H1). Any compensator supported on p ≤ P^α, α < 1, cannot couple off-line mirror sectors. N.3.3 is a form-level
  analogue with α = 1 − O(1/log P). It is **unconditional and zero-free** (no H1): it says Q itself loses positivity when the primes
  in (e^{2L−2C}, e^{2L}) are removed. It does not say how an operator compensator must be supported, so it does not discharge H1 or T2.
R1 per window, uniform constant | R2 this is the Euler-product half only. The FE half (Poisson R-term) is borrowed | R3 — | R4 pole
partner | R5 §P P.6 R-line (P_X bounded, F_θ bounded only for α > 1 + θ, G_θ′ passes).

## N.6 R5 summary for §N

| fake | N.1 MF | N.2 ladder | N.3.3 bounded-multiplier failure | N.4 ψ-fluctuation | where it breaks in §N |
|---|---|---|---|---|---|
| **P_X** (plain) | same MF (convention) | = ζ for 2L ≤ log p⁺(X) | **fails**, L ≤ ½ log X + O(1); CERTIFIED X = 23 at 1.7105 | E_Λ − μ_miss | past the ladder (no lattice / no FE; G.4a) |
| **DMV** (glued) | same MF | = ζ on ladder (`finite_windows_dmv.md`) | silent (full tail) | small dlVP fluctuation | **GAP** past the ladder (CONJECTURE fails) |
| **F_θ** | MF ≥ 0 (G.3c) | differs at n = 2 | silent (unbounded) | 2cosh(θu)-weighted; Astra's odd +2M sinh sinh | not in §N; at G.1/G.4 (polarization; §P) |
| **DH** | no MF (entire) | differs at n = 2 | silent (no pole) | no Λ | not in §N; G.1/G.4 |
| **H₋ₜ** | pole kept | not a prime-side system | silent | — | G.2/G.4 |
| **Epstein** | MF ≥ 0 | differs | silent | two Euler products | G.1d/G.4b |
| **G_θ′** (control) | twisted MF ≥ 0 | differs, but Q_G = 2Q_ζ(f cos) + 2Q_ζ(f sin) | silent: **not falsely failed** | cos-twisted | passes |

## N.7 Bearing on the shared GAP G.4d (Fl_L ≥ −MF_L uniformly in L ⇔ RH)

What the P_X limit **does** contribute (PROVED unless marked):
1. On the ladder, P_X and ζ give the same G.4d instance (N.2), so the finite-S problem with S = S_L is exactly ARP's window problem.
   RH ⇔ every P_X is positive on its full ladder (N.2″).
2. Past the ladder, **G.4d is false for every finite Euler product**, within O(1) of ½ log X (N.3.3), and r^{P_X}(L) → −∞ like
   −e^{L} (N.3.4c). **So the gap cannot be closed by any finite set of primes, and not even by all primes below e^{2L − 2C}.** The
   compensation needs the prime tail up to the window edge, i.e. the unbounded ℚ^*-sum (N.5, §P P.6).
3. Mechanism located: the failure is the **unspent Ward identity** (N.4b). G.3a's pole − PNT-mean cancellation is exact only with the
   complete PNT mean. Truncation leaves the bare hyperbolic pole form. This is fake-blind and carries no information on ζ's zeros.
4. With §S S(c)4: ζ's r(L) is nonincreasing with lim ≤ 0, and every P_X tracks it on the ladder and then plunges. So there is no
   MF-relative margin for ζ (PROVED, §S). The stronger R6 requirement ("the closing step must be an exact identity, a norm or trace
   equal to Q") is **CONJECTURE (well supported)**, not PROVED. §S rules out only an MF-relative margin, not Q_L ≥ c·N_L for some
   other norm degenerating on the dilated bumps. The P_X data are consistent with R6: any inequality whose slack is paid for by a
   finite prime set fails past the ladder.

What it **does not** contribute:
- Nothing about Fl_L^ζ ≥ −MF_L. On the ladder the statement for P_X *is* the statement for ζ, and past it P_X is a different form.
- Nothing about the FE (R2). P_X borrows W_∞ and the pole, so the FE's role (the match of W_∞ + pole with the full prime side) is
  untested by P_X. DMV (full tail, no FE) is the system that would test it, and it is a GAP (N.3.5).
- No effective L* for ζ, and no new route. **The sharp GAP is unchanged: G.4d, equivalent to RH.** §N's sharpening is that the
  compensator for Fl_L must use the primes in the last O(1) of log-scale below e^{2L}, together with the FE match. A mechanism that
  works for a finite S, or that loses a constant, is ruled out (the first is PROVED here, the second is PROVED relative to MF by §S
  S(c)4 and CONJECTURE (R6) in general).

**GAPs stated sharply.**
- **G.4d** (shared, ⇔ RH): Fl_L ≥ −MF_L on C_c^∞([−L, L]) for all L, with Fl_L the E_Λ-fluctuation.
- **GAP N.G1 (DMV past the ladder):** does the glued DMV form Q^{P_DMV}_L (borrowed W_∞, full prime tail, no FE) acquire negative
  index for some L > ½ log X? Not implied by DMV's off-line zeros (no explicit formula).
- **CONJECTURE N.C1:** onset L*(X) = ½ log p⁺(X) + o(1). **CONJECTURE N.C2:** ind_−(Q_L^{P_X}) ~ (L/π)|{m_X < 0}|.

## N.8 Notes for the synthesis (conflicts / wording)
- ARP G.3c says "P_X, DMV: same pole, same PNT mean". For plain P_X the PNT mean is the same **only by convention** (N.1 remark).
  The true mean of μ_Λ^{P_X} is O(log x), and the difference is the unspent Ward term of N.4b. Suggested wording fix, not edited here.
- ARP G.4e row "DMV / P_X ... breaks at G.4a" is consistent. §N adds that plain P_X breaks *provably* past the ladder (N.3.3), while
  glued DMV is a GAP.
- §P (1c46f59, line 153, item 4 "Nothing local bounds the GAP") cites `GppWeilSupportLadder.primeSide_eq_truncation`. The source
  namespace is `GppWeilLadder`, which ARP §0.3 on main already uses. This is a one-word fix for §P's owner. Not edited here (read-only).
- `LEAN_INVENTORY.md` does not exist. ARP §0.3 is used as the inventory.
- ARP §0.1 R1 says DMV/P_X "agree with ζ on every bounded window". The precise statement is per X: 2L ≤ log p⁺(X) (N.2), with
  every bounded window covered by some X. That is consistent, and only a wording sharpening.

## N.9 Alignment with main 8079b30 (S103c steering)

- **(1) Settled sign of the shared GAP (ARP G.4d).** The positive object is a trace over coker E = ker E* (the deficit), and the GAP
  is lim_Λ Δ_Λ = Weil, i.e. δ_Λ → 0 (RH-equivalent). "Additive positivity transfers through E" is WRONG-SIGN and withdrawn. The
  2 log Λ cancellation is Connes 1999 §VIII eq. (16) and the proof of Thm 5, eqs. (17)–(26) (cited from ARP G.4d / Queen
  `arp/sign_check.md`, not re-read). §N works entirely on the Weil-form side (Q_L = MF_L + Fl_L) and uses no positivity "through E".
  For P_X there is no E (no lattice ℚ^*, Queen §Q OS1 fails), so the deficit picture does not apply, and N.3.3 is not a statement about
  δ_Λ. Nothing in §N conflicts with the settled sign.
- **(2) μ_ρ.** §N does **not** use the residual formula lim Δ_Λ = Σ_ρ m_ρ ∫ĝ dμ_ρ. That formula is PROVED for function fields
  (Connes 1999 §VIII Lemma 3) and is CONJECTURE for ℚ (ARP G.4d). P_X has no zeros attached to the borrowed W_∞ (no explicit formula),
  so μ_ρ is not defined for it. N.3.3 is zero-free.
- **(3) Queen's OS axioms (§Q).** OS0–OS4 are PROVED except OS2b = G.4d, THE GAP. §Q's table row for P_X reads "OS2b ✗, first
  failure OS1". **N.3.3 supplies an explicit, unconditional proof of the OS2b failure for plain P_X** (L ≤ ½ log X + O(1); CERTIFIED
  at X = 23, L = 1.7105). It is consistent with §Q and sharpens it. On OS4: §Q says OS4 (pole/boundary) is fully used up by MF. For P_X
  that holds *by Convention N.0*, which keeps ζ's PNT mean in MF. N.4b shows the bookkeeping cost: Fl^{P_X} carries the unspent Ward
  term ∫_{log X}^∞ e^{u/2} g_s. Equivalently, with P_X's own (truncated) mean, OS4 would not be used up, and N.3.3's witness is exactly
  the unconsumed hyperbolic pole direction. §Q's "P_X satisfies OS2′ (MF ≥ 0)" (§Q line 206) is N.1. Wording: §Q R1 cites the simple
  ladder L ≤ (log X)/2. The sharp ladder is 2L ≤ log p⁺(X) (N.2). This is consistent, not a conflict.
- **(4) §C C.2d (KMS/primon positivity) passes P_X. Consistent, and §N makes it quantitative (PROVED).** §C C.2d proves that BC/primon OS
  positivity at β = 1 is equivalent to w_n ≥ 0 (Widder; Berg–Christensen–Ressel), and lists P_X among the fakes that pass (§C, R line
  for C.2: "F_θ, DMV, P_X and Epstein … pass"). P_X's weights (log p) p^{−k/2}, p ≤ X, are ≥ 0, so it passes, while N.3.3 shows Q^{P_X}
  fails Weil positivity for L ≥ L₁(X). So P_X is an unconditional, explicit instance of §C's "R5 FAILS": a positivity that P_X passes
  cannot imply OS2b. It also matches §C's sign statement: C.2d's bound for f ≥ 0 is an *upper* bound on Fl_L. For P_X the lower
  bound fails, and the failing witness f_c of N.3.3 is not ≥ 0 (it is odd). §C's superselection row "P_X: finite product, F > 0,
  no superselection; it breaks by finiteness, not by zeros" is N.3.3's mechanism: a bounded Euler operator (N.5).

## N.10 The two Tate boundary directions in the semi-local S_X picture (Golden Eye's pole lead)

Setting. S = S_X = {∞} ∪ {p ≤ X}, 𝔸_S = ℝ × ∏_{p∈S} ℚ_p, O_S^* = {±∏_{p∈S} p^{n_p}} (the S_X-units, replacing ℚ^*), ℤ[1/S] ⊂ 𝔸_S
diagonally, Tate measures (dx_p with vol ℤ_p = 1, d^*x_p with vol ℤ_p^* = 1, d^*x_∞ = dx/|x|). For f ∈ 𝒮(𝔸_S), put
Z_S(f, s) := ∫_{𝔸_S^*} f(x)|x|^s d^*x (Re s > 0). In Tate's global case the Mellin transform of Σ_{q∈ℚ^*} f(qx) has residues −f(0)
at s = 0 and f̂(0) at s = 1, and 𝒮(𝔸)₀ removes both (ARP G.4a; Golden Eye).

**N.10 Proposition (PROVED, elementary; pure tensors f = ⊗_{v∈S} f_v, Re s > 1 for the unfoldings, then meromorphic continuation).**
(i) *S_X-units only.* ∫_{𝔸_S^*/O_S^*} Σ_{u∈O_S^*} f(ux)|x|^s d^*x = Z_S(f, s) = ∏_{v∈S} Z_v(f_v, s). This is **holomorphic at s = 1**
  (it converges for Re s > 0). So **the f̂(0) boundary direction is absent.** At s = 0 it has a pole of order ≤ |S| = π(X) + 1, with
  equality when all f_v(0) ≠ 0, and the coefficient of s^{−|S|} is 2f_∞(0) ∏_{p∈S} f_p(0)/log p. It also has poles at
  s = 2πik/log p (k ≠ 0, p ∈ S, f_p(0) ≠ 0) on Re s = 0. In (ii) these are cancelled by the zeros of ∏_{p∈S}(1 − p^{−s}).
  For unramified f (f_p = 1_{ℤ_p}, f_∞ = e^{−πx²}), Z_S = π^{−s/2}Γ(s/2) P_X(s), the "completed P_X": no pole at 1 and no FE.
(ii) *Full S-integer sum.* Every ξ ∈ ℤ[1/S] ∖ {0} is uniquely u·n with u ∈ O_S^* and n ≥ 1 coprime to S, and |n|_S := ∏_{v∈S}|n|_v = n.
  Hence ∫_{𝔸_S^*/O_S^*} Σ_{ξ∈ℤ[1/S]∖0} f(ξx)|x|^s d^*x = ζ^{(S)}(s)·Z_S(f, s), with ζ^{(S)}(s) = Σ_{(n,S)=1} n^{−s} = ζ(s)∏_{p∈S}(1 − p^{−s}).
  This has **residue f̂(0) at s = 1 and −f(0) at s = 0**, exactly Tate's two boundary directions.
*Proof.* (i) is unfolding plus the local Tate integrals. Z_p(f_p, s) = f_p(0)(1 − p^{−s})^{−1} + h_p(s), with h_p entire (f_p is
locally constant with compact support). Z_∞(f_∞, s) has a simple pole at 0 with residue 2f_∞(0), and is holomorphic on Re s > 0.
(ii) The factorization ξ = u·n is unique: n = |ξ|_∞·∏_{p∈S}|ξ|_p, up to the sign absorbed in ±1 ∈ O_S^*. Substituting x ↦ x/n
gives the factor n^{−s}. At s = 1, ζ^{(S)} has residue ∏_{p∈S}(1 − p^{−1}) and Z_S(f, 1) = ∏_{p∈S}(1 − p^{−1})^{−1} ∫_{𝔸_S} f dx,
so the residue is ∫ f dx = f̂(0). At s = 0, (1 − p^{−s})Z_p(f_p, s) → f_p(0), and ζ(s)Z_∞(f_∞, s) has residue 2ζ(0)f_∞(0) = −f_∞(0),
so the residue is −∏_v f_v(0) = −f(0). ∎
**Reading (PROVED as stated; interpretation is a remark).**
- Replacing ℚ^* by the S_X-units deletes the s = 1 direction outright. It does not delete the s = 0 direction but degrades it: a
  simple pole becomes an order-|S| pole plus the imaginary companions 2πik/log p. The codimension-2 condition f(0) = f̂(0) = 0 that
  defines 𝒮(𝔸)₀ no longer removes the semi-local boundary. f(0) = 0 kills only the leading Laurent coefficient at s = 0, and f̂(0) = 0
  imposes a zero Z_S(f, 1)·Res ζ^{(S)} = 0 rather than removing a pole.
- The f̂(0) residue at s = 1 is carried **entirely by ζ^{(S)}, i.e. by the integers coprime to S, i.e. by the primes outside S_X**.
  This is the adelic form of N.4b's "unspent Ward identity". Convention N.0 re-inserts both of ζ's boundary terms (ĝ(±i/2)) by hand
  without the coprime-to-S tail that produces them, and N.3.3 shows they are then uncompensated. Equivalently (§P P.6): the bounded
  finite Euler operator L_{S_X} is the S_X-unit sum, and the missing unbounded factor ζ^{(S)} is the ℚ^*-sum over the complement.
- ℤ[1/S] is a cocompact lattice in 𝔸_S (standard; e.g. Weil, *Basic Number Theory*, theorem number unverified), so Poisson
  summation and an FE hold for the full sum in (ii). They do not hold for the S_X-unit sum in (i). This is the adelic version of
  "P_X has the Euler product but no FE" (R2).
R1 per S | R2 (ii) has both halves, (i) only the Euler half | R3 not used | R4 the s = 1 boundary (pole state) lives in the outside-S tail |
R5 DMV: no ℤ[1/S] lattice (Beurling primes), so (ii) has no analogue. F_θ, G_θ′: products of shifted copies of (ii). DH, Epstein:
not of this form. H₋ₜ: not adelic.

## References (verified this session unless marked)
- Diamond, H. G., Montgomery, H. L., Vorhauer, U. M. A., "Beurling primes with large oscillation", Math. Ann. 334 (2006) 1–36,
  doi:10.1007/s00208-005-0638-2, Theorem 1 (bibliographic data and Theorem 1 statement checked on the web this session).
- Weil, A., Comm. Sém. Math. Lund (1952); Bombieri, E., Rend. Mat. Acc. Lincei 11 (2000) 183–233 (via §S; not re-read).
- Yoshida, H., "On Hermitian forms attached to zeta functions", Adv. Stud. Pure Math. 21 (1992) 281–325, doi:10.2969/aspm/02110281.
- Connes, A., Selecta Math. 5 (1999) 29–106 (Thm 4 per ARP G.4c; §VIII eq. (16), Thm 5 eqs. (17)–(26), Lemma 3 per ARP G.4d; not re-read); Meyer, R., Duke Math. J. 127 (2005) (Lemma 5.3
  per §P, not re-read).
- Connes, A., Consani, C., "Weil positivity and trace formula, the archimedean place", Selecta Math. 27 (2021) art. 77,
  doi:10.1007/s00029-021-00689-4, arXiv:2006.13771.
- Landau, H. J., Widom, H., "Eigenvalue distribution of time and frequency limiting", J. Math. Anal. Appl. 77 (1980) 469–481
  (indicator-symbol case only; general-symbol theorem number unverified).
- Tate, J., thesis (1950), in Cassels–Fröhlich (1967): local zeta integrals (§2) and the global FE (§4), per ARP G.1a/G.4a.
- Widder (1941) Ch. VI; Berg–Christensen–Ressel (1984), via §C C.2d (not re-read).
- Suzuki, M., arXiv:2206.03682 (explicit-formula normalization, via §S).
- Team: ARP.md §0, §G (main 8079b30); §Q `team/queen/ARP_queen.md` (OS table, DMV/P_X); §C `team/celestial-cft/ARP_celestial-cft.md` (C.2d); §S `team/spectrum-measure/ARP_spectrum-measure.md`; §P `team/principal-series/ARP_principal-series.md`
  @1c46f59; bridge `work/2026-10-27_astra_prime_pnt_os_ward_identity.md`, `work/2026-10-27_astra_k2_dom_bulk_sewing.md`,
  `work/2026-10-27_astra_offline_zero_counterterm_control.md` (Astra's F_θ control: counterterm closure does not fix real parts; §N
  uses no residue/counterterm argument, and N.3.3's hypothesis excludes F_θ explicitly); `finite_windows_dmv.md`, `px_positivity.md`.
- GPPVerify `GppVerify/RiemannHypothesis/WeilSupportLadder.lean` (`namespace GppWeilLadder`): `convolution_hasSupportIn`,
  `primeSide_term_eq_zero`, `primeSide_eq_truncation`, `primeSide_eq_zero_of_support_lt_log_two`.
