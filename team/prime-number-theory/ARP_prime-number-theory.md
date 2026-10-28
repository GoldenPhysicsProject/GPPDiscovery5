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
  **GAP (DMV past the ladder), now CLOSED by N.11.3 (PROVED, ineffective):** is Q_L^{P_DMV} non-PSD for some L? No explicit formula ties the borrowed W_∞ to DMV's zeros (no FE),
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
| **DMV** (glued) | same MF | = ζ on ladder (`finite_windows_dmv.md`) | silent (full tail); **fails by N.11.3** (off-line zero ⇒ T unbounded) | small dlVP fluctuation | past the ladder, PROVED (ineffective, N.11) |
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
- **GAP N.G1 (DMV past the ladder), CLOSED in N.11 (PROVED, ineffective; effective onset is GAP N.G1′):** does the glued DMV form Q^{P_DMV}_L (borrowed W_∞, full prime tail, no FE) acquire negative
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

## N.11 GAP N.G1 attacked: glued DMV fails Weil positivity past the ladder (S103d)

**Outcome. N.G1 is RESOLVED, but ineffectively.** The glued DMV form Q^{P_DMV} (full Beurling prime tail, borrowed pole + W_∞, no FE)
has ind_− ≥ 1 for all sufficiently large L (PROVED, N.11.3; monotone in L by test-space inclusion). Under RH the onset is past the
ladder ½ log p⁺(X), by N.2. Unconditionally, only "for large L" is proved. The proof is a Laplace-transform (Landau-type) argument: any off-line zero of
ζ_P forces the smoothed ψ_P-fluctuation to be unbounded, and an even/odd two-bump witness turns that into a negative value of Q.
**No explicit L₁(DMV) is obtained** (new GAP N.G1′, N.11.6). The single-zero "quartet" witness is blocked, and N.11.5 isolates the
term that blocks it. This is the easy direction of Weil's criterion (positivity ⇒ no zeros in Re s > ½), done without the FE. It gives
nothing toward RH. **RH is OPEN.**

### N.11.1 Setting (definitions)
A *borrowed-W_∞ system* is a real measure μ on (0, ∞) with ∫ e^{−σu} d|μ|(u) < ∞ for σ > ½ and μ((0, u₀)) = 0 for some u₀ > 0, with form
  Q^μ_L(f) := ĝ(i/2) + ĝ(−i/2) − ∫_{(0,∞)} g_s dμ + W_∞(g),  f ∈ C_c^∞([−L, L]) real, g = f ∗ f̃ (g is even, so g_s = 2g).
Its Dirichlet-Laplace transform is D_μ(z) := ∫ e^{−zu} dμ(u) − 1/(z − ½), holomorphic on Re z > ½. Here −1/(z − ½) is the Laplace
transform of the subtracted mean e^{u/2}du, so D_μ is the transform of E^μ := μ − e^{u/2}1_{u>0} du. Examples: ζ (μ = μ_Λ, u₀ = log 2,
D = −ζ′/ζ(½ + z) − 1/(z − ½)); plain P_X (§N.0); glued DMV P_{X} (`finite_windows_dmv.md` §3, μ = Σ_{g-primes,k} (log p)p^{−k/2}δ_{k log p},
u₀ = log 2, D = −ζ_{P_X}′/ζ_{P_X}(½ + z) − 1/(z − ½)).
Two-bump witnesses: φ ∈ C_c^∞((−η, η)) even, φ ≥ 0, φ ≢ 0, with 2η < u₀; Φ(z) := ∫ φ(v)e^{zv} dv (entire, even); ψ_φ := φ ∗ φ ≥ 0
(support [−2η, 2η], ∫ ψ_φ(w)e^{zw} dw = Φ(z)²). For c > 2η and ε ∈ {±1}, f_{c,ε} := φ(· − c) + ε φ(· + c) ∈ C_c^∞([−(c + η), c + η]).
Smoothed fluctuation: T_μ(U) := ∫ ψ_φ(u − U) dE^μ(u), U ∈ ℝ (T_μ = 0 for U < −2η; T_μ is locally bounded).

### N.11.2 Lemma (PROVED, elementary): the two-bump identity
For c > 2η: **Q^μ(f_{c,ε}) = K(c, ε) − 2ε T_μ(2c)**, with |K(c, ε)| ≤ K_φ := 8‖φ‖² + 4ω̃_φ + 2e^{η}(∫φ)², ω̃_φ := (1/2π)∫ |φ̂|²|Ω|,
independent of c and ε. (Ω is as in N.3.2.)
*Proof.* g = 2ψ_φ + εψ_φ(· − 2c) + εψ_φ(· + 2c), and the pieces have disjoint supports since 2c > 4η. Write Q^μ = MF + W_∞ − ∫ g_s dE^μ
(N.1, same algebra). Then (i) MF(f) = (1/2π)∫ |f̂|²/(t² + ¼) ≤ 4‖f‖² = 8‖φ‖². (ii) W_∞(g) = (1/2π)∫ |f̂|²Ω with
|f̂|² = |φ̂|²|e^{ict} + εe^{−ict}|² ≤ 4|φ̂|², so |W_∞| ≤ 4ω̃_φ. (iii) On (0, 2η) ⊂ (0, u₀) the measure μ vanishes, so the central piece is
−∫_0^{2η} 4ψ_φ dE^μ = 4∫_0^{2η} ψ_φ e^{u/2} du ∈ [0, 2e^{η}(∫φ)²]. (iv) The u ≈ 2c piece is −∫ 2εψ_φ(u − 2c) dE^μ = −2εT_μ(2c). ∎
*Consequence (PROVED).* If sup_{U≥0} |T_μ(U)| = ∞, then for every K there are c and ε with Q^μ(f_{c,ε}) < 0. So **ind_−(Q^μ_L) ≥ 1 for
every L ≥ c + η, with c = U/2 for any U where |T_μ(U)| > K_φ/2**, and r^μ(L) := inf Q^μ_L/MF_L → −∞ (MF ≤ 8‖φ‖², and r^μ is
nonincreasing by inclusion of test spaces).
R1 it uses all L | R2 no FE | R3 no weight signs (|μ| only enters through convergence) | R4 for ε = +1 the witness is even and nodeless,
the soft mode's shape. For ε = −1 it is the odd pole partner (N.3.3) | R5 universal (any borrowed-W_∞ system).

### N.11.3 Theorem (PROVED, elementary; Laplace transform plus identity theorem): off-line poles force failure
Let μ be a borrowed-W_∞ system whose D_μ continues meromorphically to Re z > a for some a < ½. Suppose D_μ has a pole z₀ with
Re z₀ > max(a, 0), and choose φ with Φ(z₀) ≠ 0 (possible when η|Im z₀| < π/2, since Re Φ(z₀) = ∫ φ(v)e^{v Re z₀}cos(v Im z₀) dv > 0).
Then **T_μ is unbounded on [0, ∞), and more precisely T_μ(U) ≠ O(e^{(Re z₀ − δ)U}) for every δ > 0. Hence Q^μ is not PSD for
large L.**
*Proof.* For Re z > ½, Fubini (absolute convergence: ∫ e^{−u Re z} d|E^μ| < ∞) gives
  ∫_ℝ T_μ(U)e^{−zU} dU = ∫ e^{−zu} [∫ ψ_φ(w)e^{zw} dw] dE^μ(u) = Φ(z)² D_μ(z).
If |T_μ(U)| ≤ C e^{(Re z₀ − δ)U} on [0, ∞), the left side is holomorphic on Re z > Re z₀ − δ (the piece over [−2η, 0] is entire).
By the identity theorem on the connected set {Re z > max(a, Re z₀ − δ)} minus the poles of D_μ, Φ²D_μ would have a removable
singularity at z₀. But Φ(z₀)² ≠ 0 and D_μ has a pole there. Contradiction. With δ = Re z₀ this is unboundedness. Then apply N.11.2. ∎
*Corollaries.*
- **(a) Glued DMV, N.G1 (PROVED).** By `finite_windows_dmv.md` §3, ζ_{P_X} = ζ_P·E_X with E_X holomorphic and zero-free on σ > 0, and
  ζ_P analytic on σ ≥ θ except a simple pole at 1 (DMV, Math. Ann. 334 (2006), Thm 1). So D = −ζ_{P_X}′/ζ_{P_X}(½ + z) − 1/(z − ½) is
  meromorphic on Re z > θ − ½ (θ < 1, θ − ½ < ½). The s = 1 pole cancels: residue +1 of −ζ′/ζ against 1/(z − ½). The poles are at
  z_j = ρ_j − ½ for DMV's zeros ρ_j on σ = 1 − a/log t (this a is DMV's constant, not the continuation abscissa a of N.11.3), so Re z_j = ½ − a/log t_j > 0 for t_j > e^{2a}. The hypothesis of N.11.3 is Re z_j > max(θ − ½, 0) (continuation abscissa θ − ½);
  for θ > ½ this still holds for large t_j, since ½ − a/log t_j → ½ > θ − ½ (θ < 1), i.e. for t_j > e^{a/(1−θ)} (wording fix, Golden Eye S103e). Hence **for every X, the glued
  DMV form has ind_− ≥ 1 for some L**. That L is necessarily > ½ log p⁺(X), by N.2 and Weil for ζ, *given RH*; unconditionally it is
  just "some L". Since Re z_j → ½, **T(U) ≠ O(e^{(½−δ)U}) for every δ > 0**. So along a sequence L_k → ∞,
  Q^{P_DMV}_{L_k}(f) ≤ K_φ − 2e^{(1−2δ)(L_k−η)}, i.e. r^{P_DMV} → −∞ at the same exponential rate as plain P_X (N.3.4c), up to e^{−2δL}.
  (For each z_j, Φ_λ(z_j) = 0 for only a discrete set of dilations λ of φ, so a generic dilation works for all j at once.)
- **(b) Plain P_X (consistency with N.3.3, PROVED).** −P_X′/P_X(½ + z) is holomorphic on Re z > 0, so D_{P_X} has the pole z₀ = ½ from
  −1/(z − ½), with Φ(½) > 0. That is the unspent mean (N.4b). Explicitly T_{P_X}(U) = −e^{U/2}Φ(½)² + O_X(1), so ε = −1: the
  **odd-pole bump of N.3.3 is the ε = −1 member of the same family.** N.3.3 is the effective case. N.11 adds nothing for P_X.
- **(c) ζ (PROVED; classical direction of Weil's criterion, here without explicit formula and without FE).** Q^ζ ≥ 0 for all L ⇒ ζ has
  no zeros in Re s > ½. RH ⇒ positivity is *not* obtained (N.11.4).
- **(d) Fakes (R5).** F_θ: the same proof works with its own pole terms at s = 1 ± θ and two subtracted means e^{(½±θ)u}
  (proof sketch, Grok audit S103d: N.11.1 must be widened to two pole terms and to convergence only for σ > ½ + θ, then N.11.2–N.11.3 run unchanged on Re z > θ; routine but not written out, so label GAP-minor until written; its MF is a bounded multiplier, as in G.3c). D then has poles at z = θ + iγ (zeros ρ + θ of ζ(s − θ)), so it
  fails ✓, even under RH. DH and Epstein: their log-derivative Dirichlet series have poles in σ > 1 (Davenport–Heilbronn-type zeros in
  σ > 1; theorem numbers unverified), so they violate the hypothesis ∫e^{−σu}d|μ| < ∞ for σ > ½. **N.11.3 is silent on them**, and
  they fail elsewhere (§Q OS3, ARP G.4e). G_θ′: poles of D at ρ ± iθ′ lie on Re z = 0 under RH, so the theorem is silent and G_θ′ is **not falsely
  failed** ✓. H₋ₜ: no prime side, so it is silent.
R1 an R1 statement (needs U → ∞) | R2 no FE anywhere | R3 not reliant on weight signs (F_θ, weights > 0, fails; G_θ′, signed, is not
failed) | R4 witness family contains the soft-mode shape | R5 as in (d).

### N.11.4 What the FE adds (R2), stated precisely (PROVED unless marked)
- **Failure direction: nothing.** N.11.3 uses the prime measure, the pole cancellation at s = 1 and the bare boundedness of W_∞ on a
  fixed φ. It applies identically to systems with or without FE.
- **Positivity direction: everything.** For a borrowed-W_∞ system, "D_μ has no poles in Re z > 0" does **not** imply Q^μ ≥ 0.
  Counterexample (PROVED; a modified system without pole terms): plain P_X has no zeros in σ > 0 at all (finite Euler product), and its only pole in Re z > 0 is the
  uncancelled z = ½. Replace the subtracted mean by P_X's own: the form without pole terms, Q̃(f) = W_∞(g) − ∫ g_s dμ_{P_X}. Then
  D̃ = −P_X′/P_X(½ + z) has no poles in Re z > −½ (its poles are at Re z = −½), but Q̃ is not PSD: by N.3.2 its multiplier m_X has m_X(0) = Ω(0) − A_X < 0, and
  a dilated bump concentrated at t ≈ 0 gives Q̃ < 0 (the M-part argument of N.3.4a, without the pole). So "no off-line poles" ⇏
  positivity.
- **The exact missing identity.** For ζ, shifting the contour in Q^ζ(f) = (1/2πi)∫ ĝ(−i(s − ½))·(…) across the strip and using
  ξ(s) = ξ(1 − s) turns W_∞(g) + (left line integral) into the sum over zeros: Weil's explicit formula Q^ζ = Σ_ρ ĝ(γ_ρ). Without the
  FE, the remainder
    R_θ(f) := W_∞(g) + (1/2πi)∫_{(θ)} [ĝ(−i(s − ½)) + ĝ(i(s − ½))] (ζ_P′/ζ_P)(s) ds
  stays (given polynomial bounds for ζ_P′/ζ_P on σ = θ, which is GAP N.G1″ below). So Q^P = Σ_{Re ρ>θ} m_ρ[ĝ(−i(ρ−½)) + ĝ(i(ρ−½))]
  + R_θ, with no reason for R_θ ≥ 0. In words (PROVED direction, for ζ by Weil: FE ⇒ W_∞ is the archimedean term of μ's own explicit formula; it has no analogue for DMV.
  The converse, explicit formula with W_∞ ⇒ FE, is a Hamburger-type statement, not proved or cited here; Grok audit S103d narrowed the
  original "exactly"): **the FE makes W_∞ the archimedean term of μ's own explicit formula**, i.e. R_θ collapses to the zeros with Re ρ ≤ θ. It is
  needed only for the direction "zeros on the line ⇒ positivity", which is the RH-relevant one.

### N.11.5 Why the single-zero quartet witness φ_L = e^{−iγ₀x} sinh(β₀x) is blocked (isolated, labelled)
For ζ under ¬RH, Weil's explicit formula plus the FE quartet {ρ₀, ρ̄₀, 1 − ρ₀, 1 − ρ̄₀} lets a (1,1) vector like φ_L make the quartet
term 4Re ĝ(γ₀ − iβ₀) dominate and go negative. For DMV this needs Q^P's zero-side expansion, and three terms block it:
1. **No quartet is forced (PROVED).** Without FE the zeros are only guaranteed to come in conjugate pairs ρ, ρ̄. The (1,1) quartet structure is an FE artefact, so
   φ_L is not adapted to a symmetric zero set.
2. **R_θ (GAP N.G1″).** Q^P = Σ_{Re ρ>θ}(…) + R_θ needs bounds on ζ_P′/ζ_P along a line σ = θ₁ ∈ (θ, ½ + β₀). These are available via
   Landau's lemma from |ζ_P(s)| ≪ |t|/(σ − θ) (partial summation from N(x) = κx + O(x^θ)) together with σ = 2 lower bounds: standard,
   but not written out here, so GAP. Even then, R_θ(φ_L) = O(e^{(2θ₁−1)L}·poly(L)) and is not sign-controlled.
3. **Larger off-line zeros dominate (heuristic obstruction, not a theorem: it presupposes the zero-side expansion of N.G1″; relabelled by Grok audit S103d).** For f with f̂ of fixed shape, the pair term of ρ_j is ≍ e^{2β_j L}|Φ-type(z_j)|², and DMV
   has β_j = ½ − a/log t_j → ½ > β₀. So no single zero dominates. The zero sum has oscillating phases e^{2iγ_j c} and the size
   e^{L(1−o(1))} of the pole term. A sign for a single ρ₀ cannot be read off.
N.11.3 bypasses all three. It never expands Q into zeros: the zero enters only as a singularity of a Laplace transform. The cost is
ineffectiveness.
- **Principal series `dh_signed/` vectors** (DH cos-basis Ritz/witness vectors at L = 1.069…, 1.711, K ≤ 360): not used. They are
  DH-specific (no pole, complex weights) and are NUMERICAL witnesses at fixed L, whereas N.11.3 is analytic in U → ∞.

### N.11.5a Prior art (Grok audit S103d)
- N.11.3(c) is the easy direction of Weil's criterion (Weil 1952; Bombieri 2000, cited in References). No novelty is claimed for it.
- The N.11.3 mechanism (if the fluctuation is O(e^{σU}), its Laplace transform is holomorphic on Re z > σ, so a singularity at z₀ forces
  an Ω-result of size e^{(Re z₀−δ)U}) is the classical Phragmén–Landau/Ingham oscillation argument (Ingham, *The Distribution of Prime
  Numbers*, Cambridge 1932, Ch. V; chapter-level citation, theorem numbers unverified). Beurling-prime oscillation of this kind is also
  classical, and DMV 2006 Thm 1 already gives ψ_P − x = Ω_±. Only the Ω (not Ω_±) form is needed here: the choice of ε in N.11.2
  absorbs the sign, so Landau's theorem on Dirichlet integrals with non-negative coefficients is not used.
- What N.11 adds is local to the ARP: the two-bump identity that turns the smoothed Ω-result into ind_− ≥ 1 of the borrowed-W_∞ form,
  and the R5 classification in N.11.3(d). Scope: "every prime-side fake with an off-line zero fails" holds only for borrowed-W_∞ systems
  that meet N.11.3's hypotheses (∫e^{−σu}d|μ| < ∞ for σ > ½, meromorphic continuation to Re z > a, a pole z₀ of D_μ with Re z₀ > max(a, 0), a the continuation abscissa of the Laplace identity ∫T e^{−zU} = Φ²D_μ). DH and Epstein are
  outside them, P_X fails through its uncancelled pole rather than a zero, and F_θ is covered only by the sketch in N.11.3(d).

### N.11.6 GAPs and labels after N.11
- **N.G1: CLOSED (PROVED, ineffective).** Glued DMV past its ladder: ind_− ≥ 1 for some L, and r → −∞ along a sequence. (Via poles z_j with Re z_j > max(θ − ½, 0), which holds for large t_j for every θ < 1.)
- **GAP N.G1′ (effective onset).** Give explicit L₁(DMV; X, θ, a), i.e. an effective Ω-result for the smoothed fluctuation T(U) of
  DMV's ψ_P. Route: an effective oscillation theorem (Turán power sums / Pintz-type), fed with explicit heights of DMV zeros. DMV
  Thm 1 as transcribed gives zeros on σ = 1 − a/log t with no explicit heights (theorem-internal constants unverified).
- **GAP N.G1″ (zero-side expansion for DMV).** Polynomial bounds for ζ_P′/ζ_P on a zero-free line σ = θ₁ ∈ (θ, 1). Standard
  (Landau lemma) but not written. Needed only for the quartet-type route.
- **CONJECTURE N.C3.** ind_−(Q^{P_DMV}_L) → ∞ (expected by combining k separated pairs at U's where |T| is large. The cross terms
  T(c_i ± c_j) are not controlled).
- Main G.4d is unchanged: N.11.3 is the RH-silent direction, valid without FE.

## N.12 Prime side of the p = 2 step for Celestial's §C.13 (IF): edge-strip estimates, and where the H¹ comparison fails (S103e)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN. Nothing in N.12 proves or implies RH.**
Track tag **[L]** throughout: support extension at fixed archimedean cutoff Λ = 1, half-width L past ½ log 2 (Celestial §C.13,
origin/team/celestial-cft 407715d). Nothing here is on the [Λ, §S] track.

**Outcome.** No analytic L₀ > ½ log 2 is proved. The route "Cauchy–Schwarz on the strips, then an H¹ comparison with 𝒮_L" fails at
one named step (N.12.4, step S2): S1 is sharp, and S2 implies archimedean-only positivity, which Connes–Consani report fails inside
the range (NUMERICAL). N.12.3 gives an exact parity/reflection splitting of the (IF) left side. With it, (IF) is equivalent to the two
parity conditions (K2) and (K3), with explicit constant 2c₂ (PROVED as a reformulation). (K1) is the high-band input. Their truth on
(log 2, log 3] is CERTIFIED (consistency only) and analytically GAP (N.G2).

### N.12.0 Notation (matches §C.13.0; Celestial confirmed)
- a := log 2, c₂ := Λ(2)·2^{−1/2} = (log 2)/√2 = 0.49013 (so 2c₂ = √2·log 2 = 0.98026), Q := −∂² + ¼.
- ξ ∈ C_c^∞(−L, L) real, h = ξ⋆ξ*, h(b) = ∫ ξ(x)ξ(x + b) dx, g = Qh. ξ̂(t) = ∫ ξ e^{itx} dx as in N.0.
- Penetration ℓ := 2L − a ∈ (0, log(3/2)]. The two strips are S₋ := [−L, L − a] and S₊ := S₋ + a = [a − L, L], each of length ℓ.
  They are disjoint, since ℓ < a ⇔ 2L < 2 log 2.
- J: x ↦ −a − x reflects S₋ about its midpoint −a/2. P_± := ½(1 ± J) on L²(S₋).
- ‖u‖²_{H(I)} := ∫_I (u′² + ¼u²). E := S₋ ∪ S₊.
- R_L(ξ) := 𝒮_L(ξ) + 2ε′(1⁺)(‖ξ‖² − ⟨ξ, K_Lξ⟩) = W_∞(Qh) is the archimedean form on the pole-free class (§C.13 T3–T5).
- B(ξ) := ∫_{S₋} [ξ′(x)ξ′(x + a) + ¼ξ(x)ξ(x + a)] dx = g(a) = Qh(a). The (IF) left side is 2c₂B(ξ).

### N.12.1 The exact prime side on (log 2, log 3] (PROVED; tail vanishing PROVED-LEAN) [L]
- **Pole-free class.** For 2L ∈ (log 2, log 3], Q_L(Qh) = R_L(ξ) − 2c₂B(ξ).
  - Only n = 2, k = 1 is present: n = 3 contributes g(log 3) = 0 at 2L = log 3, and n = 4 enters at 2L > log 4.
  - Lean names: `GppWeilLadder.primeSide_eq_truncation`, `primeSide_term_eq_zero`, `primeSide_eq_zero_of_support_lt_log_two`.
  - B = g(a) is h″/h algebra: h″(b) = −∫ ξ′(x)ξ′(x + b) dx (§C.13 T7).
- **ζ-type versus P_X-type.** On this range ζ ≡ P_X for every X ≥ 2 (N.2, p⁺(2) = 3).
  - Any failure here would be ζ-type, i.e. a zero off the line: impossible under RH, and excluded by the CERTIFIED window.
  - P_X-type failure (missing primes) starts only at 2L > log 3, where P_2 lacks the n = 3 atom.
  - N.3.3's odd-pole bump for P_2 needs L ≳ ½ log 2 + 2.8 (NUMERICAL constant), so it is far-field here.
- **Semi-local reading (PROVED, from N.10 and N.2).** The window prime side is exactly S = {∞, 2}.
  - By N.10(i), the S-unit sum ±2^ℤ supplies the f(0) side (s = 0, the local pole).
  - The f̂(0) direction at s = 1 is carried only by ζ^{(S)}, i.e. by the primes ≥ 3, none of which is visible in the window.
  - On the pole-free class both directions are removed once, at §C.13 T4. So the p = 2 line has to be absorbed by archimedean data,
    with no pole direction left over. This is §C.13 C.13.2's "pole fact", seen from the prime side.
  - In Q.8.E language (§C.13.4), the q = 2^{±1} summands of E are the only source of δ_{±a}.
- R1 per window | R2 FE invisible (ζ = P_2 here, N.2′) | R3 the sign of B is the content, not weight signs | R4 see N.12.3 |
  R5 ζ, P_X (X ≥ 2) and glued DMV coincide here. F_θ has weight (2^θ + 2^{−θ})/2 times c₂ at n = 2, and G_θ′ has a cos(θ′a)
  modulation: both differ from ζ here. DH has no Euler line at a (§C.13 R5).

### N.12.2 Step S1: edge-strip H¹ bound (PROVED, elementary, sharp) [L]
- **Bound.** |B(ξ)| ≤ ‖ξ‖_{H(S₋)}‖ξ‖_{H(S₊)} ≤ ½‖ξ‖²_{H(E)}. Hence **|LHS(IF)| ≤ c₂‖ξ‖²_{H(E)}**.
  - *Proof.* Cauchy–Schwarz in H(S₋) for the pair (ξ|_{S₋}, ξ(· + a)|_{S₋}), then AM–GM. ∎
- **Comparison with the Fourier form.** ‖ξ‖²_{H(E)} ≤ ‖ξ‖²_{H(ℝ)} = (1/2π)∫(¼ + t²)|ξ̂|².
  - Localization therefore halves the multiplier constant: c₂ in place of the 2c₂|cos(ta)| of (IF-F).
- **Sharpness (PROVED).** Take φ ∈ C_c^∞(int S₋) and ξ = φ ± φ(· − a). Then B = ±‖φ‖²_{H(S₋)} and ‖ξ‖²_{H(E)} = 2‖φ‖²_H.
  - So the constant c₂ is attained, with either sign, for every frequency content of φ.
  - This is §C.13 T8 (no L² bound) in its natural H¹ form. Hence **no loss occurs in S1**.

### N.12.3 Exact splitting of the (IF) left side by parity and the strip reflection J (PROVED, elementary) [L]
Write ξ = ξ_e + ξ_o for the parity parts.
- **Parity-diagonal.** Both R_L (even multiplier) and B are parity-diagonal.
  - For B, the cross terms ∫ξ_e(x)ξ_o(x + a) + ∫ξ_o(x)ξ_e(x + a) cancel under x ↦ −x − a, and likewise for the derivative term.
- **Even ξ.** ξ(x + a) = (Jξ)(x) and ξ′(x + a) = −(Jξ′)(x) on S₋. Since d/dx swaps the J-parity, this gives
  **B(ξ_e) = ‖P₊ξ_e‖²_{H(S₋)} − ‖P₋ξ_e‖²_{H(S₋)}**.
- **Odd ξ.** **B(ξ_o) = ‖P₋ξ_o‖²_{H(S₋)} − ‖P₊ξ_o‖²_{H(S₋)}**.
- **Exact (IF).** (IF) ⇔ for σ ∈ {e, o}:
  > R_L(ξ_σ) + 2c₂‖P_{−s(σ)}ξ_σ‖²_{H(S₋)} ≥ 2c₂‖P_{s(σ)}ξ_σ‖²_{H(S₋)},  with s(e) = +, s(o) = −.
- **Harmful cones** (where the p = 2 term is negative for Q): (even, J-symmetric on the strip) and (odd, J-antisymmetric on the strip).
- **Edge expansion** (PROVED; ξ in the energy closure, edge slope s = ξ′(−L⁺)).
  - Even nodeless profile: B(ξ_e) = −s²ℓ + O(ℓ²), so the p = 2 term **helps**.
  - Odd profile: B(ξ_o) = +s²ℓ + O(ℓ²), so it **hurts at first order in ℓ**.
  - Uniform version: ‖P₊ξ_e‖²_{H(S₋)} ≤ ¼ℓ³(‖ξ′‖²_∞ + ‖ξ″‖²_∞). The harmful even mass is cubic in ℓ for C² edges; the harmful odd mass is not small.
  - Sign flip relative to f-coordinates: in k = (½ + ∂)ξ coordinates the nodeless even k is harmed. The (IF) statements are all in ξ.
- R4: the soft-mode-shaped (even nodeless) direction lies in the *helped* part at leading order. The odd pole partner is in the harmed part.

### N.12.4 The suggested H¹ route and its failing step (PROVED as stated; the onset is NUMERICAL prior art) [L]
**S2 (the comparison).** R_L(ξ) ≥ c₂‖ξ‖²_{H(E)} for all ξ. Together with S1 this would give (IF).
- **(i) S2 ⇒ R_L ≥ 0.** S2 implies archimedean-only positivity on the pole-free class at L (PROVED: the right side is ≥ 0).
- **(ii) S2 overpays by exactly 4c₂‖P_{−s}ξ‖²_{H(S₋)} per parity** (PROVED, from N.12.3). It charges the helped component at the wrong sign.
- **(iii) Prior art: S2 fails inside the range.** Connes–Consani arXiv:2106.01715 §2.2 (Figs 5–7; their L is our 2L, μ = e^{2L}).
  The following are NUMERICAL, text- and figure-level:
  - The archimedean even matrix has smallest eigenvalue ≈ 0.00133 at 2L = log 2.
  - It ceases to be positive "in the upper part of" 2L ∈ [log 2 − 0.2, log 2 + 0.2], so μ_a < e^{0.893} ≈ 2.44. The plotted range of Fig 7 is 2.05–2.30; the crossing value is not stated in the text.
  - The p = 2 term first lowers λ_min on μ ∈ (2, 2.27), then restores positivity.
  - λ_min < 6×10^{−8} at μ = 3.
  - Hence S2, and with it any route that bounds |LHS| before comparing, gives at best **2L₀ < log μ_a < log 2.44**. It never reaches log 3.
  - *Grok audit (scope of (iii); labels tightened).*
    - CC's §2.2 matrices are −W_ℝ (possibly with the rank-one w_{0,2} of CC Lemma 2.7) on the **full** even/odd test space, not on the pole-free class.
    - Loss of positivity there does not by itself give R_L < 0 on the pole-free class: with a single negative eigenvalue, the codimension-one restriction can still be ≥ 0, by Cauchy interlacing.
    - So "S2 fails inside the range" and the cap 2L₀ < log 2.44 are **NUMERICAL-indicative**, not established. (i) and (ii) are PROVED. The cap would become a statement once R_L < 0 on the pole-free class is shown at some 2L < log 2.44.
    - S2 is much stronger than R_L ≥ 0 anyway, by (ii). So the route is wasteful in any case.
- **(iv) What S2 would need on (log 2, log μ_a)** (PROVED as a reformulation): the archimedean edge coercivity
  m_L := inf_ξ R_L(ξ)/‖ξ‖²_{H(E)} ≥ c₂ = 0.4901.
  - Low band. For fixed smooth ξ, ‖ξ‖²_{H(E)} = 2s²ℓ + O(ℓ²) for even/odd edge slope s. So m_L ≥ c₂ needs an L²-type margin of R_L of size ≳ 2c₂s²ℓ.
  - No analytic lower bound on R_L for L > ½ log 2 exists in CC 2021: their T6 step is a sign at L = ½ log 2, with a computer-assisted finite-rank part. **GAP.**
  - High band. On edge pairs at frequency m, R_L ≈ σ_S(m)‖ξ‖²_{H(ℝ)} (heuristic), so S2 needs σ_S(m) ≥ c₂.
  - If σ_S(t) ∼ log(|t|/2π) (§C.13 heuristic, not checked here), the threshold is |t| ≳ 2πe^{c₂} ≈ 10.26. The non-localized (IF-F) needs 2πe^{2c₂} ≈ 16.7.
- **Replacement (exact, no loss).** (K2) ∧ (K3) ⇔ (IF) (PROVED, by N.12.3). (K1) locates the high band of any proof (heuristic).
  - (K1) high band, localized: **σ_S(t) ≥ c₂ for |t| ≥ T**. This is the per-mode constant on the sharp edge pairs of N.12.2, where R_L ≈ σ_S(m)‖ξ‖²_{H(ℝ)} (heuristic). It is not a proof step.
  - (K2) odd: **R_L(ξ_o) + 2c₂‖P₊ξ_o‖²_{H(S₋)} ≥ 2c₂‖P₋ξ_o‖²_{H(S₋)}**. First-order cost 2c₂s²ℓ on odd profiles.
  - (K3) even: **R_L(ξ_e) + 2c₂‖P₋ξ_e‖²_{H(S₋)} ≥ 2c₂‖P₊ξ_e‖²_{H(S₋)}**. The harmful mass is O(ℓ³) for C² edges, and the helped term is ≈ 2c₂s²ℓ. This is the mechanism behind CC's restoration of the even matrix.
  - **Named obstruction N.G2 (GAP): no analytic lower bound for R_L on the odd harmful cone (K2) or the even harmful cone (K3) for 2L > log 2.** The constant needed is 2c₂ = 0.98026 against the H(S₋)-mass of P_{s}ξ, equivalently c₂ = 0.49013 per strip against H(E).
- **Truth status of (K1)–(K3) on (log 2, log 3].**
  - CERTIFIED: positivity up to L = 1.42 (main 2576c89, `team/spectrum-measure/AUDIT_284.md`). Used as consistency only, not as an analytic L₀; conditionality as labelled on main.
  - NUMERICAL: CC 2106.01715.
- R1 a fixed-window statement, RH-silent | R2 the FE enters only through W_∞ in R_L | R3 the J-sign of B, not weight signs |
  R4 soft-mode shape helped, odd partner harmed | R5 fake-blind in truth (DH positive here); breaks structurally for DH (no a-line).

### N.12.5 First test of Conj Q.8.E (consistent with §C.13.4)
- Under (I1), D_L must carry −c₂S_Λ(V(2) + V(2)*)S_Λ. Here V(2) is translation by a in u = log|x|.
- (I2) on this window ⇔ (IF) ⇔ (K2) ∧ (K3).
  - *Grok audit (precision).* For a **given** D_L satisfying (I1), "(I2) ⇒ (IF)" is PROVED.
  - The converse holds only existentially: (IF) ⇒ some D_L ≥ 0 with (I1) exists abstractly (GNS, as in Q.8). That is not a statement about a specified D_L.
  - "(IF) ⇔ (K2) ∧ (K3)" is PROVED (N.12.3).
- The test Celestial can run directly: an explicit D_L must dominate the harmful cones of N.12.3 with constant 2c₂ on H(S₋). In particular, on the odd sector it needs first-order-in-ℓ margin 2c₂s²ℓ.
- Requirement: PROVED (from (I1)). Existence of D_L: CONJECTURE (Q.8.E).

### N.12.6 Prior art
- Connes–Consani 2021 (Selecta 27, art. 77, arXiv:2006.13771) Thm 1: support [2^{−1/2}, 2^{1/2}], i.e. our L = ½ log 2, pole-free class. N.12 begins exactly where the p = 2 atom enters.
- Connes–Consani arXiv:2106.01715 §2.2–2.4: this range NUMERICALLY (figures; theorem numbers not applicable).
- Yoshida 1992: small-support positivity (via ARP G.2b; theorem number unverified).

### N.12.7 GAPs and labels after N.12
- **N.G2 (GAP, new; a question of mechanism, since positivity here is CERTIFIED conditional on Mori, main 992a011; refined in N.13 into N.G2-m₀ and N.G2-bulk):** analytic L₀ > ½ log 2 for (IF). Reduced to (K2)/(K3) with explicit constants. The S1 + S2 route is capped at 2L₀ < log μ_a (NUMERICAL, μ_a < 2.44).
- N.G1′ PARKED (unchanged).
- G.4d unchanged.

## N.13 Mechanism of the p = 2 step: conditional perturbative L₀ = ½ log 2 + δ(m₀), and the p = 2 help on cross-edge profiles (S103f)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN. Nothing in N.13 proves or implies RH.**
Track **[L], Λ = 1** throughout. Notation is §C.13 and §C.13.B (team/celestial-cft 325fae5/1ccf45e) and N.12.0, in ξ-coordinates:
- a = log 2, c₂ = (log 2)/√2, ℓ = 2L − a, S₋ = [−L, −L + ℓ], S₊ = [L − ℓ, L], E = S₋ ∪ S₊.
- ‖u‖²_H = (1/2π)∫(¼ + t²)|û|², and ‖·‖_{H(I)} is its local form on an interval I.
- R(ξ) := W_∞(Qh) = (1/2π)∫(¼ + t²)Ω(t)|ξ̂(t)|² dt, with Ω = 2θ′ = Re ψ(¼ + it/2) − log π.
- Q_L(Qh) = R(ξ) − 2c₂B(ξ) on 𝒯_L := C_c^∞(−L, L) (N.12.1).

**Framing (lead steer, main 992a011).** Positivity on (½ log 2, ½ log 3] is already CERTIFIED, conditional on Mori's Prop 3.1, Lemma 3.2 and Thm N, via the 1.42 Arb replay. So (IF) and N.G2 are **questions of mechanism**: which analytic step carries positivity across the entry of the p = 2 line. They are not open questions about the truth of positivity.

**Outcome.**
- **(A) Conditional entry theorem.** If R is coercive at L = ½ log 2 with an explicit margin m₀ > 0 in H, then (IF) holds for 2L ≤ log 2 + ℓ*(m₀), with ℓ*(m₀) = 6.98×10⁻⁴·m₀·e^{−21.07/m₀} (Thm N.13.A, PROVED given m₀).
  - The odd-sector m₀ is new content: CC Thm 1 gives only ≥ 0 on the odd sector. It stays **GAP N.G2-m₀**.
  - Qualitatively, m₀ > 0 ⇔ R has no radical in the energy closure (PROVED, N.13.3).
- **(B) Prime-side p = 2 help on smooth cross-edge profiles.** Decompose ξ = ξ_b + ξ_s, with ξ_b a ramp on the strips and ξ_s in the strip interiors. Then the p = 2 term on ξ_b is the exact rank-two boundary form +2c₂α₋α₊(1/ℓ − ℓ/24) in the edge values α_± (N.13.4, PROVED).
  - The derivative cross term with ξ_s vanishes identically.
  - On ξ_b, S1 is sharp to within a factor 1 − ℓ²/8. So in the even sector the p = 2 line *pays back* the full strip energy of the ramp, and in the odd sector it charges it in full.
  - **Mechanism:** the odd entry must be carried by the archimedean margin m₀; the even entry is carried by the p = 2 line itself.

### N.13.1 Window growth: R_L has no L-dependence except its domain (PROVED) [L]
- The pole-free Weil form is R restricted to 𝒯_L, where R is a translation-invariant multiplier form independent of L (§C.13.B (3a), (3c)). The split 𝒮_L + 2ε′(1⁺)(‖·‖² − ⟨·, K_L·⟩) is a rewrite of R, so **the change of K_L with L has no effect on R itself**.
- As L grows, the only change is that the test space gains the strips. The new boundary mass is exactly the strip mass ‖ξ‖²_{H(E)}, handled in N.13.2.
- Celestial (3b): the K_I/ε part is L²-bounded and gives no H-level margin. Consistent with this, **no margin in N.13 is taken from K_L**. Every margin comes from the Ω multiplier, or from the hypothesis m₀, which is a statement about R.

### N.13.2 Explicit Ω facts (PROVED, elementary; decimals single-step NUMERICAL check `sN/n13_constants.py`) [L]
- **(Ω1) Value at 0 (closed form).** Ω(0) = ψ(¼) − log π = −γ − π/2 − 3 log 2 − log π = −5.37218.
- **(Ω2) Monotone in |t|.** Ω is even and increasing in |t|: Re ψ(x + iy) = ψ(x) + Σ_{n≥0} y²/(c_n(c_n² + y²)), with c_n = n + x.
- **(Ω3) Lower bound.** With x = ¼, y = t/2, keep the n = 0, 1 terms exactly and bound the decreasing tail n ≥ 2 below by ∫_{9/4}^∞. This gives
  > Ω(t) ≥ Ω(0) + 4t²/(¼ + t²) + 0.8t²/(6.25 + t²) + ½ log(1 + t²/20.25).
  - Hence **Ω(t) ≥ c₂ for |t| ≥ t₀ = 12.8**. This is the high-band constant, and it is analytic.
  - The value from §C.13.B, ≈ 10.27, is NUMERICAL.
  - The cruder form Ω(t) ≥ Ω(0) + log(2|t|) holds for |t| ≥ ½ (from ∫_{1/4}^∞).
  - On the σ_S form: σ_S = Ω + ε̂. A pointwise σ_S ≥ c₂ needs bounds on ε̂, which is an interface request to Spectrum & Measure (prolate side). N.13 does not need it, because it works with R = the Ω-form directly.
- **(Ω4) Dilation slope.** h(y) := y∂_y Re ψ(¼ + iy) = Σ 2y²c_n/(c_n² + y²)² satisfies h ≤ min(1 + 0.6495/y, 2ζ(3, ¼)y²), with ζ(3, ¼) ≤ 64.84.
  - The first bound is unimodal sum ≤ integral + max; the second is c_n ≥ ¼.
  - Hence h* := sup h ≤ 4.5, and **Ω(λt) − Ω(t) ≤ h* log λ** for λ ≥ 1.
- **(Ω5) Low-band strip mass.** For a band ξ^{<T} (with ξ̂ restricted to |t| < T),
  > ‖ξ^{<T}‖²_{H(E)} ≤ (2ℓT/π)·‖ξ^{<T}‖²_H.
  - This follows from |E| = 2ℓ and Cauchy–Schwarz on sup|u′| and sup|u| over |t| < T.
  - **This is the explicit ℓ-rate for the bulk**: smooth (low-band) profiles carry strip mass O(ℓ).

### N.13.3 Theorem N.13.A (conditional entry; PROVED given (M₀^σ)) [L]
**Hypothesis (M₀^σ) (GAP N.G2-m₀).** For the parity σ ∈ {e, o}: R(η) ≥ m₀‖η‖²_H for all η ∈ 𝒯_{½ log 2} of parity σ, with 0 < m₀ ≤ 1.
- **Status of (M₀^σ).**
  - Odd: ĝ(0) = 0 automatically, so CC Thm 1 gives exactly m₀ ≥ 0 (lead and §C.13.B.1). An explicit m₀^o > 0 is new content.
  - Even: positivity on the even pole-free class without ĝ(0) = 0 is not CC Thm 1 (it is CC Thm 6.11/eq. (5) territory, per §C.13.6).
  - CC's λ_min ≈ 0.00133 (2106.01715 Fig 5) is NUMERICAL, on the full test space, in their basis normalization. It is **not** an m₀ in H.
- **Qualitative m₀ (PROVED, compactness).** m₀^σ > 0 ⇔ R ≥ 0 on the energy closure V^σ with no nonzero radical vector.
  - *Proof sketch, complete.* Take ‖ξ_n‖_H = 1 with R(ξ_n) → 0, and ξ_n ⇀ ξ.
    - On |t| < T, the ξ̂_n converge boundedly and pointwise, since the supports are fixed. So the low band converges strongly.
    - The Ω₊ part is weakly lower semicontinuous. Hence R(ξ) ≤ 0, so ξ is radical, so ξ = 0.
    - Then the low band tends to 0, and R(ξ_n) ≥ Ω(T)‖ξ_n^{>T}‖² − o(1) → Ω(T) > 0. Contradiction. ∎
  - Strict positivity of R on V^σ at λ = √2 is asserted in CC 2106.01715's introduction ("known to be strictly positive"; source not re-verified at the closure level). So m₀^σ > 0 qualitatively is **CONJECTURE-strong / GAP-minor**, and its **value** is GAP.
- **Step 1 (dilation, PROVED).** For ξ ∈ 𝒯_L, ξ_λ := ξ(λ·) ∈ 𝒯_{½ log 2}, with λ = 2L/a = 1 + ℓ/a, and the parity is preserved. Changing variables and using (Ω1), (Ω4) and A_λ := ¼λ^{−2} + s² ≤ ¼ + s²:
  > R(ξ) ≥ λ^{−1}R(ξ_λ) − κ(λ)‖ξ‖²_H,  κ(λ) = h* log λ + (1 − λ^{−2})|Ω(0)|.
  - Also ‖ξ_λ‖²_H ≥ λ^{−1}‖ξ‖²_H. So R ≥ m‖ξ‖²_H on 𝒯^σ_L, with m = m₀λ^{−2} − κ(λ) ≥ m₀ − 24.88ℓ.
  - So **m ≥ m₀/2 for ℓ ≤ m₀/50**. This is the window-growth continuity, with explicit rate.
- **Step 2 (band split, PROVED).** Let m′ = m₀/2, θ = 2|Ω₀|/(2|Ω₀| + m′), and T = ½exp(|Ω₀| + 2c₂ + 4c₂|Ω₀|/m′). Write R = θR + (1 − θ)R.
  - Apply Step 1 to θR, and (Ω2)/(Ω3) by band to (1 − θ)R.
  - Apply S1 with (a + b)² ≤ 2a² + 2b² to |2c₂B| ≤ c₂‖ξ‖²_{H(E)}, and (Ω5) to the low band.
  - Then Q_L(Qh) ≥ [θm′ + (1 − θ)Ω(T) − 2c₂]‖ξ^{>T}‖²_H + [θm′/2 − 4c₂ℓT/π]‖ξ^{<T}‖²_H.
  - The first bracket is ≥ 0 by the choice of T. The second is ≥ 0 iff ℓ ≤ πθm′/(8c₂T), and θ ≥ ½.
- **Result (PROVED given (M₀^σ)).** (IF) holds on parity σ for
  > **2L ≤ log 2 + ℓ*(m₀),  ℓ*(m₀) = (π/(16c₂e^{6.3524}))·m₀·e^{−21.0645/m₀} = 6.98×10⁻⁴·m₀·e^{−21.06/m₀},**
  - i.e. **δ(m₀) = ℓ*(m₀)/2**, and L₀ = ½ log 2 + δ(min(m₀^e, m₀^o)).
  - The constants are exact closed forms: 21.0645 = 8c₂|Ω(0)|, and 6.3524 = |Ω(0)| + 2c₂.
- **Remarks.**
  - The e^{−c/m₀} rate is the cost of the crude low-band bound Ω ≥ Ω(0) in Step 2. It is not intrinsic.
  - Celestial's Lemma C.13.B (strip interiors, 2L ≤ log 2 + 0.119, NUMERICAL until certified) shows that on strip-localized vectors no m₀ is needed at all.
  - The exponential loss sits entirely in the cross-edge (bulk) part, which N.13.4 isolates.
  - This theorem is **bound-first**, and uses the p = 2 term only through |B|. It does not identify the mechanism. N.13.4 does.
- R1 a fixed range, RH-silent | R2 FE through Ω only | R3 |B| only, no weight signs | R4 odd/even treated separately | R5 fake-blind in truth. DH has no line at a, so the theorem is not about it.

### N.13.4 Deliverable for Celestial: the p = 2 help on smooth cross-edge profiles (PROVED, elementary) [L]
**Decomposition.** For ξ ∈ 𝒯_L (or its H¹_0 closure), let α₋ := ξ(−L + ℓ) and α₊ := ξ(L − ℓ), the values at the inner strip edges.
- ξ_b := ξ on [−L + ℓ, L − ℓ]; ξ_b(x) := α₋(x + L)/ℓ on S₋; ξ_b(x) := α₊(L − x)/ℓ on S₊. These are ramps, continuous, with ξ_b(±L) = 0.
- ξ_s := ξ − ξ_b = φ₋ + τ_aφ₊, with φ_± ∈ H¹_0(int S₋) (Dirichlet at both ends of each strip).
- **Lemma C.13.B applies to ξ_s.** It is stated for C_c^∞(int S); it extends to H¹_0 by density in the energy norm.
- Parity: for even ξ, α₊ = α₋ = α; for odd ξ, α₊ = −α₋ = −α.

**Exact p = 2 identity.** With B(u, v) the polarization of B:
> **2c₂B(ξ) = 2c₂B(ξ_s) − 2c₂α₋α₊(1/ℓ − ℓ/24) + (c₂/2)X,**
> X := ∫_{S₋}[ξ_b(x)ξ_s(x + a) + ξ_s(x)ξ_b(x + a)] dx,  |X| ≤ (ℓ^{3/2}/(π√3))(|α₋|‖φ₊′‖ + |α₊|‖φ₋′‖).
- *Proof.*
  - On S₋, ξ_b(x) = α₋y/ℓ and ξ_b(x + a) = α₊(ℓ − y)/ℓ, with y = x + L. So B(ξ_b) = ∫_0^ℓ[−α₋α₊/ℓ² + ¼α₋α₊y(ℓ − y)/ℓ²] dy = α₋α₊(ℓ/24 − 1/ℓ).
  - In 2B(ξ_b, ξ_s), the derivative terms are (α₋/ℓ)∫_{S₋}ξ_s′(x + a)dx − (α₊/ℓ)∫_{S₋}ξ_s′(x)dx = 0, because ξ_s vanishes at all strip endpoints.
  - The remaining ¼-term is ¼X, bounded with ‖ξ_b‖_{L²(S)} = |α|√(ℓ/3) and the Dirichlet Poincaré inequality ‖φ‖ ≤ (ℓ/π)‖φ′‖. ∎

**Sharpness on the ramp (PROVED).** ‖ξ_b‖²_{H(S₋)} = α₋²(1/ℓ + ℓ/12), and the same holds on S₊ with α₊. Hence on the ramp
> −σB(ξ_b) / (½‖ξ_b‖²_{H(E)}) = (1 − ℓ²/24)/(1 + ℓ²/12) ≥ 1 − ℓ²/8,  σ = +1 (even), −1 (odd).
- **The bulk edge trace of a smooth profile is an almost-extremal vector of S1.**
  - Even sector: the p = 2 line returns c₂(1 − ℓ²/8)·‖ξ_b‖²_{H(E)} = 2c₂α²(1/ℓ − ℓ/24) as **help**.
  - Odd sector: the same amount is **cost**. This is N.12.3's 2c₂s²ℓ, with α = sℓ + O(ℓ²).
  - There is no smaller cost to find in the odd sector. The 2c₂α²/ℓ is exact, so it must be paid by R.

**Interface form for the {∞, 2} D_L test on (K2)/(K3).**
> Q_L(Qh) = [R(ξ_s) − 2c₂B(ξ_s)] + [R(ξ_b) + 2σc₂α²(1/ℓ − ℓ/24)] + [2R(ξ_b, ξ_s) − (c₂/2)X].
- **First bracket:** ≥ β_ℓ(‖φ₋‖²_H + ‖φ₊‖²_H), with β_ℓ > 0 for ℓ ≤ 0.119 (Lemma C.13.B, NUMERICAL until certified).
- **Second bracket:** the bulk, with the q = 2^{±1} line acting as a **rank-two boundary form in (α₋, α₊)**, namely −2c₂(1/ℓ − ℓ/24)·α₋α₊ subtracted.
  - In D_L language, the 2^{±1} summands of E restricted to the edge traces are the matrix −c₂(1/ℓ − ℓ/24)·[[0, 1], [1, 0]] in (α₋, α₊).
  - Its eigenvectors are the even (help) and odd (cost) edge pairs.
- **Third bracket:** (c₂/2)|X| ≤ c₂[ε'α²/ℓ + ℓ⁴‖φ′‖²/(12π²ε')] for any ε' > 0 (AM–GM).
  - This is absorbed by a fraction ε' of the second bracket's α²/ℓ and by the first bracket's β_ℓ, as soon as ℓ⁴/(12π²ε') ≪ β_ℓ.
  - What remains is **2R(ξ_b, ξ_s)**: the nonlocal archimedean coupling between the bulk ramp and the strip interior. This is the irreducible part of **N.G2-bulk (GAP)**, now isolated as one bilinear term.
- **Reduced conditions.**
  - **(K3)-bulk:** R(ξ_b) + 2R(ξ_b, ξ_s) ≥ −(2 − ε')c₂α²(1/ℓ − ℓ/24) − (β_ℓ-absorbable terms). The budget is **positive**, supplied by the p = 2 line.
  - **(K2)-bulk:** R(ξ_b) + 2R(ξ_b, ξ_s) ≥ (2 + ε')c₂α²(1/ℓ − ℓ/24) − (…). The budget is **negative**, so it must come from R, i.e. from m₀^o (Thm N.13.A Step 1 gives m ≥ m₀/2 against ‖ξ_b‖²_H ≥ ‖ξ_b‖²_{H(E)} ≈ 2α²/ℓ).
  - Conditional (K2)-bulk (PROVED given (M₀^o)): if R ≥ m‖·‖²_H with m ≥ (2 + ε')c₂ on the span, the odd bulk closes. Since m₀ ≪ c₂ is expected, this is **not** how the odd entry closes for small m₀. The Thm N.13.A route (band split, exponential δ) is the one that works.
- **Mechanism statement** (PROVED as a reformulation; it answers "which step carries positivity"):
  - **Even entry:** carried by the p = 2 line, through the sign of B on the ramp. Archimedean data only need to control R(ξ_b) + 2R(ξ_b, ξ_s) from below by a *negative* budget.
  - **Odd entry:** carried by the archimedean margin, through R's coercivity at ½ log 2 and its dilation continuity (N.13.3 Step 1). The p = 2 line is a pure cost there, at exactly 2c₂α²/ℓ.
  - This matches §C.13.B.1's side remark: at ℓ = 0 the odd (K2) is CC Thm 1, so the odd entry is the CC-margin question.
  - It also matches CC 2106.01715 Fig 7: p = 2 restores the even matrix (NUMERICAL, full space).
- R1 one rule per ℓ | R2 FE via Ω in R | R3 the sign of B is used (σ-dependent), not |weights| | R4 even nodeless is helped; the odd pole partner is charged | R5 DH: no line at a.

### N.13.5 Interface requests (no duplication)
- **Spectrum & Measure (archimedean/prolate):** an explicit m₀^o > 0 with R(η) ≥ m₀‖η‖²_H on odd η ∈ 𝒯_{½ log 2} (strict CC Thm 1 in H), and the even analogue. Either gives δ by Thm N.13.A. A frequency-resolved version (a margin on |t| < T only) would remove the e^{−21/m₀} loss.
- **Celestial (C.14):** a bound on the single coupling term 2R(ξ_b, ξ_s) in terms of α²/ℓ and β_ℓ‖φ‖²_H. With it, (K3)-bulk closes from the p = 2 budget alone (expected; CONJECTURE).

### N.13.6 Labels after N.13
| item | label |
|---|---|
| R_L = R restricted to 𝒯_L; K_L dependence is internal | PROVED |
| (Ω1)–(Ω5); t₀ = 12.8 for Ω ≥ c₂ | PROVED (decimals NUMERICAL single-step) |
| Thm N.13.A, δ(m₀) = ½·6.98×10⁻⁴·m₀e^{−21.06/m₀} | PROVED given (M₀^σ) |
| m₀ > 0 ⇔ no radical (compactness) | PROVED; strictness on the closure GAP-minor |
| value of m₀^σ | **GAP N.G2-m₀** |
| N.13.4 decomposition, exact boundary form, ramp sharpness 1 − ℓ²/8 | PROVED |
| Lemma C.13.B (cited) | NUMERICAL until certified |
| coupling 2R(ξ_b, ξ_s) | **GAP N.G2-bulk** (irreducible term) |
| (K3)-bulk closes from the p = 2 budget | CONJECTURE |
| positivity on the range itself | CERTIFIED (conditional on Mori Prop 3.1, Lemma 3.2, Thm N; 1.42 Arb replay; main 992a011) |

## References (verified this session unless marked)
- Diamond, H. G., Montgomery, H. L., Vorhauer, U. M. A., "Beurling primes with large oscillation", Math. Ann. 334 (2006) 1–36,
  doi:10.1007/s00208-005-0638-2, Theorem 1 (bibliographic data and Theorem 1 statement checked on the web this session).
- Weil, A., Comm. Sém. Math. Lund (1952); Bombieri, E., Rend. Mat. Acc. Lincei 11 (2000) 183–233 (via §S; not re-read).
- Yoshida, H., "On Hermitian forms attached to zeta functions", Adv. Stud. Pure Math. 21 (1992) 281–325, doi:10.2969/aspm/02110281.
- Connes, A., Selecta Math. 5 (1999) 29–106 (Thm 4 per ARP G.4c; §VIII eq. (16), Thm 5 eqs. (17)–(26), Lemma 3 per ARP G.4d; not re-read); Meyer, R., Duke Math. J. 127 (2005) (Lemma 5.3
  per §P, not re-read).
- Connes, A., Consani, C., "Weil positivity and trace formula, the archimedean place", Selecta Math. 27 (2021) art. 77,
  doi:10.1007/s00029-021-00689-4, arXiv:2006.13771.
- Connes, A., Consani, C., "Spectral triples and ζ-cycles", arXiv:2106.01715, §2.2–2.4 (read this session; figure-level numerics).
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
