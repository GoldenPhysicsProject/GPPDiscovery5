# S12 — Candidate axiom X (local PD-preservation + product formula ⇒ reflection positivity?) — verdict: OBSTRUCTION

**Axiom X (candidate).** L(s)=∏_v L_v(s) with (X1) every finite local factor PD-preserving: r ↦ L_p(σ+ir)/L_p(σ) positive definite for σ>0 (⇔ unitary Satake, S9/S10, PROVED); (X2) archimedean factor Γ_R with its Lévy/digamma jump law (S-Lévy, classical); (X3) product formula/FE; (X4) single pole. Desired conclusion: Weil form W(f*f̃) ≥ 0 on C_h = {f*f̃ : supp f ⊂ [−h,h]} for some h > ½log2 (past Yoshida/Connes–Consani window, where prime 2 enters).
Falsifiers: F_θ violates X1 (S9/S10), DH violates X3 Euler product. X is satisfied by ζ. So X separates — the question is only whether X ⇒ positivity.

**Lemma S12.1 (PROVED, elementary): local PD is destroyed by log-differentiation.** The Weil prime term at p has Fourier symbol S_p(r)= −2 log p Σ_m p^{-m/2} cos(m r log p): an almost-periodic series with zero mean, so it is sign-indefinite (min at r=0: −2log p·p^{-1/2}/(1−p^{-1/2}); max at r=π/log p: +2 log p·p^{-1/2}/(1+p^{-1/2})). PD of the multiplicative factor L_p does not pass to −(log L_p)′; every local Weil term is indefinite on all of L².

**Lemma S12.2 (PROVED): local (Gershgorin) bound.** For a ≥ h, sup_{f∈L²[−h,h]} |⟨f,τ_a f⟩|/‖f‖² ≤ 1/2. Hence W ≥ 0 on C_h follows from λ₀(h) := λ_min(pole+arch on C_h) ≥ Σ_{n<e^{2h}} Λ(n) n^{-1/2}.

**S12.3 (NUMERICAL, float Ritz, M=14, output.txt).** That sufficient condition fails by >2 orders and fails structurally: λ₀(h) (pole+arch alone) is itself NEGATIVE for h ≥ 0.40 (odd) and h ≥ 0.45 (even): −1.2e-2, −1.4e-1, −2.5e-1, −3.4e-1 (odd, h=.40,.45,.50,.54), while the full ζ form stays positive (+3.1e-2 … +1.0e-4). On the bottom vector of pole+arch the prime term contributes POSITIVELY (−Σ Λ n^{-1/2}F(log n) = +5.8e-2…+4.1e-1).

**Obstruction (precise).** Past the prime-2 threshold the archimedean+pole part is not coercive — it is indefinite — and the Euler-product primes are what restore positivity. So no argument of the form "local PD at each place + pole/arch dominance" can work on any C_h with h ≳ 0.40; X ⇒ positivity would need the primes to act as a *positive corrector aligned with the archimedean negative directions*, i.e. exactly Turbo's open "prime-coupled, zero-blind positive pairing". X as stated (local data + product formula) does not encode that alignment. OPEN whether a strengthened X (archimedean–prime compatibility) is provable without RH strength.

Does NOT claim: anything about RH (OPEN); Ritz values are upper bounds, so pole+arch negativity is a genuine witness (float) but full-form positivity is not proved; classical-adjacent (cf. Connes–Consani: archimedean positivity confined to small window).
