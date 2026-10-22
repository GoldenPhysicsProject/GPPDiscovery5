## H56b: why the knockout margin is not c/log p (edge-layer lemma)

Setting. Q_ζ(f) = arch + pole − Σ_n w_n Re g(log n), with g = f ⋆ f̃ (autocorrelation), w_n = 2Λ(n)/√n, f ∈ L²[−L, L].
Knocking out p means removing the n = p term: Q_{−p}(f) = Q_ζ(f) + w_p Re g(log p). Write L = ½log p + δ with 0 < δ < ½log p, so p² is invisible.

(1) PROVED (elementary). Edge-layer identity. g(log p) = ∫_{L−2δ}^{L} f(x) f̄(x − log p) dx involves only the two edge layers of width 2δ. If h(y) = f(L − y) on [0, 2δ], then
- odd f: g(log p) = −∫_0^{2δ} h(y) h̄(2δ − y) dy = −(h ⋆ h)(2δ);
- even f: g(log p) = +(h ⋆ h)(2δ).
So the knockout hurts odd f whenever the layer profile keeps its sign (automatic for small δ), and it hurts even f only when h changes sign inside a layer of width 2δ. That is why the odd sector always fails first, and why the even margin is larger (S54/S55 data).

(2) PROVED (pen-and-paper, constants not made explicit, Lean-open). δ_p > 0 strictly whenever Q_ζ is strictly positive slightly past ½log p:
- (a) Layer bound: for any set E ⊂ ℝ of measure |E| < 1/e, ∫_E |f|² ≤ (C/log(1/|E|))·(A_log(f) + ‖f‖²), where A_log(f) = ∫|f̂(t)|² log(1+|t|) dt/2π. Proof: split f̂ at |t| = N. The low part has ‖·‖_∞² ≤ (N/π)‖f‖². The high part has ‖·‖² ≤ A_log/log(1+N). Take N = 1/(|E| log(1/|E|)).
- (b) Coercivity: arch ≥ A_log − C₁‖f‖² (digamma symbol ~ log|t|). The pole term and the finitely many visible prime terms are bounded on L²[−L, L]. So Q_ζ ≥ A_log − B‖f‖².
- (c) By Cauchy–Schwarz, |g(log p)| ≤ ∫_E |f|², where E is the two layers, |E| = 4δ. Combine with Q_ζ ≥ λ‖f‖², write Q_ζ = θQ_ζ + (1−θ)Q_ζ, and take θ = w_p ε, with ε = C/log(1/(4δ)). The result is Q_{−p} ≥ [(1−θ)λ − θ(B+1)]‖f‖², which is > 0 once
  δ < ¼·exp(−C w_p (B+1+λ)/λ).
- Hence δ_p ≥ ¼ exp(−C w_p (B+1)/λ_min) > 0. For p ≤ 13 the input λ_min ≥ 1.2e−69 is CERTIFIED (S55: full complex up to L = 1.3505 > ½log 13). So δ_p > 0 is rigorous there, but the lower bound is astronomically small, about exp(−10^70).

(3) Consequence. The true margin lies in (exp(−C w_p B/λ_min), Ritz upper bound]. Because λ_min of the Weil form is itself super-exponentially small, the lemma permits δ_p to be essentially zero. In that case primes would be "barely necessary": the moment the support passes log p, removing p breaks positivity, up to a margin controlled by λ_min. This is a CONJECTURE. φ absent.
