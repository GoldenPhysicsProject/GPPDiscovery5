# C.18 shared entries: spec for reuse (Temple route and Schur / Thm S routes)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Credits: Principal (P.17 Temple lemmas), Prime (N.17 kernel/G1 lemmas), Golden Eye (the Temple idea), S&M (Thm S, N_N and n_perp). RH is OPEN.

**Files.**
- `c18_shared_entries_odd_log3.json`: 2L = log 3, primes {2}. At 2L = log 3 exactly, the n = 3 atom sits at the window edge and contributes 0.
- `c18_shared_entries_odd_log3p03.json`: 2L* = log 3 + 0.03, primes {2, 3}.
- Generator: `c18_temple_arb.py`, which uses the engine `c18_entries.py`. Run `python3 c18_temple_arb.py` (about 1 min on 1 core) to regenerate both files.

**Encoding.** Every number is a pair of strings `[mid, rad]`, meaning |x − mid| ≤ rad. These are Arb balls written by `arb.str(30, radius=True)`, which is itself rigorous. Read one back as `arb(mid) + arb(0, rad)`, or as `arb("[mid +/- rad]")`.

**Basis.** This is the odd sector:
- f_k = sin(w_k x)/(√L n_k) for k ≥ 0, with w_k = (k+1)π/L and n_k = √(w_k² + 1/4).
- The f_k are H-orthonormal in H¹₀(−L, L), with ‖ξ‖²_H = ∫ ξ′² + ¼ξ².
- Prime's basis is e_n = sin(w_k(x+L))/√L with k = n−1. Since sin(w_k(x+L)) = (−1)^{k+1} sin(w_k x), we get e_n = (−1)^{k+1} n_k f_k, and Prime's Gram is G_N = diag(n_k²). So Prime's matrices are D·(ours)·D with D = diag((−1)^{k+1} n_k).

**Keys.**

| key | content |
|---|---|
| `T`, `n_square` = 64, `M_rows` = 800 | cutoff T = 150; square block size; rectangular row count |
| `Omega_T`, `q_T`, `primes[].c_n` | Ω(T); q_T = Ω(T) − 2Σc_n; c_n = log p/√p |
| `NOm_rows` (800 × 64) | N^Ω_{jk} = (1/π)∫₀^T (¼+t²)(Ω(T)−Ω(t)) f̂_j f̂_k dt, for rows j < 800 and columns k < 64. This is **S&M's N_N** (m_neg Gram). |
| `B_rows[p]` (800 × 64) | B^p_{jk} = (1/π)∫₀^T (¼+t²)(1−cos(t log p)) f̂_j f̂_k dt |
| `A_square`, `M_square` (64 × 64) | A = q_T I − N^Ω + Σ 2c_p B^p (the Gram of the comparison form Q′_T); M = N^Ω + Σ 2c_p B^p (the m_abs Gram) |
| `eps[K].abs`, `eps[K].neg`, K ∈ {64, 128, 256, 800} | ε^g_K = (1/π)∫₀^T (¼+t²) g K_{W_K} dt, with g = m_abs or m_neg, and K_{W_K} = K_o − Σ_{k<K} f̂_k², where K_o is Lemma C.16.K. **S&M's n_⊥ (for N = K) equals eps[K].neg.** |
| `temple` | The certificate data: c_dyadic_int (c_j = int·2^−48), ρ, s_v, r2_head = Σ_{k<800} r_k², τ², ν₂, δ, s_V, μ, certified. |

Here f̂_k(t) = √L[sinc((w_k−t)L) − sinc((w_k+t)L)]/n_k.

**Enclosure method.** The method is PROVED and described in the docstring of `c18_entries.py`:
- Panel Gauss–Legendre: unit panels, 40 nodes per panel, with nodes and weights enclosed by `arb.legendre_p_root`.
- Error bound: (h/2)·8·M_p ρ^{−80}/(ρ−1), where ρ = (b + √(b²+h²/4))/(h/2) and b = 0.4.
- M_p comes from acb ball evaluation on the boundary of a rectangle that contains the Bernstein ellipse. This uses the maximum modulus principle; the integrands are holomorphic for |Im t| < 1/2.
- The worst entry radius is below 2e−20.

**Validation.**
- `c18_replay_c16.out`: at T = 40, N = 96, 2L = log 3, all 2·96² entries of A and M, plus ε^abs and ε^neg, overlap the independent `acb.integral` entries in `c16_odd_arb_entries.json` (257807 CPU-s). The max |mid diff| is 5e−30.
- Even ℓ = 0.1, T = 60: A_jk for j, k < 4 agrees with `c15_entries_ref.out` (mpmath, 12 digits) in every printed digit.
- q_T(150) agrees with `p15_p16_check.log` part 7.

**How to reuse.**
1. **Temple (C.18).** Read `temple`. Recompute ρ = cᵀAc/|c|² from `A_square`. Compute r_k from rows k < 800, where A_rows = −NOm_rows + Σ 2c_p B_rows[p], plus q_T on the diagonal for k < 64. Then τ² ≤ s_v·eps[800].abs, and μ = ρ − (Σr² + τ²)/(ν₂ − ρ).
2. **Schur (Thm C.15.A / C.16.A).** For N = 64, use A_square, M_square and eps[64]. The result is a valid but weak certificate. For larger N, regenerate the square blocks with `Engine(...).gram('Om', range(N), range(N))` (seconds) and use eps[128] or eps[256].
3. **S&M Thm S.** N_N = NOm (square block), n_⊥ = eps[N].neg, and the prime lines are B^p. All of them are cut at T. **REMARK/GAP for G1 (Prime N.17):** the full-Q quantities (the [T, ∞) Ω-tail of A_N, Prime's W, E with κ_o = coth(L/2), and μ_⊥) are **not** in these files. This route caps q at q_T for t ≥ T and integrates nothing beyond T. A full-Q route needs its own tail enclosure.

**Label.** The entries are CERTIFIED as enclosures. They are inputs only; the theorems they feed are stated in ARP C.18.
