# Round 34 — attacking R_c by ratio classes (Heart of the Cosmos / D5)

Daniel's go: "Do your thing" (2026-11-21). Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN.
No zeros are proof inputs. Weil positivity is an output, never the target. **RH is OPEN.**

**Builds on Codex's filing** (credited, not duplicated):
- GPPDiscovery2 `8ef8374`, `experiments/short_reciprocal/prime_ratio_control.js`, "Regroup Fejer off-diagonals by reduced ratio and stress-test partial inversion";
- the bridge entry "Codex — Prime-ratio contribution to Heart R33";
- Turbo's audit (PR #15 comment 7062427004, bridge `70ac6ad`).

Everything below was recomputed independently in C/Python from the definitions. The JS was read, not run.

## 0. Verdict

- **Codex's (a)–(d) all stand**, with Turbo's Ω wording for (c). (b) and (d) reproduce to every printed digit.
- **The regrouping goes one level finer.** Sort the ratio classes v/u into **layers k = ω(u)+ω(v)**.
  - Every layer has the fixed sign (−1)^k (PROVED).
  - The prime class is exactly layer 1.
  - The layer generating polynomial E_X(z) = Σ_k E_k z^k is an average of Fejér squares over random multiplicative signs, so E_X(z) ≥ 0 on [−1, 1] (PROVED).
  - E_X(1) = G is real μ; E_X(0) = H is the Haar (counted-once) part; E_X(−1) = G_{|μ|} is the Codex/Turbo translate.
- **The payback is X-scale, not √X/log X.** At X = 10⁶ the layers are:
  - E_0 = +155.9 (H)
  - E_1 = −4172.0 (prime class)
  - E_2 = +5.94e5
  - E_3 = −2.67e6
  - E_4 = +5.04e6
  - E_5 = −5.08e6
  - …
  - E_11 = −0.4

  Their alternating sum is G = 5.43, while Σ|E_k| = 1.75e7. Codex's "all other ratios +4021.58" is itself the difference of two layer totals, each 8.7e6. The prime class is 0.05% of the negative mass.
- **The identity that forces the payback (PROVED, Euler-product algebra):**
```
Σ_{m,n} μ(m)μ(n) z^{k(m,n)} m^{-a} n^{-b} = ζ(a)^{-z} ζ(b)^{-z} ζ(a+b)^{1-z²} ζ(2a)^{z(1-z)/2} ζ(2b)^{z(1-z)/2} · H_z(a,b),
```
  H_z is holomorphic and nonzero for Re a, Re b > 1/3, and H ≡ 1 at z = ±1. At z = 1 every non-integer singular exponent disappears: what is left is 1/(ζ(a)ζ(b)), and that is analytic at a, b = 1. This is Σ_{d|n} μ(d) = [n = 1] at the level of classes.
  - The X-scale mass of each layer comes from the Selberg–Delange singularity of ζ^{−z} at s = 1, carrying the factor 1/Γ(−z)², which vanishes at z = 0 and z = 1.
  - Numerically, E_X(z) tracks L·4X·A(z)²J(z)(log X)^{−2z−2}/Γ(−z)² (ratio 0.88 at z = −1 and 0.96 at z = −0.5 at X = 10⁷; it converges more slowly near the zeros at z = 0, 1).
  - **The payback is the vanishing of all Selberg–Delange main terms at z = 1.**
- **Which identity the twist breaks.** The Codex twist (flip μ(p) on (X/2, 3X/4]) breaks the inverse identity at 18,700 primes, so it is no longer the inverse of a fixed Dirichlet series:
  - Meissel: Σμ′(n)⌊X/n⌋ = 37,401 instead of 1;
  - the σ = ½ prefix F(1) goes from 0.207 to 47.78;
  - Σμ′(n)/n = 0.0609, which still passes the elementary |Σμ/n| ≤ 1 but fails the PNT-level O(1/log² X).
  - **PNT-level information already kills this twist**; RH-level information is not needed.
- **Controls (each with its own μ, X = 3·10⁴).**
  - The same payback appears for ζ², F_0.1 and Beurling (multiplicative, so their layers alternate too), and approximately for EPS and DH.
  - **The payback is generic to Dirichlet inverses.** What separates ζ is the residual after the payback, which is R33's GAP: G grows exponentially for F_θ, EPS and DH (R33 table).
- **GAP-34** (below) is the exact alternating-layer statement. It is ⇔ RH in the ε-form and strictly weaker than (implied by) uniform-in-shift Chowla. The class route reproduces the Vinogradov–Korobov bound unconditionally.
  - Its first unattained milestone is E_X(1) ≪ X^{1−δ}. That is a zero-free strip Re s > 1 − δ/2 (quasi-RH), which the Selberg–Delange expansion around z = 1 cannot supply, because its error terms are the zeros.

## 1. Verification of Codex's (a)–(d)

| claim | status | my recompute |
|---|---|---|
| (a) Reduced-ratio grouping: weight depends only on \|log v/u\|; μ(gu)μ(gv) = μ(g)²μ(u)μ(v) for (g, uv) = 1, else 0; the n = 1 centering is the single term −2c(L − log v)μ(v)/√v on u = 1; reflection keeps the sign | **PROVED** (unique factorization). CHECKED: the brute-force pair loop (`r34_pairs.c`) equals the fast gcd/d-decomposition (`r34_layers.c`) at X = 3·10⁴ to all printed digits (E_0 94.32, E_1 −1013.70, R −89.22, G 5.1062). Sign check on 4.4e5 pairs: 0 mismatches | Holds for any multiplicative μ with sgn μ(p^e) = (−1)^e (ζ, ζ², F_θ, Beurling). **Fails for EPS** (non-multiplicative): 36,585 mismatches, e.g. μ(4)μ(6) = 2 but μ(2)μ(3) = 0. DH: no mismatch in the tested range, but not proved |
| (b) X = 10⁶: prime pairs −4172.01, others +4021.58, total −150.42 | **CHECKED exactly** | E_1 = −4172.0076, R = −150.42, H = 155.85, G = 5.4304. Also X = 10⁴: −639.69 / −71.62; X = 10⁵: −1656.23 / −107.49 (= Turbo's recompute) |
| (c) \|prime subtotal\| ≳ √X/log X while R ≥ −H, H ≲ log² X, so compensation is necessary | **PROVED as Ω(√X/log X)** (Turbo's wording; ">>" is not proved). Keep g = 1 and primes in (X/2, 3X/4]; Chebyshev suffices. \|E_1\|/(√X/L) = 58.9, 60.3, 57.6, 55.8, 53.8 (X = 10⁴ … 10⁷): flat-to-slowly-falling | **Refinement (PROVED):** the compensation is of order X, not √X/log X. Σ_k\|E_k\| = E(−1) = G_{\|μ\|} ≥ Q_{\|μ\|} ≫ X, and E(1) = O(G) is small, so the even and odd layer totals are each ≍ X (8.74e6 at 10⁶). Since there are ≤ 2ω_max(X)+1 ≪ L/log L layers, some layer is ≫ X·log L/L |
| (d) Flip μ(p) on (X/2, 3X/4]: H unchanged; R = +30902.18 at 10⁶; eventually positive | **CHECKED exactly; proof PROVED** (Turbo-audited: M(y) = O_A(y/log^A y), A = 2) | R_twist = +340.87, +3527.53, +30902.18, +264548.69 (X = 10⁴, 10⁵, 10⁶, 10⁷); H identical. Consistency: G_twist ≈ L·F′(1)² = 13.8·47.78² ≈ 3.15e4 vs 3.106e4 |

**The prime class has a closed form (PROVED; CHECKED).**
```
E_1 = −2 Σ_{p≤X} (L − log p) p^{−1/2} (S_p(X/p) − c)
S_q(y) = (6/π²) Π_{p|q} p/(p+1) · [log y + γ − 2ζ′(2)/ζ(2) + Σ_{p|q} log p/(p+1)] + O(y^{−1/2})
```
- S_q(y) is the residue of ζ(1+w)/ζ(2+2w)·Π_{p|q}(1+p^{−1−w})^{−1}·y^w/w at w = 0. CHECKED: the error × √y stays ≤ 0.06 for q = 2, 3, 30, 210.
- With the main term only, the sum gives −4163.20 against the exact −4172.01 at 10⁶.
- Leading asymptotic: E_1 ≈ −2(√X/L)[(6/π²)(16 + 4C₀) − 4c], C₀ = γ − 2ζ′(2)/ζ(2) = 1.7171. Its 1/L corrections are large: the ratio to the leading term is 1.62–1.81.

## 2. The class/layer structure (PROVED)

**(I34.1) Class form (Codex).** With a_n = μ(n)n^{−1/2} − c[n = 1] and K = L − |log(n/m)|:
```
R_c = Σ_{coprime u<v} 2 (L − log(v/u)) μ(u)μ(v) (uv)^{−1/2} S_{uv}(X/v)  −  2c Σ_{v≥2} (L − log v) μ(v) v^{−1/2},
S_{uv}(y) = Σ_{g≤y, (g,uv)=1} μ(g)²/g .
```

**(I34.2) Layers.**
- Put k(m,n) = ω(u)+ω(v) = ω(m Δ n) and E_k = Σ_{k(m,n)=k} a_m a_n K(m,n). Then E_0 = H_c and R_c = Σ_{k≥1} E_k.
- **Sign:** for squarefree-supported multiplicative μ with c < 1, every pair in layer k has sign (−1)^k, so λ_k := (−1)^k E_k ≥ 0.
  - Each class contributes μ(u)μ(v)·(positive) = (−1)^k·(positive).
  - The g = 1 term of class (1, v) carries a_1 = 1 − c > 0.
  - CHECKED for X = 10⁴ … 10⁷, with DFT noise ≤ 3·10⁻³ absolute on the vanishing top layers at 10⁷.

**(I34.3) Random-sign representation.** Let ε be a random completely multiplicative ±1 function with independent ε_p of mean z ∈ [−1, 1]. Then 𝔼[ε(m)ε(n)] = z^{k(m,n)}, so
```
E_X(z) = Σ_k E_k z^k = 𝔼_ε ‖f_c(τ_ε ·)‖²_{ν_L} ≥ 0   on [−1, 1],
E_X(1) = G (real μ),   E_X(0) = H (Haar / Rademacher random multiplicative model),   E_X(−1) = G_{|μ|} = Σ_k λ_k.
```
- The Codex sparse twist is a single deterministic ε. On [−1, 1] the observed minimum of E_X is at z = 1, for every X tested.

**(I34.4) Exact computation** (`r34_layers.c`):
```
E_X(z) = Σ_{d squarefree} (1 − z²)^{ω(d)} d^{−1} B_d(z) + centering,
B_d(z) = Σ_{m',n' ≤ X/d, (m'n',d)=1} b_{m'} b_{n'} K,   b_n = μ(n) z^{ω(n)} n^{−1/2},
```
- This uses z^{−2ω(g)} = Σ_{d|g} Π_{p|d}(z^{−2} − 1).
- B_d is evaluated as a bilinear Fejér step integral at 48 roots of unity, and E_k comes from the inverse DFT (imaginary residue ≤ 1e-7).
- Runtime: 10 s at X = 10⁷.

**(I34.5) Dirichlet-series form** (stated in §0; PROVED by per-prime algebra).
- Per prime, 1 − z(x+y) + xy = (1−x)^z(1−y)^z(1−xy)^{−(1−z²)}(1−x²)^{−z(1−z)/2}(1−y²)^{−z(1−z)/2}·(1 + O(third order)).
- Exact special cases: z = 1 gives 1/(ζ(a)ζ(b)); z = 0 gives ζ(a+b)/ζ(2a+2b), the Haar diagonal; z = −1 gives ζ(a)ζ(b)/(ζ(2a)ζ(2b)), the |μ| translate.
- **Where the layer masses come from:**
  - ζ(a+b)^{1−z²} is the counted-once diagonal pole. It is cancelled exactly at z = ±1.
  - The X-scale layer masses come from the s = 1 singularity of ζ(s)^{−z}: Selberg–Delange gives Σ_{n≤x} μ(n)z^{ω(n)} = x(log x)^{−z−1}(A(z)/Γ(−z) + O(1/log x)), uniformly for |z| ≤ A.
  - Every term of the expansion carries 1/Γ(−z−j), all of which vanish at z = 1 (classical; Tenenbaum II.5).
  - Main-term model: E_X(z) ≈ L·4X·A(z)²J(z)(log X)^{−2z−2}/Γ(−z)², with A(z) = Π(1−z/p)(1−1/p)^{−z} and J(z) = Π(1 − 2z/p + p^{−2})/(1 − z/p)². This is HEURISTIC; for the bilinear form it is NUMERICAL.
  - Ratio of E_X to the model at X = 10⁵/10⁶/10⁷: z = −1: 0.83/0.86/0.88; z = −0.5: 0.95/0.96/0.96; z = 0.5: 2.9/2.3/2.0; z = 0.95: 24.8/6.9/3.1. It converges as X grows.
- **This is the main term that cancels between the prime class and the rest.** The layers are the z-Taylor coefficients of a function whose main terms all vanish at z = 1. Individually they are ≍ X; their alternating sum G is O(1) in range.

## 3. Why the sparse twist breaks it, and what real μ needs

The twist is not μ_b for any fixed b: it is X-dependent. It violates the inverse identity at the flipped primes (Σ_{d|p}μ′(d) = 2) and nowhere else up to X.

**Ledger at X = 10⁶** (`r34_results.json` → ledger):

| quantity | μ | twist | elementary bound / PNT level |
|---|---|---|---|
| Σ_{n≤X} μ(n)⌊X/n⌋ (Meissel) | 1 | 37,401 | = 1 exactly |
| Σ_{n≤X} μ(n)/n | 0.00020 | 0.0609 | elementary \|·\| ≤ 1 still holds; PNT O(1/L²) = 0.0052 is violated |
| F(1) = Σ μ(n)n^{−1/2} − c | 0.207 | 47.78 | shift = 2Σ p^{−1/2} = 47.576 ≍ √X/log X |
| M(X) | 212 | 37,612 | PNT o(X/log² X) is violated |

**What this means:**
- In layer language the twist moves the Fejér measure from the identity to τ_ε. G_twist ≈ L·F′(1)² regains a √X/log X-sized lump at the base point, which no Selberg–Delange cancellation can remove.
- **The quantitative property of real μ that the balance needs** is sup_{y≤X} |Σ_{n≤y} μ(n)n^{−1/2} − 1/ζ(½)| ≤ X^ε. This is the σ = ½ Littlewood form, since G ≤ 2Q_c + L F² + 2|F||m_c| and Q_c is controlled by the sup.
- The twist test itself is defeated by PNT-level o(√X/log X). A family of twists with Σ_{flip} p^{−1/2} ≍ X^θ would test the property at level X^θ.

## 4. GAP-34 (exact) and its position

**Statement:**
```
GAP-34:  ∀ε>0 ∃B:   Σ_{k≥0} (−1)^k λ_k(X)  ≤  B X^ε   for all X,
         where λ_k(X) = (−1)^k E_k(X) ≥ 0 are the ratio-class layer masses (PROVED nonnegative),
         Σ_k λ_k = G_{|μ|}(X) ≍ X·L (PROVED ≫ X),   λ_0 = H ≤ L(L+1+(1−c)²),  λ_1 = |prime class| = Ω(√X/log X).
Equivalently:  R_c(X) ≤ −H_c(X) + B X^ε   (the sharpest form; R_c ≥ −H_c is automatic).
```

**Position:**
- GAP-34 ⇔ GAP-33 (ε-form) ⇔ RH. This follows from G = 2Q_c + L F² − 2F m_c, Q_c ≤ G, and R33 Thm 33.2.
- **It is weaker than the uniform Chowla requirement of R33.** Uniform-in-shift log-Chowla (Σ_{h<X}|Σ_n μ(n)μ(n+h)w/n| ≪ X^ε) bounds every shift separately and implies GAP-34. GAP-34 needs only the one signed total over all classes/layers; RH is not known to imply uniform Chowla.
- **Target and attack (not "equivalent so OPEN"):**
  1. *Layer by layer*: impossible in principle. The layers are ≍ X and alternate, so any argument bounding |E_k| separately loses a factor X (PROVED by the λ-sum lower bound).
  2. *Through z, the attack built here*: E_X(z) is a polynomial of degree ≤ 2ω_max(X), nonnegative on [−1, 1].
     - Selberg–Delange with the classical zero-free region gives E_X(z) = main(z) + O(X·exp(−c√log X)) uniformly near z = 1, with main(1) = 0 to all orders.
     - This **reproduces unconditionally** G ≪ L·X·δ(X)² (Vinogradov–Korobov δ). It kills the Codex twist (PNT level, A = 2).
  3. *The first unattained milestone*: E_X(1) ≪ X^{1−δ} for one δ > 0. By R33 Thm 33.2(a) this already gives the zero-free strip Re s > 1 − δ/2 (quasi-RH), which is open.
     - The z-route cannot supply it: for non-integer z, ζ^{−z} has branch points at every zero ρ, so the Selberg–Delange error terms near z = 1 *are* the zeros.
     - What the route would need is a **power-saving Selberg–Delange theorem at z = 1** (Σ_{n≤x} μ(n)z^{ω(n)} with error ≪ x^{1−δ} uniformly for |z − 1| ≤ r). That is a zero-free strip statement, so it is the precise thing to attack.
- **CONJECTURE-34 (NUMERICAL support):** min_{z∈[−1,1]} E_X(z) = E_X(1) for all X. It holds at X = 10⁴ … 10⁷ (tested on 81–201-point grids). It would mean real μ minimizes the random-sign-averaged Fejér square. It does not by itself bound E_X(1).

## 5. Controls (X = 3·10⁴, each model with its own μ; `r34_controls.py`, brute-force pairs)

| model | E_0 = H | E_1 (prime layer) | E_2 | max\|E_k\| (k) | Σ\|E_k\| | G = E(1) | G/Σ\|E_k\| | layer signs alternate? | (a) sign check |
|---|---|---|---|---|---|---|---|---|---|
| ζ | 94.32 | −1013.70 | +2.70e4 | 1.27e5 (4) | 3.70e5 | 5.106 | 1.4e-5 | yes (PROVED) | 0 / 4.4e5 |
| ζ² | 2054.4 | −1.91e4 | +4.10e5 | 5.89e6 (5) | 1.96e7 | 65.22 | 3.3e-6 | yes | 0 / 6.1e5 |
| F_0.1 | 2636.1 | −2.35e4 | +6.31e5 | 8.57e6 (5) | 2.85e7 | 85.72 | 3.0e-6 | yes | 0 / 6.1e5 |
| EPS | 6512.8 | −2355.2 | +3.81e4 | 1.21e5 (4) | 4.15e5 | 190.6 | 4.6e-4 | **no** (E(−1) = 209 ≪ Σ\|E_k\|) | **fails**: 36,585 / 1.3e5 |
| DH | 415.3 | −291.5 | −759.2 | 759 (2) | 3337 | 2.977 | 8.9e-4 | **no** (E_2 < 0) | no mismatch found (not proved) |
| Beurling | 7502.5 | −5056.0 | +1.60e4 | 1.60e4 (2) | 4.84e4 | 7034 | 0.15 | yes | holds (multiplicative) |

**What breaks for each control:**
- **ζ², F_0.1 and Beurling** are multiplicative Dirichlet inverses. They show the same alternating X-scale layers and the same payback: the z = 1 exponent is an integer for each of them (ζ² ↔ ζ^{−2z}, F_θ ↔ (ζ(s−θ)ζ(s+θ))^{−z}). **None breaks the payback.** They differ only in the post-payback residual, which is governed by their zeros: G ∝ X^{2θ} for F_θ, ~L³ for ζ² (R33).
- **Beurling's G/Σ|E| = 0.15** is the R32/R33 transient from c = 1/ζ_B(½) = −25.9 (a_1 = 26.9 dominates H). It is not a layer effect.
- **EPS and DH** are non-multiplicative. For them the ratio-class regrouping loses its sign structure ((a) fails for EPS), the layers are not Euler-product Taylor coefficients, and E(−1) is not the |μ| translate. Their G still shows the payback at this X; their failure is again the R33 exponential growth of G at larger X.
- **Conclusion (NUMERICAL + PROVED structure):** the ratio-class/layer balance is a property of *every* Dirichlet inverse, not of ζ. The zeros enter only beyond all Selberg–Delange orders. This is exactly where GAP-34 sits.

## 6. Lean candidates (proposals; claim on the bridge first; NOT built)

- **(L11) Ratio-class regrouping identity** (finite, Finset/Nat.Coprime): for f : ℕ → ℝ and X,
  Σ_{m,n≤X} μ(m)μ(n) f(m/gcd, n/gcd) = Σ_{(u,v) coprime} μ(u)μ(v) f(u,v) Σ_{g≤X/max(u,v), Coprime g (uv)} μ(g)².
- **(L12) Layer sign:** for squarefree m, n, μ(m)μ(n) = (−1)^{ω(m/gcd)+ω(n/gcd)}. Elementary from ArithmeticFunction.moebius.
- Spare **(L13)**: the random-sign identity 𝔼_ε[ε(m)ε(n)] = z^{ω(mΔn)} for independent ±1 coordinates of mean z, giving E_X(z) ≥ 0 on [−1, 1].

## 7. Files

| file | contents |
|---|---|
| `r34_layers.c` | Fast exact layers for real μ / twist; `l_X.json`, `lt_X.json` for X = 10⁴ … 10⁷ |
| `r34_pairs.c` + `r34_controls.py` | Brute-force layers for any model, plus the (a) sign test: `pairs_*.json`, `r34_controls.json`, `r34_controls_summary.json` |
| `r34_analysis.py` → `r34_results.json` | Layer table, E_X(z) curves vs the Selberg–Delange model, prime-class main term, S_q check, twist ledger |
| `r34_figs.py` → **`fig1_ratio_classes.png`** (main) | (a) \|E_k\| vs k; (b) E_X(z) with the 1/Γ(−z)² model; (c) payback scales vs X; (d) controls |

## 8. Labels and credit

- **PROVED:** I34.1–I34.5; the layer sign; E_X ≥ 0 on [−1, 1]; the Ω(√X/log X) prime class; the X-scale compensation; the twist mechanism.
- **CHECKED:** (b) and (d) to all digits; fast = brute; the S_q main term.
- **NUMERICAL:** the layer tables, the Selberg–Delange shape, the controls.
- **HEURISTIC:** the bilinear Selberg–Delange model.
- **CONJECTURE:** CONJECTURE-34.
- **GAP:** GAP-34.
- **OPEN:** RH, the quasi-RH strip, a power-saving Selberg–Delange at z = 1.

**Credit:**
- **Codex** (D2): the reduced-ratio grouping, the prime-ratio sign and Ω bound, and the sparse-twist control, GPPDiscovery2 `8ef8374`. This round's starting point.
- **Turbo** (D4): audits (PR #15 7062427004; the Ω-not-≫ wording; R33 relay).
- Classical: Selberg, Delange, Tenenbaum, Landau, Littlewood, Vinogradov–Korobov, Tao, Matomäki–Radziwiłł, Bohr.
- Project: Alpöge, Claude Code.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
