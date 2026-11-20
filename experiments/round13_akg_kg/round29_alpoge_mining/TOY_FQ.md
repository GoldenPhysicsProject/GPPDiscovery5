# Round 29 (B): the √q multiplicity trade-off, run on a toy (Stepanov's method for y² = f(x) over F_p)

**Brief:** Daniel t525u ("All three"), item B. Goal: show in the simplest case the mechanism mined from Alpöge's tree
(NOTES.md §3–§4): **forced divisibility ≤ size** for a nonzero object, a multiplicity amplifier, and error ≈ q/ℓ + ℓ with
optimum **ℓ ≈ √q**.
**Files:** `toy_fq.py` (self-contained, numpy only, about 65 s), `toy_fq.json`, `toy_fq.png`.
**Relation to Round 16:** R16 (`../round16_stepanov/stepanov_ff.py`) ran *Bombieri's* version at a single parameter choice
over F_{p²}. This toy runs *Stepanov's original* hyperelliptic version over prime fields, **sweeps the multiplicity ℓ**, and
records the ASK/MAIN bookkeeping explicitly. Labels: PROVED (classical), CHECKED (computed exactly and verified), NUMERICAL.

## 1. The classical theorem (PROVED, classical)
Stepanov (1969). The version followed here is Iwaniec–Kowalski, *Analytic Number Theory*, §11.6 (Prop. 11.18, Thm 11.13); see
also Bombieri (1973, Sém. Bourbaki 430). Let f ∈ F_q[X] of degree m, not a square, c = (q−1)/2, g = f^c, and
S_a = {x : f(x) = 0 or g(x) = a}. For q > 8m and m < ℓ ≤ q/8 there is a **nonzero**
r = f^ℓ Σ_{j<J} (r_j + s_j g) X^{jq} with deg r_j, deg s_j ≤ c − m, vanishing to order ≥ ℓ at every point of S_a, and
deg r < cℓ + 2mℓ(ℓ−1) + mq. Hence

  **|S_a| ≤ c + 2m(ℓ−1) + mq/ℓ**  (the per-multiplicity cost 2mℓ plus the shared trivial loss mq/ℓ),

and ℓ ≈ √q gives |S_a| < q/2 + 4m√q, so |N − q| < 8m√q for N = #{(x,y): y² = f(x)}.
The lower bound comes from S_{−1}, because X^q − X = X(X^c − 1)(X^c + 1).

## 2. The dictionary to Alpöge's mechanism (the point of this toy)
| Stepanov over F_q | Alpöge over ℤ (NOTES.md) |
|---|---|
| nonzero polynomial r (Lemma 11.19 needs f not a square) | nonzero algebraic integer (period determinant) |
| **ASK** = Σ_{x∈S_a} ord_x(r) ≥ ℓ·\|S_a\| (forced: x^q = x and g(x) = a kill the low Hasse derivatives) | `ASK` = Σ_{p inert} log p · Σ_i (a/p^i − w)_+ (forced: D^p ≡ 0 mod inert p) |
| **MAIN** = deg r | `Main` / `DIV` (archimedean size, Stirling and Chowla–Selberg) |
| ASK ≤ MAIN: zeros of a nonzero polynomial ≤ degree (the product formula on P¹: Σ_v ord_v r = 0, ∞ carries −deg r) | ASK − Main ≤ error: product formula / slope inequality, \|N\| ≥ 1 |
| multiplicity ℓ via Hasse derivatives | m theta lifts (no derivatives on ℤ) |
| shared trivial loss mq/ℓ | x/m (the a_E term) |
| per-multiplicity cost 2mℓ | about 1 per lift (the 2n·m² term) |
| ℓ* ≈ √q | m* ≈ √x |

## 3. Results
**Part 1, explicit construction (CHECKED).** For p ∈ {53, 101}, f₃ = x³+2x+1 (genus 1) and f₅ = x⁵+2x³+x+3 (genus 2),
a = ±1, and every ℓ from 1 to ⌊p/8⌋, the script:
- builds the I–K linear system σ^{(k)} = 0 for k < ℓ;
- takes the **smallest J** with a nonzero solution;
- constructs r explicitly;
- verifies r ≠ 0 and ord_x(r) ≥ ℓ at **every** x ∈ S_a (by repeated exact division by (X − x));
- verifies ASK ≤ MAIN.

**All 72 runs pass** (`vanish_ok`, `ASK_le_MAIN`). Example, p = 101, f₃, a = 1, |S_1| = 46:

| ℓ | J | deg r = MAIN | ASK | ℓ\|S\| | bound ⌊deg r/ℓ⌋ on \|S\| |
|---|---|---|---|---|---|
| 1 | 1 | 153 | 48 | 46 | 153 |
| 3 | 2 | 264 | 140 | 138 | 88 |
| 6 | 4 | 479 | 279 | 276 | 79 |
| 9 | 5 | 621 | 418 | 414 | 69 |
| 12 | 7 | 827 | 558 | 552 | 68 |

The resulting bound on |N_aff − p| for p = 101 (genus 1) falls from 204 (ℓ = 1) to 34 (ℓ = 12). The true value is 10, and
Hasse–Weil gives 2√p ≈ 20.1. At such small p the ℓ ≤ p/8 window ends before the constants settle, so the toy does not beat
Hasse here. That is expected: Stepanov's constant is 8m, not 2g. (CHECKED)

**Part 2, bound vs ℓ for large q (NUMERICAL, exact I–K dimension count).** For q ∈ {1009, 10⁴+7, 10⁵+3, 10⁶+3}, take the
smallest J with A = 2(c−m+1)J > B = Σ_{k<ℓ}(c−m+k(m−1)+J), deg bound (11.24), and E(ℓ) = 2⌊deg/ℓ⌋ − q.
- **The curves E(ℓ)/√q against ℓ/√q collapse onto one curve** for all four q (figure, left panel).
- Least-squares shape (staircase from integer J aside): E ≈ 3.45·q/ℓ + 3.1·ℓ for genus 1, E ≈ 5.45·q/ℓ + 5.2·ℓ for genus 2.
- Optimum: **ℓ* = 0.80–0.89·√q**, with **E(ℓ*)/√q ≈ 4.9** (genus 1) and **≈ 9.0** (genus 2), constant in q.
  That is the O(√q) law. (NUMERICAL)

**Part 3, true counts (CHECKED).** N_aff for f₃ and f₅ at the four large q, by Legendre symbols, all lie within E(ℓ*):

| q | g | N_aff − q | E(ℓ*) | 2g√q |
|---|---|---|---|---|
| 1009 | 1 | +50 | 159 | 63.5 |
| 1009 | 2 | −16 | 297 | 127.1 |
| 10007 | 1 | +20 | 505 | 200.1 |
| 10007 | 2 | −120 | 919 | 400.1 |
| 100003 | 1 | +4 | 1555 | 632.5 |
| 100003 | 2 | +157 | 2845 | 1264.9 |
| 1000003 | 1 | +1196 | 4907 | 2000.0 |
| 1000003 | 2 | −2381 | 8957 | 4000.0 |

## 4. Figure
`toy_fq.png`. Left: E(ℓ)/√q against ℓ/√q for four q (collapse, U-shape, minimum near ℓ ≈ 0.8√q). Middle: bounds from the
explicitly built r for p = 53 and 101, against the true |N−p| and Hasse. Right: MAIN = deg r, ASK = Σ ord_x r, and the
forced ℓ|S| against ℓ for p = 101. The gap MAIN − ASK is the room the nonzero polynomial leaves.

## 5. What this toy says for us (HEURISTIC)
- The mechanism is: **(nonzero object) × (forced high-order divisibility at the special points) × (size bound) ⇒ count ≤ size/ℓ**.
  The trade-off q/ℓ + ℓ is generic: the trivial loss is paid once and the multiplicity cost is paid per order. That is exactly the
  x/m + m of Alpöge's `errorBracketSizeWindows`.
- Over F_q, multiplicity is free: Hasse derivatives exist, and x^q = x makes them algebraic on the points. **Over ℤ there is no
  derivative in the integer itself.** Alpöge's substitute is theta derivatives at CM points, where D^p ≡ 0 mod an inert p.
  See `OUR_INTEGER.md` for what we might use.
- Nothing here concerns W ≥ 0. Positivity of the Weil form for curves is a *consequence* (Weil, Bombieri), never an input.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
