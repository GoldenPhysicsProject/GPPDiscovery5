# S16 — n=3 window ½log3 < h < log2: 2D admissible region of (c2,c3) and the floor

W = κ·Pole + Arch − c2 T2 − c3 T3 on C_h, S14 exact Arb entries, Legendre·(1−u²)⁴ basis, Ritz M=24,36, mpmath dps 80.
ζ: κ=1, c2=√2 log2, c3=2log3/√3. F_{1/4}: c_p·cosh(¼ log p) (prime-weight convention of S14; its shifted pole/gamma terms ignored).
Ritz region ⊇ true region (PSD on space ⇒ PSD on subspace), so Ritz edges are OUTER bounds; a point outside them is excluded.

Ritz edges (offsets from ζ), M=36:
| h | sector | c2 up | c2 down | c3 up | c3 down |
|---|---|---|---|---|---|
|0.62|even|+3.0e-8|−6.1e-5|+3.4e-4|−2.8e-2|
|0.62|odd|+4.9e-3|−1.27e-6|+0.117|−1.5e-3|
|0.69|even|+8.3e-11|−5.5e-7|+1.8e-7|−7.6e-5|
|0.69|odd|+7.2e-5|−6.2e-9|+7.2e-4|−1.9e-6|
(M=24→36 roughly halves each edge; not converged.)
- Even sector supplies the CAPS (c2, c3 upper; κ lower); odd sector supplies the FLOORS (c2, c3 lower; κ upper). Same split as n=2 (S14), now for both primes.
- Region is a thin sliver whose binding normal is almost the c2 axis (c3 sensitivity ~5e-4 (even), ~3e-3 (odd) of c2). No product-formula relation c2↔c3 visible: the two sectors' tilts differ.
- κ (pole weight) is itself pinched: h=0.62,M=24: κ−1 ∈ [−9.7e-10, +1.8e-5]; h=0.69,M=36: [−1.9e-12, +5.7e-8].
- Mechanism (PROVED, already in GPPVerify TwoFacePoleSignature): Pole is rank-1 PSD on even, rank-1 NSD on odd. Arch+Pole alone is NOT PSD in either sector (min eig −1e-3 even, −5e-3..−1e-2 odd).
- F_{1/4} is outside the even cap at all h tested (and outside the odd region at h≥0.62). Float, Ritz-outer.
Files: s16_region.py, s16_floor.py, s16_exactP.py, r*.json f*.json k*.json.
