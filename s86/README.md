# S86 (Grok, 2026-10-22): audit of Codex 0f238bc/38c633a + H86a prime-path expansion of the pole-sewing Schur scalar. RH OPEN.
- audit_codex_threshold.py → .json: exact (sympy, eps=1, {2,3}^3 box) inverse-diagonal = coefficient dual norms; float (K^-1)_11 on {n<=N} vs zeta(2eps) for eps 0.75/0.6/0.5; Mobius packet norm at eps=1/2.
- h86_prime_path.py → h86_prime_path_results.json/.log: parity forms Q = A - P ± poles for ZETA, DH, F_theta (theta 0.1, 0.25; kernel cosh(theta r) W_zeta, two poles); inertia, rho(A^-1 P), Neumann prime-path terms, k-pole Schur matrix.
- h86_inertia_N.py, h86_arch_crossings.py, h86_inertia_count.py: N-robustness of inertia of A vs R; L where arch-only A gains negative directions; inertia vs prime-power count (float64 valid only L<=0.85).
