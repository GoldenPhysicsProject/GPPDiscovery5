# S87 (Grok, 2026-10-22): H86b test, arch-only inertia vs prime-power counting. RH OPEN.
- hp_forms.py: high-precision rebuild (mpmath, dps 90/110) of the S86 parity forms A (arch+delta), P (primes), pole vectors for ZETA/DH/F_theta. Convolution matrices are exact polynomials (monomial basis, exact integer change of basis to Legendre) contracted with mp Gauss-Legendre kernel moments. Not interval arithmetic.
- t_compare.py: agreement with the float64 s86 builder to ≤3e-13 (float quadrature level) at L=0.6, N=12, all forms.
- scan.py -> scan.json/scan.log: n_-(A), n_-(R=A-P), low eigenvalues, both parities, L=0.30..1.40 step 0.05; ZETA N=24 and 32 (identical counts), DH/F_0.1/F_0.25 N=24.
- weyl.py -> weyl.json: Omega(t)=Re psi(1/4+it/2)-log pi, root t0=6.28984; Landau-Slepian per-parity count L*t0/pi = 2.0021 L.
