# S51 (Grok): audit of Codex finite xi Hankel (Grommer) bounds; detection depth
- `s51_zeros.py` — 3000 Arb-certified zeta zeros (python-flint acb.zeta_zeros) -> zeros_arb_3000.txt
- `s51_audit.py` -> `audit_out.txt` — H51a: bracket certification (Arb + Hardy Z sign change), exact-Fraction Lagrange/Frobenius majorants vs true Arb values, tail constants, zero-free h_n from Phi-kernel moments vs zero sums, det H3, lambda_min(H_N) N=3..10, CERTIFIED H_N > c I for N<=6 (Sylvester in Arb, tail S4 via N(t)<=(t/2pi)log t).
- `s51_depth.py` -> `depth_*.txt` — H51b(2): planted-quartet detection depth N* (first negative LDL pivot). Args: dps mode(zeta|lattice:c) K NMAX targets deltas.
- `s51_negeig.py` -> `negeig_out.txt` — magnitude/count of negative eigenvalues near N*.
- NOTE: dps-300 values for N≳50 and dps-800 values for N≳115 are precision-contaminated (spurious early negative pivots); k=100 values from dps 1600 (not cross-checked higher).
- `depth_truncK_k10.txt` — N* vs number K of retained zeros (Kronecker finite-atom regime N*=K+1, then saturation).
Labels: CERTIFIED (Arb/Fractions, modulo Platt–Trudgian, Rosser, Hadamard), NUMERICAL otherwise. RH OPEN. phi absent.
