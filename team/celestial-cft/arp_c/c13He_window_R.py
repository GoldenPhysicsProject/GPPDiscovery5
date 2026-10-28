"""Cross-check of Golden Eye's even Galerkin numbers: min R/||.||_H^2 on even xi in H^1_0(-L, L), L = (a+ell)/2
(archimedean part only, no p = 2 term). Hat Galerkin, M elements: upper bounds (NUMERICAL).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np
from scipy.linalg import toeplitz, eigh
from c14K_gamma_test import toeplitz_P, Om0, a
from c13He_even_bandsplit import sector_basis
for ell in (0.0, 0.028, 0.05, 0.06, 0.07, 0.1):
    L = (a+ell)/2; row = []
    for M in (80, 160, 320):
        h = 2*L/M; n = M-1; col, _ = toeplitz_P(M, h, 14.0, nGL=8)
        hc = np.zeros(M+1); hc[0] = 2/h + h/6; hc[1] = -1/h + h/24
        Rm = toeplitz(col[:n]) + Om0*toeplitz(hc[:n]); H = toeplitz(hc[:n])
        out = []
        for s in ('even', 'odd'):
            B = sector_basis(n, s); out.append(eigh(B.T@Rm@B, B.T@H@B, eigvals_only=True)[0])
        row.append(f"M={M}: even {out[0]:+.5f}, odd {out[1]:+.5f}")
    print(f"ell={ell:.3f}: " + "; ".join(row), flush=True)
