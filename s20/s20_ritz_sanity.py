"""S20 [D] NUMERICAL sanity: Ritz lambda_min of the S17 matrices (M=12) must exceed the certified c."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "s17")); sys.path.insert(0, os.path.join(HERE, "..", "s14"))
from s17_build import build, forms, arb
import numpy as np, scipy.linalg as la
M = 12; h = arb(3).log()/2; B = build(M, h, prec=256); Fm = forms(B, h)
N = np.array([[float(B["N"][i, j].mid()) for j in range(M)] for i in range(M)])
for nm in Fm:
    Q = np.array([[float(Fm[nm][i, j].mid()) for j in range(M)] for i in range(M)])
    print(f"{nm:32s} Ritz lambda_min(Q,N) = {la.eigh(Q, N, eigvals_only=True)[0]:.6f}")
