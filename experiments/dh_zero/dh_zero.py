"""Locate the Davenport-Heilbronn off-line zero near 0.8085+85.70i. NUMERICAL (float, mpmath 30 dps).
Also checks the functional equation Lambda(s)=Lambda(1-s) at a sample point.
Run: python experiments/dh_zero/dh_zero.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from mpmath import mp, mpc, findroot
from common import dh, dh_completed, write_result
mp.dps = 30
z = findroot(dh, mpc(0.8085, 85.70))
s0 = mpc(0.3, 7.1)
fe = abs(dh_completed(s0) - dh_completed(1 - s0))
print("zero", z, "|f(z)|", abs(dh(z)), "FE residual", fe)
write_result(os.path.join(os.path.dirname(__file__), "result.json"),
             {"experiment": "dh_zero", "label": "NUMERICAL", "zero": str(z), "abs_f_at_zero": str(abs(dh(z))),
              "fe_residual_at_0.3+7.1i": str(fe)}, kind="float", method="mpmath findroot (secant) on dirichlet() combination")
