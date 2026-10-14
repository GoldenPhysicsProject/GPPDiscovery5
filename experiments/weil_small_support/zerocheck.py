"""Normalization check: Q(g) from weil_ritz.form vs direct sum over zeta zeros (first 300, tail ~ negligible for smooth g)."""
import sys, numpy as np
sys.argv = ["x", "14", "0.5", "0.5", "0.1"]
import importlib.util
spec = importlib.util.spec_from_file_location("w", "weil_ritz.py")
src = open("weil_ritz.py").read().split("cases = [")[0]
ns = {"__file__": "weil_ritz.py"}; exec(src, ns)
from mpmath import zetazero
h = 0.5; Q, G = ns["form"](h, "Z", 0.0)
c = np.zeros(14); c[0] = 1; c[2] = 0.3
xq, wq = ns["xq0"]*h, ns["wq0"]*h; B = ns["basis"](xq, h); g = c @ B
gam = np.array([float(zetazero(k).imag) for k in range(1, 301)])
gh = (g*wq) @ np.exp(1j*np.outer(xq, gam))
print("form:", c @ Q @ c, " zero-sum (300 zeros):", 2*np.sum(np.abs(gh)**2))
