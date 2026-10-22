"""Grok S87: Weyl/Landau-Slepian count of negative directions of the archimedean multiplier
Omega(t)=Re psi(1/4+it/2)-log pi on test functions supported in [-L,L]: N_W(L)=(2L)(2 t0)/(2 pi), per parity ~ N_W/2."""
import mpmath as mp, json
t0 = mp.findroot(lambda t: mp.re(mp.digamma(0.25 + 0.5j * t)) - mp.log(mp.pi), 6)
print("t0 =", t0, " Omega(0) =", mp.digamma(0.25) - mp.log(mp.pi))
out = {"t0": float(t0), "per_parity_slope": float(t0 / mp.pi),
       "table": {f"{L:.2f}": float(L * t0 / mp.pi) for L in [0.3 + 0.05 * i for i in range(23)]}}
json.dump(out, open("weyl.json", "w"), indent=1); print(out["per_parity_slope"])
