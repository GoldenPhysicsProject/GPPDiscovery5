"""Named candidate formulas for 1/alpha, evaluated at 50 digits (mpmath).
Daniel Toupin / Golden Physics Project. Seat: Heart of the Cosmos / Grok (GPPDiscovery5).
Reference values (all low-energy, Thomson limit q->0):
  CODATA 2022      137.035999177(21)   (NIST SP 961, 2024)
  Fan+ 2023 (g-2)  137.035999166(15)   PRL 130, 071801
  Morel+ 2020 (Rb) 137.035999206(11)   Nature 588, 61
  Parker+ 2018(Cs) 137.035999046(27)   Science 360, 191
"""
from mpmath import mp, mpf, sqrt, pi, gamma, factorial
import json
mp.dps = 50
phi = (1 + sqrt(5)) / 2
REF = {
    "CODATA2022": (mpf("137.035999177"), mpf("0.000000021")),
    "Fan2023_ge": (mpf("137.035999166"), mpf("0.000000015")),
    "Morel2020_Rb": (mpf("137.035999206"), mpf("0.000000011")),
    "Parker2018_Cs": (mpf("137.035999046"), mpf("0.000000027")),
}
C = {}
C["Eddington 1929: 136"] = mpf(136)
C["Eddington 1930s: 137"] = mpf(137)
C["golden angle 360/phi^2"] = 360 / phi**2
C["Heyrovska: 360/phi^2 - 2/phi^3"] = 360 / phi**2 - 2 / phi**3
C["Sherbon/Pellis: 360/phi^2 - 2/phi^3 + (3phi)^-5"] = 360 / phi**2 - 2 / phi**3 + (3 * phi) ** -5
# Wyler 1969: alpha = 9/(8 pi^4) * (pi^5/(2^4 5!))^(1/4)
C["Wyler 1969"] = 1 / (mpf(9) / (8 * pi**4) * (pi**5 / (2**4 * factorial(5))) ** (mpf(1) / 4))
C["4pi^3 + pi^2 + pi (folklore)"] = 4 * pi**3 + pi**2 + pi
# GGC paper III, App. B (Daniel's Drive, ggc_iii_sm_cosmology.tex, eq. B2)
lam = phi**-4
c1, c2, c3 = pi / 3 * lam**3, pi / 8 * lam**4, pi / 668 * lam**5
C["GGC bare 360/phi^2"] = 360 / phi**2
C["GGC 1 corr: 360/phi^2 (1 - pi/3 phi^-12)"] = 360 / phi**2 * (1 - c1)
C["GGC 2 corr: ... - pi/8 phi^-16"] = 360 / phi**2 * (1 - c1 - c2)
C["GGC 3 corr (as printed, claims 137.03599908413)"] = 360 / phi**2 * (1 - c1 - c2 - c3)
C["GGC 3 corr with 300 in place of 668"] = 360 / phi**2 * (1 - c1 - c2 - pi / 300 * lam**5)
C["GGC Paper II: 4pi/phi^2"] = 4 * pi / phi**2
# Singh 2022 (EPJ Plus 137:664) quoted value before the Karolyhazy correction: alpha=0.00729713
C["Singh 2022 (quoted, uncorrected)"] = 1 / mpf("0.00729713")

def row(name, v):
    t, s = REF["CODATA2022"]
    d = v - t
    rel = abs(d) / t
    digits = -mp.log10(rel) if rel > 0 else mpf(99)
    out = {"name": name, "value": mp.nstr(v, 16), "diff_vs_CODATA22": mp.nstr(d, 4),
           "rel": mp.nstr(rel, 3), "sig_digits": float(digits)}
    for k, (tv, sv) in REF.items():
        out["pull_" + k] = float((v - tv) / sv)
    return out

if __name__ == "__main__":
    rows = [row(k, v) for k, v in C.items()]
    for r in rows:
        print(f"{r['name']:55s} {r['value']:>20s}  d={r['diff_vs_CODATA22']:>11s}  digits={r['sig_digits']:5.2f}  "
              f"pull(C22)={r['pull_CODATA2022']:+.3g}  pull(Fan)={r['pull_Fan2023_ge']:+.3g}  "
              f"pull(Rb)={r['pull_Morel2020_Rb']:+.3g}  pull(Cs)={r['pull_Parker2018_Cs']:+.3g}")
    # GGC: what third denominator would reproduce each reference?
    for k, (tv, _) in REF.items():
        need = pi * lam**5 / (1 - c1 - c2 - tv * phi**2 / 360)
        print(f"GGC: denominator n3 needed to hit {k}: {mp.nstr(need, 10)}")
    json.dump(rows, open("candidates.json", "w"), indent=1)
