"""'On the Origin of Spacetime' (Daniel's Drive, GPP folder, May-June 2025): every alpha formula, recomputed.
Daniel Toupin / Golden Physics Project. Seat: Heart of the Cosmos / Grok (GPPDiscovery5). 2026-11-21.
Drafts are NOT copied here; formulas are transcribed from the PDFs / AppendixG.tex (Origin of Spacetime-Newest.zip).
'log' is ambiguous in the drafts; both ln and log10 readings are evaluated.
"""
from mpmath import mp, mpf, sqrt, log, log10, nsum, inf, findroot
import json
mp.dps = 50
phi = (1 + sqrt(5)) / 2
REF = {"CODATA2022": (mpf("137.035999177"), mpf("0.000000021")),
       "Morel2020_Rb": (mpf("137.035999206"), mpf("0.000000011")),
       "Parker2018_Cs": (mpf("137.035999046"), mpf("0.000000027"))}
T = REF["CODATA2022"][0]
rows = []
def add(version, formula, inv_alpha, printed, note=""):
    r = {"version": version, "formula": formula, "printed": printed, "note": note}
    if inv_alpha is None:
        r.update(computed="n/a"); rows.append(r); return
    r["computed_1/alpha"] = mp.nstr(inv_alpha, 15)
    d = inv_alpha - T; r["diff_vs_C22"] = mp.nstr(d, 4)
    r["rel"] = mp.nstr(abs(d) / T, 3)
    for k, (tv, sv) in REF.items(): r["pull_" + k] = mp.nstr((inv_alpha - tv) / sv, 4)
    rows.append(r)

# V1 (May 24, 2025; OnTheOriginOfSpacetime.pdf / ...GG.pdf), App. G.4:  alpha ~ 1/(theta*varphi) at theta = phi
add("May 24 2025 (OnTheOriginOfSpacetime[GG].pdf) G.4", "alpha ~ 1/(theta*phi), theta=phi (reading phi_angle=phi)",
    phi * phi, "'observed value within projection uncertainty'", "only '~'; angle symbol undefined")
add("May 24 2025 G.4 alt reading", "alpha ~ 1/(theta*phi_deg), theta=phi, phi_deg = golden angle 360/phi^2",
    phi * 360 / phi**2, "same", "alternative reading of the undefined angle")

# V2 (May 25, 2025; MPNP (7),(8), Final_Copy, Information_and_Holographic_Spacetime (24)-(27)), App. G.2
def a_may25(th, L):
    return phi**-6 * th**mpf("1.472") / (1 + mpf("1.78") * L(1 / th))
for nm, L in (("ln", log), ("log10", log10)):
    add(f"May 25 2025 (MPNP (8)/Final_Copy) G.2 [{nm}]",
        f"alpha = phi^-6 * theta^1.472/(1+1.78 {nm}(1/theta)), theta=phi", 1 / a_may25(phi, L), "137.036",
        "two fitted decimals 1.472, 1.78")
# V3 (May 27 2025 (12) through v3 June 19 2025; identical in all later full drafts), App. G.3
def S(th):
    return nsum(lambda k: 1 / (phi**k * (1 + (th - phi**k) ** 2)), [2, inf])
s9 = S(phi**9)
add("May 27 - Jun 19 2025 (MPNP (12),(16),(19)-(22),(28), complete, Revised, v3, Newest.zip AppendixG.tex) G.3",
    "1/alpha = sum_{k>=2} kappa/(phi^k [1+(theta-phi^k)^2]), theta=phi^9, kappa=19.56 ('fitted')",
    mpf("19.56") * s9, "137.036", "kappa explicitly 'fitted'")
add("May 27 - Jun 19 2025 G.3 alt reading", "same sum read as alpha (not 1/alpha)", 1 / (mpf("19.56") * s9), "137.036", "alternative reading")
for nm, L in (("ln", log), ("log10", log10)):
    add(f"May 25 2025 G.2 alt reading [{nm}]", "same expression read as 1/alpha", a_may25(phi, L), "137.036", "alternative reading")
# May 20-22 2025 drafts (Information_and_Holographic_Spacetime (3)-(9), titled 'On the Origin of Spacetime(s)'):
add("May 20-22 2025 (Information_and_Holographic_Spacetime (3)-(9))", "alpha(theta) ~ 1/kappa_U(1)(theta)", None,
    "table lists observed 1/137.035999", "no computable formula")
add("May 14 2025 (On The Origin Of Spacetime V2.docx)", "none (qualitative: constants from projection)", None, "-", "")
kneed = {k: mp.nstr(tv / s9, 12) for k, (tv, _) in REF.items()}
# May 27 (11) has empty formula boxes in the PDF ("Using : Matches the observed value to high precision")
add("May 27 2025 (MPNP (11))", "formula boxes empty in PDF", None, "'matches to high precision'")

# Other constants in the May 25 draft (G.3, G.4)
other = {}
for nm, L in (("ln", log), ("log10", log10)):
    th = phi
    x = sqrt(phi**3 / phi**2) * (th ** mpf("1.32") / (1 + mpf("1.65") * L(1 / th))) ** mpf("1.07")
    # readings of the sin^2 formula: X^2/(1+X^2) with X as above, or with the 1.07 applied after the sqrt-factor
    other[f"sin2thetaW May25 [{nm}] X^2/(1+X^2)"] = mp.nstr(x**2 / (1 + x**2), 8)
    y = (sqrt(phi**3 / phi**2) * th ** mpf("1.32") / (1 + mpf("1.65") * L(1 / th))) ** mpf("1.07")
    other[f"sin2thetaW May25 [{nm}] alt Y^2/(1+Y^2)"] = mp.nstr(y**2 / (1 + y**2), 8)
    other[f"alpha_s May25 [{nm}] (1/phi)((1+log(1/theta))/theta)^0.85"] = mp.nstr((1 / phi) * ((1 + L(1 / th)) / th) ** mpf("0.85"), 8)
other["printed: sin2thetaW"] = "0.231 (obs 0.23122 MSbar at M_Z)"
other["printed: alpha_s"] = "~1"

if __name__ == "__main__":
    for r in rows:
        print(json.dumps(r))
    print("S(phi^9) =", mp.nstr(s9, 20))
    print("kappa needed to hit:", kneed)
    # sensitivity: kappa printed to 4 sig figs -> 1/alpha uncertainty from rounding (+-0.005)
    print("1/alpha range for kappa in [19.555,19.565]:", mp.nstr(mpf("19.555") * s9, 10), mp.nstr(mpf("19.565") * s9, 10))
    for k, v in other.items(): print(k, "=", v)
    json.dump({"rows": rows, "S_phi9": mp.nstr(s9, 30), "kappa_needed": kneed, "other": other},
              open("origin_of_spacetime_versions.json", "w"), indent=1)
