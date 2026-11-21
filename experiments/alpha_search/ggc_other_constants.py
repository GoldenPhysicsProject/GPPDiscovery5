"""Check the other closed-form constant claims in GGC III (Daniel's Drive, ggc_iii_sm_cosmology.tex abstract)."""
from mpmath import mp, mpf, sqrt, pi
mp.dps = 30
phi = (1 + sqrt(5)) / 2
MP_full = mpf("1.220890e19")    # GeV, Planck mass sqrt(hbar c/G) (CODATA 2022)
MP_red = MP_full / sqrt(8 * pi)  # reduced Planck mass
v_obs, mH_obs = mpf("246.22"), mpf("125.20")
print("phi^-80 =", mp.nstr(phi ** -80, 6))
print("v = M_P phi^-80   (full M_P):", mp.nstr(MP_full * phi ** -80, 6), "GeV  vs 246.22")
print("v = M_P phi^-80   (reduced) :", mp.nstr(MP_red * phi ** -80, 6), "GeV")
print("m_H = v phi^-4 with v=246.22:", mp.nstr(v_obs * phi ** -4, 6), "GeV  vs 125.20")
print("eta_B = phi^-8/(2pi):", mp.nstr(phi ** -8 / (2 * pi), 6), " vs ~6.1e-10 (Planck/BBN) [GGC text says 2.53e-10]")
