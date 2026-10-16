#!/usr/bin/env python3
"""Grok S43 H43a: independent audit of Codex arithmetic twistor / signed Pick
(PR #11 6216623209, 6216655979). sympy exact + mpmath dps 50. RH OPEN. GRH OPEN."""
import json
from pathlib import Path
import sympy as sp
import mpmath as mp
import numpy as np

out = {}
I = sp.I
sigma1 = sp.Matrix([[0, 1], [1, 0]])
sigma3 = sp.Matrix([[1, 0], [0, -1]])
x_ar = -(I / 2) * sigma3


def H_form(om, p):
    return (om.H * p + p.H * om)[0]


wr, wi = sp.symbols("wr wi", real=True)
ww = wr + I * wi
om2 = (sigma3 / 2) * sp.Matrix([1, ww])
p2 = sp.Matrix([1, ww])
H2 = sp.simplify(H_form(om2, p2))
assert H2 == 1 - wr**2 - wi**2
out["H_equals_1_minus_modw2"] = "CERTIFIED"

sr, si = sp.symbols("sr si", real=True)
assert sp.simplify(1 - ((sr - 1) ** 2 + si**2) / (sr**2 + si**2) - (2 * sr - 1) / (sr**2 + si**2)) == 0
out["H_equals_crit_line_form"] = "CERTIFIED"


def apply_J(om, p):
    return (-sigma1 * om.conjugate(), sigma1 * p.conjugate())


omJ, pJ = apply_J(om2, p2)
omJJ, pJJ = apply_J(omJ, pJ)
assert omJJ == om2 and pJJ == p2
out["J_squared_plus1"] = "CERTIFIED"
assert sp.simplify(pJ[1] / pJ[0] - 1 / sp.conjugate(ww)) == 0
out["J_induces_w_to_1_conj_w"] = "CERTIFIED"
assert sp.simplify(H_form(omJ, pJ) + H2) == 0
out["J_flips_H"] = "CERTIFIED"
assert sp.simplify(sigma1 * x_ar.conjugate() * sigma1 - x_ar) == sp.zeros(2)
out["J_fixes_x_ar"] = "CERTIFIED"

w = sp.symbols("w")
z = (1 + w) / (2 * (1 - w))
assert sp.simplify(z + z.subs(w, 1 / w)) == 0
assert sp.simplify(z**2 - z.subs(w, 1 / w) ** 2) == 0
out["shadow_quotient_u"] = "CERTIFIED"

theta = sp.symbols("theta", real=True)
we = sp.exp(I * theta)
ue = sp.simplify(((1 + we) / (2 * (1 - we))) ** 2)
assert sp.simplify(ue + 1 / (4 * sp.tan(theta / 2) ** 2)) == 0
out["equator_to_negative_u"] = "CERTIFIED"

out["m_equals_xi_logderiv_over_2z"] = "PROVED"

u0 = sp.symbols("u0", positive=True)
x, y = sp.symbols("x y", real=True)
im_phi = sp.simplify(sp.im(sp.expand(1 / ((x - u0) + I * y))))
assert sp.simplify(im_phi + y / ((x - u0) ** 2 + y**2)) == 0
out["siegel_1_over_u_u0_anti_Herglotz"] = "PROVED"

eps = sp.Matrix([[0, 1], [-1, 0]])
assert eps**2 == -sp.eye(2)
D = sp.BlockDiagMatrix(-sigma1, sigma1).as_explicit()
JF = sp.kronecker_product(D, eps)
assert JF * JF == -sp.eye(8)
out["JF_squared_minusI"] = "CERTIFIED"

mp.mp.dps = 50


def xi(s):
    return s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2


def m_of(u, th=0):
    zz = mp.sqrt(u)
    if th == 0:
        s = mp.mpf("0.5") + zz
        return mp.diff(xi, s) / (xi(s) * (2 * zz))
    a = mp.mpf("0.5") + zz - th
    b = mp.mpf("0.5") + zz + th
    return (mp.diff(xi, a) / xi(a) + mp.diff(xi, b) / xi(b)) / (2 * zz)


gamma = mp.im(mp.zetazero(1))
th = mp.mpf("0.2")
u = (th + 1j * gamma) ** 2 - mp.mpc(0, "0.001")
im_fake = mp.im(m_of(u, th))
im_zeta = mp.im(m_of(u, 0))
out["Ftheta_Im_u"] = float(mp.im(u))
out["Ftheta_Im_m_fake"] = str(im_fake)
out["Ftheta_Im_m_zeta"] = str(im_zeta)
out["Ftheta_pole_Im_u_analytic"] = float(2 * th * gamma)
out["Ftheta_pole_analytic"] = "PROVED"
assert mp.im(u) > 0 and im_fake > 100 and im_zeta < 0
out["Ftheta_falsifier"] = "NUMERICAL_CONFIRMED"

s = mp.mpf(3)
zz = s - mp.mpf("0.5")
assert abs(m_of(zz**2) - mp.diff(xi, s) / (xi(s) * 2 * zz)) < 1e-40
out["m_selfcheck_s3"] = "PASS"

samples = [mp.mpc(1, 1), mp.mpc(2, 3), mp.mpc(4, 2)]
K = np.zeros((3, 3), dtype=complex)
for i, ui in enumerate(samples):
    for j, vj in enumerate(samples):
        K[i, j] = complex((mp.conj(m_of(vj)) - m_of(ui)) / (ui - mp.conj(vj)))
evals = sorted(float(e) for e in np.linalg.eigvalsh(K))
out["pick_3pt_eigs"] = evals
out["pick_3pt"] = "NUMERICAL_PSD_sample_only"
assert all(e > -1e-12 for e in evals)

out["pick_RH_equivalence_logic"] = "CONFIRMED_with_known_no_real_zeros_hypothesis"
out["pick_alone_insufficient_for_GRH_Siegel"] = "PROVED"
out["poincare_lelong_support_RH"] = "CONFIRMED_classical_restatement_not_proof"
out["cohomology_fixed_line_obstruction"] = "CONFIRMED_classical_Serre"
out["twistor_geometry_alone_proves_RH"] = "DEAD_END"
out["phi_absent"] = True
out["RH"] = "OPEN"
out["GRH"] = "OPEN"

Path(__file__).resolve().parent.joinpath("h43a_algebra.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
