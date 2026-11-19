"""Round 18 candidate (e): arithmetic Hardy-space coherent states e_s(n) = n^{-s} (GPPVerify2 ArithmeticHardySpace, z = s - 1/2).
(1) Gram <e_w, e_s> = zeta(s + conj w) samples zeta only on Re >= 1 when both states are in (the closure of) l^2: zero-free region.
(2) Orthogonality catastrophe K_sigma(t,u) = zeta(2 sigma + i(u-t)) / zeta(2 sigma) -> 0 (t != u) as sigma -> 1/2+.
(3) Boundary kernel: zeta(1 + eps + i tau) ~ 1/(eps + i tau) -> pi delta(tau) - i PV(1/tau) + gamma: integral of Re over |tau|<1 -> pi.
(4) Joint eigenvectors of A_p f(n) = f(p n) need complete multiplicativity: Epstein r_Q(n) n^{-s} is not one (ratio r(pn)/r(n) not constant).
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from laws17 import eps_reps
mp.mp.dps = 20
OUT = {}
OUT['catastrophe'] = [dict(sigma=str(sg), K_t_minus_u_1=mp.nstr(abs(mp.zeta(2*mp.mpf(sg) + 1j)/mp.zeta(2*mp.mpf(sg))), 8),
                           vacuum_fidelity_1_over_zeta2sigma=mp.nstr(1/mp.zeta(2*mp.mpf(sg)), 8))
                      for sg in ('0.75', '0.6', '0.51', '0.501', '0.5001')]
OUT['boundary_delta'] = [dict(eps=str(e), int_Re_zeta_1_eps_itau_over_abs_tau_lt_1=mp.nstr(
    mp.quad(lambda t: mp.re(mp.zeta(1 + mp.mpf(e) + 1j*t)), [-1, -0.1, -0.01, 0, 0.01, 0.1, 1]), 8)) for e in ('0.1', '0.01', '0.001')]
OUT['boundary_delta_note'] = "-> pi + 2*gamma-type regular part: Re zeta(1+eps+i tau) -> pi delta(tau) + regular; Im part -> -PV(1/tau). Boundary vectors are delta-normalizable only up to a Hilbert-transform (PV) term: a Hardy (one-sided) structure, not an orthonormal delta-basis."
r = eps_reps(12000)
rows = []
for p in (2, 3, 29):
    rats = sorted({mp.nstr(mp.mpf(r[p*n])/r[n], 4) for n in range(1, 400) if r[n] and p*n <= 12000})
    zero_to_nonzero = [n for n in range(1, 400) if r[n] == 0 and r[p*n] != 0][:5]
    rows.append(dict(p=p, distinct_ratios_on_support=rats, n_with_r_n_0_but_r_pn_nonzero=zero_to_nonzero,
                     eigen=(len(rats) <= 1 and not zero_to_nonzero)))
OUT['eps_not_eigen'] = rows
print(json.dumps(OUT, indent=1))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r18_e_hardy.json'), 'w'), indent=1)
