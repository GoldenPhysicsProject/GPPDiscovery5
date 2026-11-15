"""Float NUMERICAL estimate of the uniform bound lam_Z on H_73 from sigma_x (route for the rigorous version).
Mass majorant (CHECKED by Cauchy-Schwarz): for g = sum_{j>=J} b_j phi_j, |G(t)|^2 <= (2/l) sin^2(tl) ||g||^2 sum_{j>=J} omega_j^2/(omega_j^2-t^2)^2, t < omega_J."""
import numpy as np, json
from scipy.integrate import quad
l = 0.25; J = 73
def dens(t):
    j = np.arange(J, J+200000); om = j*np.pi/l
    return (2/l)*np.sin(t*l)**2*np.sum(om**2/(om**2-t**2)**2)
M = lambda T: quad(dens, 0, T, limit=400)[0]
r = json.load(open('sigx_results.json'))
smin_hi = r['min_sigx_t>=163.36']; sneg = min(0.0, min(v for t,v in r['sigx_at'] if t <= 31))
M31, M163 = M(31.03), M(163.36)
lam = (smin_hi*(np.pi/2 - M163) + (-0.37)*M31)/np.pi
out = dict(M_0_31=M31, M_0_163=M163, sigx_min_t_ge_163=smin_hi, sigx_floor_0_31_used=-0.37, lamZ_estimate_g_units=lam)
print(json.dumps(out, indent=1)); json.dump(out, open('lamZ_est_results.json','w'), indent=1)
