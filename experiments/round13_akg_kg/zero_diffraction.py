"""Round 13j(4): 'diffraction' of the zeta zero set (test data: zeros_14.dat, 4520 zeros gamma <= 4999).
D(x) = sum_gamma e^{-(gamma/T)^2/2} cos(gamma x). Guinand-Weil: peaks at x = +-log p^k with height
-(Lambda(n)/ (2 sqrt n)) * T/sqrt(2pi) * ... (sign negative), on a smooth background from the zero density.
Test: peak (in units of the background rms) at log 2, log 3, log 4, log 5, vs log phi, 2 log phi, log(phi^2) = log(phi+1), log sqrt5."""
import numpy as np, json, math
from akg_kg import platt
gam = np.array(platt('/workspace/zeros/zeros_14.dat', 6000)); gam = gam[gam <= 4999]
T = 1500.0; w = np.exp(-(gam/T)**2/2)
xs = np.linspace(0.3, 2.5, 22001)
D = np.array([np.sum(w*np.cos(gam*x)) for x in xs])
# predicted peak: -(Lambda(n) n^-1/2 / 2) * (T / sqrt(2 pi)) * 2? measure it, compare ratios to Lambda(n) n^-1/2
bg = D.copy()
phi = (1+5**0.5)/2
pts = {'log2': math.log(2), 'log3': math.log(3), 'log4': math.log(4), 'log5': math.log(5), 'log7': math.log(7),
       'log_phi': math.log(phi), '2log_phi': 2*math.log(phi), '3log_phi': 3*math.log(phi), 'log_sqrt5': 0.5*math.log(5)}
def val(x): return float(np.interp(x, xs, D))
mask = np.ones_like(xs, bool)
for n in range(2, 13):
    mask &= np.abs(xs - math.log(n)) > 0.01
rms = float(np.sqrt(np.mean(D[mask]**2)))
out = dict(T=T, n_zeros=int(len(gam)), background_rms=rms,
           values={k: val(v) for k, v in pts.items()}, in_rms_units={k: val(v)/rms for k, v in pts.items()},
           predicted_ratio_log3_over_log2=float((math.log(3)/3**0.5)/(math.log(2)/2**0.5)),
           measured_ratio_log3_over_log2=val(math.log(3))/val(math.log(2)))
print(json.dumps(out, indent=1)); json.dump(out, open('zero_diffraction.json', 'w'), indent=1)
