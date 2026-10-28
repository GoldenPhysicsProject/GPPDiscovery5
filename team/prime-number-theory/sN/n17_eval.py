# N.17 evaluations (float, eval./NUMERICAL; NOT a certificate, no Ritz scan). Odd sector, 2L = log 3.
# (1) sigma_2(t) = Omega(t) - 2 c2 cos(a t) at the N.16 near-null spectral peaks.
# (2) mu_perp(N): bathtub (S&M Lemma B) for the symbol sigma_2 through its increasing minorant
#     w(t) = inf_{s>=t} sigma_2(s), against the projected Lemma A density dens_N = dens_o - sum_{n<=N} (modes of V_N),
#     V_N = first N odd Dirichlet modes e_n = sin(k_n (x+L))/sqrt(L), k_n = n pi/(2L), n = 2,4,...,2N.
import numpy as np
from scipy.special import psi
a = np.log(2); c2 = a/np.sqrt(2); L = np.log(3)/2
Om = lambda t: np.real(psi(0.25 + 0.5j*t)) - np.log(np.pi)
s2 = lambda t: Om(t) - 2*c2*np.cos(a*t)
print("sigma_2 at near-null peaks: odd t=7.54 ->", round(float(s2(7.54)), 4), "; even t=5.46 ->", round(float(s2(5.46)), 4),
      "; min over t in [0,20]:", round(float(s2(np.linspace(0, 20, 200001)).min()), 4), "; sigma_2(0) =", round(float(s2(0.0)), 4))
def dens_o(t):
    return (L - L*np.sinc(2*t*L/np.pi) - 2*np.sin(t*L)*(0.5/np.tanh(L/2)*np.sin(t*L) - t*np.cos(t*L))/(t*t + 0.25))/np.pi
def mode_dens(t, N):
    out = np.zeros_like(t)
    for n in range(2, 2*N + 1, 2):
        k = n*np.pi/(2*L); d = k*k - t*t; d = np.where(np.abs(d) < 1e-12, 1e-12, d)
        eh2 = (2*k*np.sin(t*L)/d)**2/L          # |e_n^(t)|^2 (L2-normalised odd Dirichlet mode)
        out += (0.25 + t*t)*eh2/(k*k + 0.25)/np.pi
    return out
TM = 6000.0; dt = 0.005; t = np.arange(dt/2, TM, dt)
w = np.minimum.accumulate(s2(t)[::-1])[::-1]     # inf_{s>=t} on the grid (sigma_2 >= Om(TM)-2c2 beyond, increasing)
base = dens_o(t)
for N in [8, 16, 32, 64, 128]:
    dN = np.maximum(base - mode_dens(t, N), 0.0)
    D = np.cumsum(dN)*dt
    mu = w[0] + np.sum(np.diff(w)*np.maximum(0.0, 1 - D[1:]))
    print(f"N={N:4d} k_N={2*N*np.pi/(2*L):8.2f}  D_N(TM)={D[-1]:.1f} (~ L*TM/pi - N modes, sanity)  mu_perp(N) ~ {mu:.4f}  (Omega(k_N)-2c2 = {Om(2*N*np.pi/(2*L)) - 2*c2:.4f})")
