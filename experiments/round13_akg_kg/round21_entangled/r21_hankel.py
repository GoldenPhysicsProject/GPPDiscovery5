"""Round 21, Task 1-2: the geometric thermofield double across the boundary x = 1.
The theta transfer operator (T f)(x) = int_0^oo F_M(x/y) f(y) dy/y  has kernel K(x,y) = F_M(x/y)/y and is diagonal on
primaries: T x^{-s} = M_M(s) x^{-s}.  Its cross-boundary block (y > 1 -> x < 1), identified back to the right half by J,
is in log coordinates the HANKEL operator  (Gamma g)(u) = int_0^oo phi_M(u+v) g(v) dv  on L^2(0,oo).
Psi(x,y) = K(x,y) 1_{x<1<y} is the 'two bulks entangled through x = 1' state.
  Schmidt coefficients = singular values of Gamma (automatically >= 0, no content);
  SIGNS = eigenvalues of the self-adjoint Gamma = the reflected (Osterwalder-Schrader) form <Theta f, T f>.
Nystrom with Gauss-Legendre on [0, L].  Models: ZETA, EPS (x^2+5y^2), DH, F_theta (theta = 0.2)."""
import json, os, sys, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r21_models import *
import mpmath as mp
t0 = time.time()
MODELS = {'ZETA': (phi_even(F_zeta), 70.0), 'EPS': (phi_even(F_eps), 80.0), 'DH': (phi_even(F_dh), 12.0),
          'FTH0.2': (phi_fth(0.2), 105.0)}
out = {}
# sanity: F_theta table reproduces its Mellin value at s = 1/2 (closed form M_zeta(1/2+th) M_zeta(1/2-th))
pf = MODELS['FTH0.2'][0]; uu = np.linspace(0, 105, 200001)
mell = 2 * np.trapezoid(pf(uu), uu)
closed = float(M_zeta(mp.mpf(0.7)) * M_zeta(mp.mpf(0.3)))
out['FTH_check'] = dict(mellin_at_half_grid=mell, closed=closed, rel_err=abs(mell - closed) / abs(closed))
print('FTH check', out['FTH_check'])
for name, (phi, L) in MODELS.items():
    res = {}
    w0 = np.linspace(0, 30, 3001); pv = phi(w0)
    res['phi_sign_on_(0,30)'] = dict(min=float(pv.min()), max=float(pv.max()), all_negative=bool((pv < 0).all()),
                                     all_positive=bool((pv > 0).all()))
    for N in (200, 400):
        x, w = np.polynomial.legendre.leggauss(N); u = (x + 1) * L / 2; w = w * L / 2
        A = np.sqrt(w)[:, None] * phi(u[:, None] + u[None, :]) * np.sqrt(w)[None, :]
        ev = np.linalg.eigvalsh((A + A.T) / 2)
        o = np.argsort(-np.abs(ev)); top = ev[o][:12]
        tol = 1e-13 * np.abs(ev).max()
        res[f'N{N}'] = dict(top12_eigs=[float(v) for v in top],
                            schmidt_top12=[float(abs(v)) for v in top],
                            n_pos=int((ev > tol).sum()), n_neg=int((ev < -tol).sum()),
                            trace=float(ev.sum()), trace_quad=float(np.sum(w * phi(2 * u))),
                            sum_pos=float(ev[ev > 0].sum()), sum_neg=float(ev[ev < 0].sum()))
    out[name] = res
    r = res['N400']
    print(name, 'pos', r['n_pos'], 'neg', r['n_neg'], 'top', ['%.3e' % v for v in r['top12_eigs'][:8]], flush=True)
# rigorous-style witnesses for ZETA: Gamma_phi indefinite.  phi < 0 on (0,oo) gives a negative direction (any f > 0);
# a positive direction: f orthogonal to e^{-u/2} on which Gamma_phi = Gamma_g, g = e^{w/2} theta_+(e^{2w}) > 0.
mp.mp.dps = 25
phiz = lambda w: mp.exp(w / 2) * mp_thplus(mp.exp(2 * w)) - mp.exp(-w / 2) / 2
def G(a, b):   # <e^{-a u}, Gamma_phi e^{-b v}> = int_0^oo phi(W) K_ab(W) dW, K_ab = (e^{-aW}-e^{-bW})/(b-a) or W e^{-aW}
    # phi = g - e^{-W/2}/2 with g = e^{W/2} theta_+(e^{2W}); pole part in closed form, theta part by quadrature on [0,4]
    # (g(W) < e^{-pi e^{8}} beyond W = 4)
    K = (lambda W: W * mp.exp(-a * W)) if a == b else (lambda W: (mp.exp(-a * W) - mp.exp(-b * W)) / (b - a))
    g = lambda W: mp.exp(W / 2) * mp_thplus(mp.exp(2 * W))
    return mp.quad(lambda W: g(W) * K(W), [0, 0.25, 0.5, 1, 2, 4]) - 1 / (2 * (a + mp.mpf(1) / 2) * (b + mp.mpf(1) / 2))
h = mp.mpf(1) / 2; a = mp.mpf(3); c = 1 / (a + h)       # f = e^{-3u} - c e^{-u/2} is orthogonal to e^{-u/2}
neg = G(h, h); pos = G(a, a) - 2 * c * G(a, h) + c * c * G(h, h)
out['ZETA_witness'] = dict(form_on_e_minus_u_over_2=mp.nstr(neg, 15),
                           form_on_e_minus_3u_minus_c_e_minus_u_over_2=mp.nstr(pos, 15),
                           note='two quadratic forms of the SAME Hankel operator Gamma_phi (ZETA) with opposite signs => indefinite (CHECKED, 25 digits)')
print(out['ZETA_witness'])
out['meta'] = dict(seconds=round(time.time() - t0, 1), quadrature='Gauss-Legendre Nystrom on [0,L], N=200/400')
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r21_hankel.json'), 'w'), indent=1)
