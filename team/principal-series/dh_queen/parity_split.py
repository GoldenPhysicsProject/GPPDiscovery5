"""Split Queen's Q_DH(f_Q) into its parity parts. f_Q = e^{-i g0 x} sinh(b x) on [-L,L] (half-width L).
Re f_Q = cos(g0 x)sinh(b x) (ODD), Im f_Q = -sin(g0 x)sinh(b x) (EVEN). Q real-symmetric & reflection-invariant =>
Q(f_Q) = Q(Re) + Q(Im) (no cross terms). Autocorrelations in closed form: S = G_O+G_E = Re G_Q (Queen),
D = G_O - G_E = 1/2[ int_0^T cos(g y)cosh(b y) dy - cosh(b u) sin(g T)/g ], T = 2L-u. Same Levy-Khintchine evaluation as dh_sinh.py.
NUMERICAL."""
import numpy as np, math, mpmath as mp, json
exec(open("dh_sinh_queen_copy.py").read().split('g0, b = ')[0])   # coefficient sieve only (c, lam, N)
g0, b = 85.699348, 0.308517
def S(u, L): return math.cos(g0*u)*(math.sinh(b*(2*L-u))/(2*b) - (2*L-u)*math.cosh(b*u)/2) if u < 2*L else 0.0
def D(u, L):
    if u >= 2*L: return 0.0
    T = 2*L-u; z1 = complex(b, g0); z2 = complex(-b, g0)
    I = (((np.exp(z1*T)-1)/z1 + (np.exp(z2*T)-1)/z2)/2).real
    return 0.5*(I - math.cosh(b*u)*math.sin(g0*T)/g0)
def Q(L, Gf, coef, Om0, nu):
    nf = Gf(0.0); pts = [0] + [2*L*j/4000 for j in range(1, 4001)]
    arch = mp.quad(lambda u: (nf - Gf(float(u)))*nu(u), pts) + nf*mp.quad(nu, [2*L, mp.inf])
    X = min(int(math.exp(2*L)), N); pr = 2*sum(coef[k]*k**-0.5*Gf(math.log(k)) for k in range(2, X+1))
    return dict(norm2=nf, Q_over_norm=(Om0*nf + float(arch) - pr)/nf)
Om = float(mp.digamma(0.75) - mp.log(mp.pi/5)); nu = lambda u: mp.e**(-u/2)/mp.sinh(u)
out = []
for L in [1.0, 1.5, 3421/2000, 2.0, 2.5, 3.0, 4.0]:
    E = Q(L, lambda u: (S(u, L)-D(u, L))/2, c, Om, nu); O = Q(L, lambda u: (S(u, L)+D(u, L))/2, c, Om, nu); T = Q(L, lambda u: S(u, L), c, Om, nu)
    r = dict(L=L, total=T['Q_over_norm'], even_ImQ=E['Q_over_norm'], odd_ReQ=O['Q_over_norm'], w_even=E['norm2']/T['norm2'],
             recombined=(E['Q_over_norm']*E['norm2']+O['Q_over_norm']*O['norm2'])/T['norm2'])
    out.append(r); print(json.dumps(r), flush=True)
json.dump(out, open('parity_split.json', 'w'), indent=1)
