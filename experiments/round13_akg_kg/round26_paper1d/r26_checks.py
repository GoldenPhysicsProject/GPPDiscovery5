"""Round 26 sanity checks for paper.tex (CHECKED-level identities; no zeta zeros used).
— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5)."""
import json, mpmath as mp
mp.mp.dps = 40
out = {}
# (1) 1D cut factor on the unitary line: Gamma(1/2+il)Gamma(1/2-il) = pi/cosh(pi l)
out["gamma_1d"] = max(abs(mp.gamma(0.5+1j*l)*mp.gamma(0.5-1j*l) - mp.pi/mp.cosh(mp.pi*l)) for l in [0,0.3,1,2.5,7])
# (2) 2D cut factor: Gamma(1+il)Gamma(1-il) = pi l / sinh(pi l)
out["gamma_2d"] = max(abs(mp.gamma(1+1j*l)*mp.gamma(1-1j*l) - mp.pi*l/mp.sinh(mp.pi*l)) for l in [0.3,1,2.5,7])
# (3) Gamma(D)Gamma(1-D) real iff Re D = 1/2 or Im D = 0 (sample off-line point)
D = mp.mpc(0.7, 0.5); out["offline_imag_part"] = float(mp.im(mp.gamma(D)*mp.gamma(1-D)))
# (4) Theta reflection: T(x) = sum_{n in Z} exp(-pi n^2 x^2) satisfies T(x) = x^{-1} T(1/x)
T = lambda x: mp.jtheta(3, 0, mp.exp(-mp.pi*x**2))
out["theta_reflection"] = max(abs(T(x) - T(1/x)/x) for x in [0.2,0.5,1.3,3.0])
# (5) Mellin of omega(x)=sum_{n>=1} exp(-pi n^2 x^2): int omega x^{s-1} dx = (1/2) pi^{-s/2} Gamma(s/2) zeta(s), at s=3
s = 3
om = lambda x: (T(x)-1)/2 if x>0.3 else ((T(1/x)/x)-1)/2
lhs = mp.quad(lambda x: om(x)*x**(s-1), [0,1,mp.inf])
out["mellin_theta"] = float(abs(lhs - 0.5*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)))
# (6) eta(sigma) > 0 and zeta(sigma) < 0 on (0,1)
out["min_eta_on_(0,1)"] = float(min(mp.altzeta(x) for x in mp.linspace(0.01,0.99,99)))
out["max_zeta_on_(0,1)"] = float(max(mp.zeta(x) for x in mp.linspace(0.01,0.99,99)))
# (7) Poisson kernel: (1-1/p)/|1-p^{-1/2-it}|^2 = 1 + 2 sum_k p^{-k/2} cos(k t log p)
p, t = 3, 1.7
K = (1-1/mp.mpf(p))/abs(1-mp.power(p, -0.5-1j*t))**2
S = 1 + 2*mp.nsum(lambda k: mp.power(p, -k/2)*mp.cos(k*t*mp.log(p)), [1, mp.inf])
out["poisson_kernel"] = float(abs(K-S))
# (8) xi real on the line, xi(1-conj s) = conj xi(s)
xi = lambda s: 0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
out["xi_imag_on_line_max_rel"] = float(max(abs(mp.im(xi(0.5+1j*tt)))/abs(xi(0.5+1j*tt)) for tt in [3,10,20.5,33.3]))
s0 = mp.mpc(0.8, 12.3); out["xi_mirror"] = float(abs(xi(1-mp.conj(s0)) - mp.conj(xi(s0))))
# (9) Weighted reflection unitary on L^2(x^a dx): ||Theta_a f|| = ||f||, f = x e^{-x}
for a in [0, 1]:
    f = lambda x: x*mp.exp(-x)
    n1 = mp.quad(lambda x: abs(f(x))**2*x**a, [0, mp.inf])
    n2 = mp.quad(lambda x: abs(x**(-(1+a))*f(1/x))**2*x**a, [0, 1, mp.inf])
    out[f"theta_unitary_alpha{a}"] = float(abs(n1-n2))
for k,v in out.items(): out[k] = float(v)
print(json.dumps(out, indent=1))
json.dump(out, open("r26_checks.json","w"), indent=1)
