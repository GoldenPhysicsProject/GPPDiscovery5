"""Round 22, approach 2: Jacobi matrix / Krein string of the zeta circuit (from Round 20's S-fraction, Taylor data only),
and the archimedean obstruction.
(i) Even contraction: alpha_0 = a_1, alpha_n = a_{2n} + a_{2n+1}, beta_n = a_{2n-1} a_{2n}. J = tridiag(sqrt beta, alpha, sqrt beta)
    is the Jacobi matrix of nu (moments c_j); its eigenvalues are Gauss nodes x_k of nu; under RH x = 1/gamma^2.
(ii) Archimedean density of states theta'(t)/pi, theta'(t) = (1/2) Re psi(1/4 + i t/2) - (1/2) log pi: sign and root.
(iii) Archimedean-only Fejer network (pole + Gamma, no primes): sign of rho_X^arch."""
import json, math, numpy as np, mpmath as mp
from r22_fejer import rho
mp.mp.dps = 60
d = json.load(open('../round20_circuit/r20_ladder_ZETA.json'))['ZETA']
a = [mp.mpf(x) for x in d['a']]; n2 = len(a)//2
alpha = [a[0]] + [a[2*n - 1] + a[2*n] for n in range(1, n2)]
beta = [a[2*n - 2]*a[2*n - 1] for n in range(1, n2)]
N = len(beta)
J = mp.matrix(N, N)
for i in range(N):
    J[i, i] = alpha[i]
    if i + 1 < N: J[i, i + 1] = J[i + 1, i] = mp.sqrt(beta[i])
ev = sorted([float(x) for x in mp.eigsy(J, eigvals_only=True)], reverse=True)
gam_nodes = [1/math.sqrt(x) for x in ev if x > 0]
out = dict(alpha=[mp.nstr(x, 12) for x in alpha[:N]], beta=[mp.nstr(x, 12) for x in beta],
           gauss_node_ordinates_1_over_sqrt_x=gam_nodes[:30], min_eigenvalue=min(ev))
# illustration only: compare the first Gauss-node ordinates with zeta zeros
zz = [float(mp.zetazero(k).imag) for k in range(1, 31)]
out['illustration_zero_ordinates'] = zz
out['illustration_rel_err_first25'] = [abs(g - z)/z for g, z in zip(gam_nodes[:25], zz)]
nn = np.arange(1, N); lb = np.log([float(b) for b in beta[1:]])
out['loglog_slope_beta_vs_n_(n>=10)'] = float(np.polyfit(np.log(nn[9:]), lb[9:], 1)[0])
# (ii) archimedean density sign
thp = lambda t: 0.5*mp.re(mp.digamma(0.25 + 0.5j*t)) - 0.5*mp.log(mp.pi)
tstar = mp.findroot(thp, 17.8)
out['theta_prime_root'] = float(tstar); out['theta_prime_at_0'] = float(thp(0)); out['theta_prime_at_14.13'] = float(thp(14.134725))
out['theta_min_value'] = float(mp.siegeltheta(tstar))
# (iii) archimedean-only Fejer network
g0 = np.arange(0, 100.0001, 0.01); arch = {}
for X in (10, 100, 1000, 10000):
    r = rho('ZETA', X, g0, np.zeros(X + 2))
    arch[str(X)] = dict(min=float(r.min()), at=float(g0[np.argmin(r)]), frac_negative=float(np.mean(r < 0)))
out['arch_only_fejer'] = arch
print(json.dumps({k: v for k, v in out.items() if k not in ('alpha', 'beta')}, indent=1)[:3000])
print('alpha[:5]', out['alpha'][:5], 'beta[:5]', out['beta'][:5])
json.dump(out, open('r22_jacobi.json', 'w'), indent=1)
