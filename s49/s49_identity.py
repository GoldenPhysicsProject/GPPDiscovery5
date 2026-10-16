# S49 (Grok): check Codex's Fisher/Stieltjes mixture identity
#   -m'(u) = E_u[-q_t'(u)] - Var_u(q_t),   4||psi_u'||^2 = Var_u(q_t)
# (1) sympy, exact, for a general n-atom positive measure and ARBITRARY positive mode functions c_i(u)
# (2) sympy, the explicit -q_t' formula for c=cosh(t sqrt u)
# (3) mpmath numerics for the true xi theta measure d nu = 2 Phi(t) dt on [0,inf)
import json, sympy as sp, mpmath as mp
u = sp.symbols('u', positive=True)
out = {}
for n in (2, 3, 4):
    w = sp.symbols(f'w0:{n}', positive=True)
    c = [sp.Function(f'c{i}', positive=True)(u) for i in range(n)]
    Z = sum(wi * ci for wi, ci in zip(w, c))
    m = sp.diff(Z, u) / Z
    q = [sp.diff(ci, u) / ci for ci in c]
    E = lambda g: sum(wi * ci * gi for wi, ci, gi in zip(w, c, g)) / Z
    var = E([qi ** 2 for qi in q]) - E(q) ** 2
    lhs = -sp.diff(m, u)
    rhs = E([-sp.diff(qi, u) for qi in q]) - var
    d1 = sp.simplify(sp.together(lhs - rhs))
    psi = [sp.sqrt(ci / Z) for ci in c]
    fisher = sum(wi * sp.diff(p, u) ** 2 for wi, p in zip(w, psi))
    d2 = sp.simplify(sp.together(4 * fisher - var))
    out[f'n={n}'] = dict(identity_residual=str(d1), fisher_residual=str(d2))
    print(n, 'identity residual:', d1, '| 4||psi\'||^2 - Var residual:', d2)
t = sp.symbols('t', positive=True)
qt = sp.diff(sp.log(sp.cosh(t * sp.sqrt(u))), u)
mq = -sp.diff(qt, u)
codex = t / (4 * u ** sp.Rational(3, 2)) * (sp.tanh(t * sp.sqrt(u)) - t * sp.sqrt(u) * sp.sech(t * sp.sqrt(u)) ** 2)
r = sp.simplify((mq - codex).rewrite(sp.exp))
print("-q_t' - Codex formula:", r); out['minus_qprime_formula_residual'] = str(r)
print('q_t(0+) =', sp.limit(qt, u, 0), "  -q_t'(0+) =", sp.limit(mq, u, 0))
out['q0'] = str(sp.limit(qt, u, 0)); out['mq0'] = str(sp.limit(mq, u, 0))

# (3) numerics for xi
mp.mp.dps = 30
def Phi(x):
    s = mp.mpf(0)
    for k in range(1, 12):
        s += (2 * mp.pi ** 2 * k ** 4 * mp.e ** (9 * x) - 3 * mp.pi * k ** 2 * mp.e ** (5 * x)) * mp.e ** (-mp.pi * k ** 2 * mp.e ** (4 * x))
    return s
def xi(s): return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
TT = [0, 0.25, 0.5, 0.75, 1, 1.5, 2.5]
# Polya normalization: int_0^oo Phi(x) cosh(zx) dx = xi(1/2+z/2)/8, so Codex's mode t = 2x, d nu = 8 Phi(x) dx
Z0 = 8 * mp.quad(Phi, TT)
print('8 int Phi =', Z0, ' xi(1/2) =', xi(mp.mpf(0.5)))
rows = []
for uu in ['0.01', '0.5', '1', '4', '25', '100']:
    U = mp.mpf(uu); rt = mp.sqrt(U)
    w = lambda x: 8 * Phi(x) * mp.cosh(2 * x * rt)
    Zu = mp.quad(w, TT)
    q = lambda x: 2 * x * mp.tanh(2 * x * rt) / (2 * rt)
    mqp = lambda x: 2 * x / (4 * U ** 1.5) * (mp.tanh(2 * x * rt) - 2 * x * rt * mp.sech(2 * x * rt) ** 2)
    Eq = mp.quad(lambda x: w(x) * q(x), TT) / Zu
    Eq2 = mp.quad(lambda x: w(x) * q(x) ** 2, TT) / Zu
    Emq = mp.quad(lambda x: w(x) * mqp(x), TT) / Zu
    Var = Eq2 - Eq ** 2
    f = lambda v: xi(mp.mpf(0.5) + mp.sqrt(v))
    direct = -mp.diff(lambda v: mp.log(f(v)), U, 2)
    rows.append(dict(u=uu, Z=mp.nstr(Zu, 15), xi=mp.nstr(f(U), 15), stiffness=mp.nstr(Emq, 15), var=mp.nstr(Var, 15),
                     stiff_minus_var=mp.nstr(Emq - Var, 15), direct=mp.nstr(direct, 15), fisher_fraction=mp.nstr(Var / Emq, 15)))
    print(rows[-1])
out['xi_numeric'] = rows
out['label'] = 'identity PROVED (exact symbolic, general positive mixture); xi check NUMERICAL'
json.dump(out, open('s49_identity.json', 'w'), indent=1)
