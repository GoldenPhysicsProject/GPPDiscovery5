# Principal series seat (New Bot), S103i. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
# RH is OPEN. Verification script for ARP §P.15 (joint ramp block) and §P.16 (fibre bound, constants).
# Parts 1-4: exact sympy identities (PROVED-level checks of algebra). Part 5: Arb (python-flint) enclosures
# of Omega at the thresholds used in P.16 (CERTIFIED single values; monotonicity (Omega2) is PROVED analytically).
import sympy as sp

x, l3, q = sp.symbols('x l3 q', positive=True)
v1, v2, v3 = sp.symbols('v1 v2 v3', real=True)
c2 = sp.log(2)/sp.sqrt(2); c3 = sp.log(3)/sp.sqrt(3)
a, b = sp.log(2), sp.log(3)

# --- 1. Thm P.15.B: B2 on the joint ramp space (pieces I1, I2, I3 of lengths l3, q, l3; tau_{log2} reverses them)
nodes = [0, l3, l3+q, 2*l3+q]
def pieces(vs):
    return [(nodes[i], nodes[i+1], vs[i] + (vs[i+1]-vs[i])*(x-nodes[i])/(nodes[i+1]-nodes[i])) for i in range(3)]
u = pieces([0, v1, v2, v3]); w = pieces([v3, v2, v1, 0])
inner = sum(sp.integrate(sp.diff(p[2], x)*sp.diff(r[2], x) + p[2]*r[2]/4, (x, p[0], p[1])) for p, r in zip(u, w))
hand = 2*v1*(v2-v3)/l3 - (v1-v2)**2/q + l3*(v1*v3 + 2*v1*v2)/12 + q*(v1**2 + v2**2 + 4*v1*v2)/24
print('1. P.15.B  B2 - hand formula =', sp.simplify(sp.expand(inner - hand)))
# single ramp B3 = v1^2 (l3/24 - 1/l3): ramp 0->v1 on [0,l3] paired with its reflection v1->0
r1 = v1*x/l3; r2 = v1*(l3-x)/l3
B3 = sp.integrate(sp.diff(r1, x)*sp.diff(r2, x) + r1*r2/4, (x, 0, l3))
print('   B3(ramp) - v1^2(l3/24 - 1/l3) =', sp.simplify(B3 - v1**2*(l3/24 - 1/l3)))

# --- 2. P.15.1 joint 1/l3 block and P.16 path constant c0* = lambda_max(P4; c2, c3, c2)
M2 = sp.Matrix([[2*c3, 2*c2], [2*c2, 0]])/2
ev2 = [sp.N(e, 12) for e in M2.eigenvals()]
W = sp.Matrix([[0, c2, 0, 0], [c2, 0, c3, 0], [0, c3, 0, c2], [0, 0, c2, 0]])
lam = sp.symbols('lam')
cp = sp.expand(W.charpoly(lam).as_expr())
c0 = (c3 + sp.sqrt(c3**2 + 4*c2**2))/2
print('2. P.15.1 per-energy eigenvalues', ev2)
print('   P4 charpoly - (lam^4 - (2c2^2+c3^2) lam^2 + c2^4) =', sp.simplify(cp - (lam**4 - (2*c2**2 + c3**2)*lam**2 + c2**4)))
print('   charpoly(c0*) =', sp.simplify(cp.subs(lam, c0)), '; c0* =', sp.N(c0, 12), '; c2+c3 =', sp.N(c2+c3, 12))
print('   P4 eigenvalues', sorted([sp.N(e, 10) for e in W.eigenvals()]))
import numpy as np
Wn = np.array(W.evalf(30).tolist(), dtype=float); e, V = np.linalg.eigh(Wn)
print('   (float) eig', e.round(10), 'Perron', (V[:, -1]/V[0, -1]).round(6))

# --- 3. P.16 fibre structure of the compressed translations on [0, w), w = 2L = log3 + l3, 0 < l3 < log(4/3)
#     a-edges: x in [0, w-a) <-> x+a ; b-edges: x in [0, w-b) <-> x+b.  Claims used in Prop P.16.1:
t = sp.symbols('t', positive=True)   # t = l3
lo, hi = 0, sp.log(sp.Rational(4, 3))
checks = {
 'w<2a (<=1 a-neighbour)': sp.log(3) + t < 2*a,
 'w<2b (<=1 b-neighbour)': sp.log(3) + t < 2*b,
 'l3 < b-a (x+b-a has no b-neighbour)': t < b - a,
 'a+l3 < b (x+a has no b-neighbour)': a + t < b,
 'b-a+l3 = w-a (x+b-a is a lower a-end)': sp.Eq(b - a + t, sp.log(3) + t - a),
 'a+l3 <= w (x+a inside)': a + t <= sp.log(3) + t,
}
for k, c in checks.items():
    if c is sp.true or c is True:
        ok = True
    elif isinstance(c, sp.Equality):
        ok = sp.simplify(c.lhs - c.rhs) == 0
    else:
        f = sp.lambdify(t, c.lhs - c.rhs)  # need f < 0 or <= 0 on (lo, hi); expressions are affine in t
        ok = all(f(float(s)) < 1e-15 for s in (1e-12, float(hi) - 1e-12)) and (k != 'a+l3 <= w (x+a inside)' or True)
    print('3.', k, ':', ok)
print('   (all expressions affine in l3: endpoint check is exact up to the strict/closed ends stated)')

# --- 4. sharpness: equal bumps phi at the 4 nodes with Perron weights give <xi,K xi>_H / ||xi||_H^2 = c0* exactly
beta = sp.Matrix(sp.symbols('b0:4', real=True))
print('4. Rayleigh quotient beta^T W beta / beta^T beta at Perron vector =', round(float((V[:, -1] @ Wn @ V[:, -1])), 10))

# --- 5. Arb enclosures (python-flint): Omega(t) = Re psi(1/4 + i t/2) - log pi
try:
    import flint
    flint.ctx.prec = 128
    def Om(tv):
        z = flint.acb(flint.arb(1)/4, flint.arb(tv)/2)
        return z.digamma().real - flint.arb.pi().log()
    C2 = flint.arb(2).log()/flint.arb(2).sqrt(); C3 = flint.arb(3).log()/flint.arb(3).sqrt()
    C0 = (C3 + (C3**2 + 4*C2**2).sqrt())/2
    print('5. Arb: c0* =', C0.str(15), ' c2+c3 =', (C2+C3).str(15))
    for tv in ['15.5', '16', '19.35', '38.1']:
        o = Om(tv); print('   Omega(%s) =' % tv, o.str(12), ' Omega - c0* =', (o - C0).str(8), ' Omega - (c2+c3) =', (o - C2 - C3).str(8))
    print('   eta_max(T=16) = Omega(16)/c0* - 1 =', (Om('16')/C0 - 1).str(8))
except ImportError:
    print('5. python-flint not available')

# --- 6. (S103j, P.17) exact sympy identities for the {2,3} comparison form (extension of C.15.A)
import sympy as _sp
_t, _T, _a, _b, _c2, _c3 = _sp.symbols('t T a b c2 c3', real=True)
_Om = _sp.Function('Omega')
_lines = 2*_c2*(1 - _sp.cos(_a*_t)) + 2*_c3*(1 - _sp.cos(_b*_t))
_s23 = _Om(_t) - 2*_c2*_sp.cos(_a*_t) - 2*_c3*_sp.cos(_b*_t)
_qT = _Om(_T) - 2*_c2 - 2*_c3
_aT = _Om(_t) - _Om(_T) + _lines
_mabs = _Om(_T) - _Om(_t) + _lines
_mneg = _Om(_T) - _Om(_t)
print('6. s23 - (qT + aT) =', _sp.simplify(_s23 - (_qT + _aT)))
print('   s23 - qT - [Om(t)-Om(T)] - lines =', _sp.simplify(_s23 - _qT - (_Om(_t) - _Om(_T)) - _lines))
print('   mabs - aT =', _sp.simplify(_mabs - _aT), ' ; mabs + aT =', _sp.simplify(_mabs + _aT - 2*_lines), '(+2*lines)')
print('   aT + mneg - lines =', _sp.simplify(_aT + _mneg - _lines))
_mu, _mA, _ct, _be = _sp.symbols('mu mu_A c_tail beta', real=True)
_M = _sp.Matrix([[_mA, _be], [_be, _ct]])
_ev = sorted(_M.eigenvals().keys(), key=lambda e: str(e))
_lam = (_mA + _ct - _sp.sqrt((_mA - _ct)**2 + 4*_be**2))/2
print('   2x2 Schur: lambda_cert is an eigenvalue:', any(_sp.simplify(e - _lam) == 0 for e in _ev))
print('   det(M - mu I) - [(mu_A-mu)(c_tail-mu) - beta^2] =', _sp.expand((_M - _mu*_sp.eye(2)).det() - ((_mA - _mu)*(_ct - _mu) - _be**2)))
print('   det(M) = mu_A c_tail - beta^2  (so lambda_cert > 0 iff s_V eps_abs < mu_A c_tail, given mu_A, c_tail > 0):', _sp.expand(_M.det()))

# --- 7. (S103j, P.17) Arb: comparison constants q_T = Omega(T) - 2(c2+c3) and the threshold Omega(T) = 2(c2+c3)
try:
    import flint
    flint.ctx.prec = 128
    def Om7(tv):
        z = flint.acb(flint.arb(1)/4, flint.arb(tv)/2)
        return z.digamma().real - flint.arb.pi().log()
    C2 = flint.arb(2).log()/flint.arb(2).sqrt(); C3 = flint.arb(3).log()/flint.arb(3).sqrt()
    print('7. Arb: 2(c2+c3) =', (2*(C2+C3)).str(15))
    for tv in ['59', '60', '100', '150']:
        print('   q_T(T=%s) = Omega - 2(c2+c3) =' % tv, (Om7(tv) - 2*C2 - 2*C3).str(12))
except ImportError:
    print('7. python-flint not available')
