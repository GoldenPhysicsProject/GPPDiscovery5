"""Round 17B: symbolic checks for the massless spin-2 conformal primary (PS 2017 sec. 6.3, Mellin representative):
   phi^{Delta,+}_{mu nu}(X; z) = eps+_mu eps+_nu (-q.X)^{-Delta},  q = (1+z zb, z+zb, -i(z-zb), 1-z zb), eps+ = d_z q.
(a) traceless, transverse (Lorenz), and box phi = 0 for every Delta;
(b) under a rotation z -> e^{i th} z the polarization picks up e^{-2 i th} (|J| = 2) and (-q.X)^{-Delta} is invariant,
    under a boost z -> lam z (lam > 0) the profile scales as lam^{-Delta}: the weights are (Delta, J = +-2) independently;
(c) the energy integral that sets normalizability is the scalar one: int omega^{Delta-1} eps eps e^{i omega q.X} = eps eps * (scalar Mellin).
So the spin-2 data enter only through eps eps; the Delta-dependence (and hence the PS condition Re Delta = 1) is spin-blind."""
import sympy as sp, json, os
X = sp.symbols('X0:4'); z, zb = sp.symbols('z zb'); D = sp.symbols('Delta'); th, lam, w = sp.symbols('theta lambda omega', positive=True)
eta = sp.diag(-1, 1, 1, 1)
q = sp.Matrix([1+z*zb, z+zb, -sp.I*(z-zb), 1-z*zb])
ep = q.diff(z)
dot = lambda a, b: (a.T*eta*b)[0]
qX = dot(q, sp.Matrix(X))
epl = eta*ep   # lower-index polarization
phi = [[epl[m]*epl[n]*(-qX)**(-D) for n in range(4)] for m in range(4)]
up = lambda m: eta[m, m]
res = {}
res['q.q'] = sp.simplify(dot(q, q)); res['q.eps'] = sp.simplify(dot(q, ep)); res['eps.eps'] = sp.simplify(dot(ep, ep))
res['trace'] = sp.simplify(sum(up(m)*phi[m][m] for m in range(4)))
res['div'] = [sp.simplify(sum(up(m)*sp.diff(phi[m][n], X[m]) for m in range(4))) for n in range(4)]
res['box'] = [sp.simplify(sum(up(m)*sp.diff(phi[0][n], X[m], 2) for m in range(4))) for n in range(4)]
# rotation about the 3-axis: z -> e^{i th} z, zb -> e^{-i th} zb
qr = q.subs({z: sp.exp(sp.I*th)*z, zb: sp.exp(-sp.I*th)*zb}, simultaneous=True)
R = sp.Matrix([[1, 0, 0, 0], [0, sp.cos(th), -sp.sin(th), 0], [0, sp.sin(th), sp.cos(th), 0], [0, 0, 0, 1]])
res['q_rot = R q'] = sp.simplify((qr - R*q).applyfunc(lambda e: sp.expand(e.rewrite(sp.cos)))) == sp.zeros(4, 1)
epr = ep.subs({z: sp.exp(sp.I*th)*z, zb: sp.exp(-sp.I*th)*zb}, simultaneous=True)
res['eps_rot = e^{-i th} R eps'] = sp.simplify((epr - sp.exp(-sp.I*th)*R*ep).applyfunc(lambda e: sp.expand(e.rewrite(sp.cos)))) == sp.zeros(4, 1)
# boost along 3-axis with rapidity log lam: z -> lam z
ql = q.subs({z: lam*z, zb: lam*zb}, simultaneous=True)
ch, sh = (lam+1/lam)/2, (lam-1/lam)/2
B = sp.Matrix([[ch, 0, 0, -sh], [0, 1, 0, 0], [0, 0, 1, 0], [-sh, 0, 0, ch]])
res['q_boost = lam B q'] = sp.simplify(ql - lam*B*q) == sp.zeros(4, 1)
res['eps_boost = B eps'] = sp.simplify(ep.subs({z: lam*z, zb: lam*zb}, simultaneous=True) - B*ep) == sp.zeros(4, 1)
# (c) energy Mellin of exp(-omega a) (Euclidean continuation of e^{i omega q.X}), a > 0
a = sp.symbols('a', positive=True)
res['mellin'] = sp.simplify(sp.integrate(w**(D-1)*sp.exp(-w*a), (w, 0, sp.oo), conds='none'))
out = {k: str(v) for k, v in res.items()}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r17b_spin2_sympy.json'), 'w'), indent=1)
