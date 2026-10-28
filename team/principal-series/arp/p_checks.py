# Principal series §P: symbolic/exact checks of single derivation steps (no exploratory numerics). RH OPEN.
# Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
import sympy as sp
t, th, x = sp.symbols('t theta x', real=True)
p = sp.symbols('p', positive=True)
# (P1) local factor: L_p(conj z) = conj L_p(z) for real p>1, so on Re s = 1/2:
# |gamma_p(s-th) gamma_p(s+th)| = 1 and |gamma_p(s+th)| = |L_p(1/2+it-th)/L_p(1/2+it+th)|.
Lp = lambda z: 1/(1 - p**(-z))
s = sp.Rational(1, 2) + sp.I*t
gm = Lp(1 - (s - th))/Lp(s - th); gp = Lp(1 - (s + th))/Lp(s + th)
print('(P1) |gamma(s-th)gamma(s+th)|^2 - 1 at p=7, t=13/10, th=1/5 (exact):', sp.nsimplify(sp.N((sp.Abs(gm*gp)**2 - 1).subs({p: 7, t: sp.Rational(13, 10), th: sp.Rational(1, 5)}), 50), tolerance=1e-40))
# symbolic: |1 - p^{-(a+it)}|^2 depends on a only through p^{-a}; conj pairing proof
a = sp.symbols('a', real=True)
mod2 = sp.expand(sp.expand_complex((1 - p**(-a)*sp.exp(-sp.I*t*sp.log(p)))*(1 - p**(-a)*sp.exp(sp.I*t*sp.log(p)))))
print('(P1) |1-p^{-(a+it)}|^2 =', sp.simplify(mod2), '(even in t, so |L_p(a-it)| = |L_p(a+it)|)')
# archimedean: |Gamma(conj z)| = |Gamma(z)| (Schwarz reflection) gives the same product identity for gamma_inf.
G = lambda z: sp.pi**(-z/2)*sp.gamma(z/2)
gim = G(1 - (s - th))/G(s - th); gip = G(1 - (s + th))/G(s + th)
print('(P1inf) |gamma_inf(s-th)gamma_inf(s+th)| at t=13/10, th=1/5:', sp.N(sp.Abs((gim*gip).subs({t: sp.Rational(13, 10), th: sp.Rational(1, 5)})), 30),
      ' |gamma_inf(s+th)|:', sp.N(sp.Abs(gip.subs({t: sp.Rational(13, 10), th: sp.Rational(1, 5)})), 30))
# (P2) Knapp-Stein form on K-type 2m of the GL(2,R) spherical principal series pi(|.|^th,|.|^-th):
# eigenvalue ratio c_{2m}/c_0 = prod_{j=1}^m (j-1/2-th)/(j-1/2+th) (standard; sign convention of weight vectors fixed so c_0>0)
m = sp.symbols('m', integer=True, positive=True)
ratio = lambda mm, tv: sp.prod([(j - sp.Rational(1, 2) - tv)/(j - sp.Rational(1, 2) + tv) for j in range(1, mm + 1)])
print('(P2) all ratios > 0 for th in (0,1/2)? th=0.49, m=1..200:', all(ratio(mm, sp.Rational(49, 100)) > 0 for mm in range(1, 201)))
print('(P2) th=3/4: ratio m=1 =', ratio(1, sp.Rational(3, 4)), ' (negative: form indefinite for th in (1/2,3/2))')
# p-adic spherical constant c_p(th) = zeta_p(2th)/zeta_p(1+2th) > 0 for real th > 0
cp = (1 - p**(-1 - 2*th))/(1 - p**(-2*th))
print('(P2) c_p(th) at p=7, th=1/5:', sp.N(cp.subs({p: 7, th: sp.Rational(1, 5)}), 20))
# (P3) polarization identity on the window: e^{th r} = e^{th u} e^{-th(u-r)}; Q_F = 2Q(f cosh) - 2Q(f sinh)
u, r = sp.symbols('u r', real=True)
print('(P3) multiplicativity residual:', sp.simplify(sp.exp(th*r) - sp.exp(th*u)*sp.exp(-th*(u - r))))
hp, hm = sp.symbols('h_p h_m')  # abstract real bilinear symmetric B: B(h+,h-)+B(h-,h+) = 1/2[B(h++h-)-B(h+-h-)]
B = lambda X, Y: sp.Symbol('B')*X*Y  # one-dimensional model suffices for the algebra
print('(P3) polarization algebra:', sp.simplify(B(hp, hm) + B(hm, hp) - sp.Rational(1, 2)*(B(hp + hm, hp + hm) - B(hp - hm, hp - hm))))
# (P4) 2-atom witness: PD kernel P(r) = 1 + cos r (FT of delta_0 + (delta_1+delta_-1)/2 >= 0), atoms x=(0,1/2), c=(1,-1), th=10
P = lambda rr: 1 + sp.cos(rr)
xs = [0, sp.Rational(1, 2)]
for name, tw in [('cosh (F_theta)', lambda rr: 2*sp.cosh(10*rr)), ('cos (G_theta)', lambda rr: 2*sp.cos(10*rr)), ('untwisted', lambda rr: 1)]:
    M = sp.Matrix(2, 2, lambda i, j: tw(xs[i] - xs[j])*P(xs[i] - xs[j]))
    ev = [sp.N(e, 12) for e in M.eigenvals()]
    print('(P4)', name, 'eigenvalues', ev, ' c^T M c =', sp.N((sp.Matrix([[1, -1]])*M*sp.Matrix([1, -1]))[0], 12))
