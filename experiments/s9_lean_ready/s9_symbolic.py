"""S9/S10 Lean-ready symbolic checks (sympy, exact). Grok, 2026-10-14."""
import sympy as sp
t,th,h,u,x,y,ph=sp.symbols('t theta h u x y phi',real=True)
k=sp.symbols('k',integer=True,nonzero=True)
c=sp.symbols('c',real=True)
# (a) algebraic step: cos(h(t +- i theta)) at t=2pi k/h equals cosh(h theta)
tk=2*sp.pi*k/h
for s in (1,-1):
    e=sp.expand_trig(sp.cos(h*(tk+s*sp.I*th)).rewrite(sp.exp)).rewrite(sp.cosh)
    assert sp.simplify(sp.expand(sp.cos(h*tk+s*sp.I*h*th),trig=True)-sp.cosh(h*th))==0
print("(a) cos(h(t_k +- i th)) = cosh(h th): OK")
# (b) Re (t+i th)^{-2} = (t^2-th^2)/(t^2+th^2)^2
re=sp.re(sp.expand_complex(1/(t+sp.I*th)**2))
assert sp.simplify(re-(t**2-th**2)/(t**2+th**2)**2)==0
print("(b) Re (t+i th)^-2 identity: OK")
# (c) full closed form: FT of 2cosh(th u)(1-|u|/h)_+ at xi, integrate exactly
xi=sp.symbols('xi',real=True)
hp=sp.symbols('hp',positive=True)
F=2*sp.integrate(2*sp.cosh(th*u)*(1-u/hp)*sp.cos(xi*u),(u,0,hp))
F=F.subs(xi,2*sp.pi*k/hp)
target=4*(1-sp.cosh(hp*th))*((2*sp.pi*k/hp)**2-th**2)/(hp*((2*sp.pi*k/hp)**2+th**2)**2)
val={th:sp.Rational(3,7),hp:sp.Rational(5,2),k:3}
print("(c) FT closed form residual (exact subs):",sp.N(F.subs(val)-target.subs(val),30))
# (d) sign lemma: c>1, t^2>th^2, h>0 => negative; factorization exact
expr=4*(1-c)*(t**2-th**2)/(h*(t**2+th**2)**2)
print("(d) sign lemma numerator factor 4(1-c)(t^2-th^2): each factor sign fixed; sympy:",
      sp.ask(sp.Q.negative(expr), sp.Q.positive(c-1)&sp.Q.positive(t**2-th**2)&sp.Q.positive(h)))
# (e) Codex 2x2 Gram: vector (1,-1), energy 2c... with w(0)=2, w(u)=2cosh: 4(1-cosh)
G=sp.Matrix([[2,2*c],[2*c,2]]); v=sp.Matrix([1,-1])
print("(e) Gram energy:",sp.factor((v.T*G*v)[0]))
# (f) normalization: -L'/L coefficient at p^k for alpha=p^z, beta=p^-z is (alpha^k+beta^k)log p = 2cosh(z k log p) log p
p=sp.symbols('p',positive=True); z=th+sp.I*ph
a=p**(k*z)+p**(-k*z); U=k*sp.log(p)
assert sp.simplify(a.rewrite(sp.exp)-2*sp.cosh(z*U).rewrite(sp.exp))==0
herm=sp.simplify(sp.re(sp.expand_complex(2*sp.cosh(z*U))))
print("(f) Re a_k =",herm," (= 2cosh(th u)cos(phi u) at u=k log p)")
# (g) untempered with |a(p)|<=2: z=1/10+i*pi/(2 log 2), p=2
zz=sp.Rational(1,10)+sp.I*sp.pi/(2*sp.log(2))
print("(g) |a(2)| for theta=0.1 untempered:",sp.N(sp.Abs(2*sp.cosh(zz*sp.log(2))),15))
# (d') sign lemma via positive reparametrization c=1+A, t^2=th^2+B, A,B,h,th^2>=0
A,B,H=sp.symbols('A B H',positive=True); T2=sp.symbols('T2',nonnegative=True)
e2=4*(1-(1+A))*B/(H*(2*T2+B)**2)
print("(d') sign lemma is negative:", sp.ask(sp.Q.negative(e2)), "| expr =",sp.factor(e2))
