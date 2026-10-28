# Grok S102: symbolic replays (no numerics). RH OPEN.
# (1) G_theta' twist identity (Goddess of PNT 3c7f048 / Principal series): kernel and zero-side parallelogram.
# (2) F_theta polarization: 2cosh(theta r) g(r) = h+ * h-~ + h- * h+~, h± = e^{±theta u} f.
# (3) Astra Oct 27 pole/first-prime Ward identity, split-integral form (checked by hand in ARP.md §G.3; here on an explicit even G0).
# (4) Astra Oct 26 mirror overlap identity 1 - R_p^2 = (a-b)^2/(1-ab)^2, a=p^-s, b=p^-(1-s).
# (5) Local gamma factors on Re s = 1/2: |gamma_p| = 1 for zeta; for F_theta each GL(1) constituent has |.| != 1, the product has |.| = 1.
import sympy as sp
u, r, t, th, v, w = sp.symbols('u r t theta v w', real=True)
print('(1) kernel:', sp.simplify(sp.expand_trig(sp.cos(t*u)*sp.cos(t*(u-r)) + sp.sin(t*u)*sp.sin(t*(u-r)) - sp.cos(t*r))))
A, B, C, D = sp.symbols('A B C D', real=True); a1 = A + sp.I*B; b1 = C + sp.I*D
print('(1) zero side:', sp.simplify(sp.expand(sp.Abs(a1 + sp.I*b1)**2 + sp.Abs(a1 - sp.I*b1)**2 - 2*sp.Abs(a1)**2 - 2*sp.Abs(b1)**2)))
print('(2) polarization integrand:', sp.simplify(sp.exp(th*u)*sp.exp(-th*(u - r)) + sp.exp(-th*u)*sp.exp(th*(u - r)) - 2*sp.cosh(th*r)))
# (3) with G0(w) = exp(-w^2) (even, all exponential moments); v > 0 symbolic check of each split piece
G0 = sp.exp(-w**2)
M = sp.integrate(G0*sp.exp(w/2), (w, -sp.oo, sp.oo))
vv = sp.Rational(3, 7)
lhs = M*(sp.exp(vv/2) + sp.exp(-vv/2)) - sp.integrate(sp.exp(u/2)*(G0.subs(w, u + vv) + G0.subs(w, vv - u)), (u, 0, sp.oo))
rhs = sp.integrate(G0*sp.exp(-(vv - w)/2), (w, -sp.oo, vv)) + sp.integrate(G0*sp.exp(-(w - vv)/2), (w, vv, sp.oo))
print('(3) Ward split residual (v=3/7):', sp.N(lhs - rhs, 30))
a, b = sp.symbols('a b', positive=True)
R2 = (1 - a**2)*(1 - b**2)/(1 - a*b)**2
print('(4) mirror overlap:', sp.simplify(1 - R2 - (a - b)**2/(1 - a*b)**2))
p = sp.Integer(7); tt = sp.Rational(13, 10); thv = sp.Rational(1, 5); s = sp.Rational(1, 2) + sp.I*tt
Lp = lambda z: 1/(1 - p**(-z))
gz = Lp(1 - s)/Lp(s)
g_plus = Lp(1 - s - thv)/Lp(s + thv)    # constituent |.|^{theta}: s -> s+theta
g_minus = Lp(1 - s + thv)/Lp(s - thv)
print('(5) |gamma_p zeta| =', sp.N(sp.Abs(gz), 20), ' |gamma_+| =', sp.N(sp.Abs(g_plus), 20), ' |gamma_-| =', sp.N(sp.Abs(g_minus), 20), ' |product| =', sp.N(sp.Abs(g_plus*g_minus), 20))
