"""Grok S70 review: exact sympy checks of Codex Oct 17 algebra claims (Moyal [x^3,y^3], FS antipode T_k^2=(-1)^k,
half-density mirror reversal of the Schwinger cocycle, D^2 cocycle n(n^4-1)/30)."""
import sympy as sp
x,y,h=sp.symbols('x y h')
def Pk(f,g,k):  # k-th power of Poisson bidifferential (d_x x d_y - d_y x d_x)
    return sum(sp.binomial(k,i)*(-1)**i*sp.diff(f,x,k-i,y,i)*sp.diff(g,x,i,y,k-i) for i in range(k+1))
def star(f,g,N=8): return sp.expand(sum((h/2)**k/sp.factorial(k)*Pk(f,g,k) for k in range(N)))
c=sp.expand((star(x**3,y**3)-star(y**3,x**3))/h); print("[x^3,y^3]_*/h =",c, " expected 9x^2y^2+3h^2/2:",sp.simplify(c-(9*x**2*y**2+sp.Rational(3,2)*h**2))==0)
# FS antipode: T e_j = (-1)^j e_{k-j}, antilinear with real coefficients
ok=True
for k in range(0,9):
    T=lambda j:((-1)**j,k-j)
    for j in range(k+1):
        s1,j1=T(j); s2,j2=T(j1); ok&= (j2==j and s1*s2==(-1)**k)
print("T_k^2=(-1)^k for k=0..8:",ok)
# direct check of T on polynomials
z0,z1=sp.symbols('z0 z1'); okd=True
for k in range(0,6):
    for j in range(k+1):
        e=z0**j*z1**(k-j)
        Te=sp.conjugate(e.subs({z0:-sp.conjugate(z1),z1:sp.conjugate(z0)},simultaneous=True))
        okd&= sp.expand(Te-(-1)**j*z0**(k-j)*z1**j)==0
print("T e_j = (-1)^j e_{k-j} directly, k<=5:",okd)
# cocycle Psi(A_n,f;A_-n,g)=sum_{j=1}^n f(-j)g(n-j); mirror z^n f(D)->z^-n f(-D-1); antisymmetry Psi(A_-n,f;A_n,g)=-Psi(A_n,g;A_-n,f)
D=sp.symbols('D'); a=sp.symbols('a0:5'); b=sp.symbols('b0:5')
f=sum(a[i]*D**i for i in range(5)); g=sum(b[i]*D**i for i in range(5))
Psi=lambda n,F,G: sum(F.subs(D,-j)*G.subs(D,n-j) for j in range(1,n+1))
okm=True
for n in range(1,9):
    ft=f.subs(D,-D-1); gt=g.subs(D,-D-1)
    lhs=-Psi(n,gt,ft)   # Psi(mirror A_n f, mirror A_-n g) = Psi(A_-n ft; A_n gt) = -Psi(A_n gt; A_-n ft)
    okm&= sp.expand(lhs+Psi(n,f,g))==0
print("half-density mirror reverses Psi (quartic f,g, n<=8):",okm)
print("D^2 cocycle = n(n^4-1)/30 for n<=12:",all(Psi(n,D**2,D**2)==sp.Rational(n*(n**4-1),30) for n in range(1,13)))
