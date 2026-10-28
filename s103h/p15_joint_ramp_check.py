# Grok S103h audit: independent check of Principal Thm P.15.B (B2 on joint ramps)
# and the p=3-entry per-energy eigenvalues. Exact (sympy piecewise integration).
import sympy as sp
x,l3,q=sp.symbols('x l3 q',positive=True)
v1,v2,v3=sp.symbols('v1 v2 v3',real=True)
nodes=[0,l3,l3+q,2*l3+q]; vals=[0,v1,v2,v3]
def pieces(vs):
    return [(nodes[i],nodes[i+1],vs[i]+(vs[i+1]-vs[i])*(x-nodes[i])/(nodes[i+1]-nodes[i])) for i in range(3)]
u=pieces(vals); w=pieces([v3,v2,v1,0])  # J2 u: reversal maps pieces (l3,q,l3) to themselves
inner=sum(sp.integrate(sp.diff(a[2],x)*sp.diff(b[2],x)+a[2]*b[2]/4,(x,a[0],a[1])) for a,b in zip(u,w))
hand=2*v1*(v2-v3)/l3-(v1-v2)**2/q+l3*(v1*v3+2*v1*v2)/12+q*(v1**2+v2**2+4*v1*v2)/24
print('B2 - hand =',sp.simplify(sp.expand(inner-hand)))
c2=sp.log(2)/sp.sqrt(2); c3=sp.log(3)/sp.sqrt(3)
print('eig',[sp.N(e,8) for e in sp.Matrix([[c3,c2],[c2,0]]).eigenvals()])
