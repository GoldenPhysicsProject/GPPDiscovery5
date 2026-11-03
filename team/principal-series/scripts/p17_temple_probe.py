# NUMERICAL (float64; quadrature to Tmax, no remainders): Temple-Lehmann budget for Q_L, odd/even, 2L = log 3 (+l3).
# v = Ritz vector on the first n modes; residual measured in H* via Parseval over the H-orthonormal Dirichlet basis
# (modes n..M-1 only; the k >= M tail is NOT included -- that is the remainder Celestial must bound).
from p17_probe import *
def qmat(par,L,M,Tmax=3000.,ng=300000):
    t=np.linspace(0,Tmax,ng+1)[1:]; dt=t[1]-t[0]; F=modes(par,M,L,t); w=(0.25+t*t)/np.pi*dt
    s=Om(t)-2*c2*np.cos(a*t)-(2*c3*np.cos(b*t) if 2*L>b+1e-12 else 0)
    Q=np.real((F*(w*s))@F.conj().T); G=np.real((F*w)@F.conj().T)
    # Loewdin-orthonormalise the quadrature Gram (removes the ~4e-4 truncation defect)
    e,U=eigh(G); S=U@np.diag(e**-0.5)@U.T
    return S@Q@S
print('# p17_temple_probe.py -- NUMERICAL. RH OPEN.')
for par,l3 in [('odd',0.0),('odd',0.03),('even',0.0)]:
    M=160; Q=qmat(par,(b+l3)/2,M)
    lam=eigh(Q,eigvals_only=True)
    for n in [10,20,40]:
        ev,V=eigh(Q[:n,:n]); rho=ev[0]; nu2=ev[1]; c=np.zeros(M); c[:n]=V[:,0]
        r=Q@c-rho*c; res=np.linalg.norm(r)
        tb=rho-res**2/(lam[1]-rho)
        print('%s l3=%.2f n=%d rho=%.4e lam1(M=%d)=%.4e lam2=%.4f ||r||_H*(modes<M)=%.3e Temple(nu=lam2)=%.4e Temple(nu=0.44)=%.4e'%(par,l3,n,rho,M,lam[0],lam[1],res,tb,rho-res**2/(0.44-rho)),flush=True)
