# Check: the Galerkin Rayleigh quotient of the probe vector f = e^{-i g0 x} sinh(b x) (projected onto the basis) -> the exact probe value.
import sys, math, numpy as np; sys.path.insert(0,'.')
from galerkin import matrix, phis
from gen_coeffs import atoms
from numpy.polynomial import legendre as Lg
Q=(2,2,3); x0,y0=2.0794141048108840,7.6288675877798635; L=2.0
us,ws,Om0,m=atoms(Q,L)
for N in [24,32,40,48]:
    A=matrix(L,N,us,ws,Om0); gx,gw=Lg.leggauss(200); xx=L*gx; wx=L*gw; F=phis(xx,L,N)
    fr=np.cos(y0*xx)*np.sinh((x0-.5)*xx); fi=-np.sin(y0*xx)*np.sinh((x0-.5)*xx)
    cr=(F*wx)@fr; ci=(F*wx)@fi
    num=cr@A@cr+ci@A@ci; den=cr@cr+ci@ci
    print(N, num/den, '(probe exact: -168.3)', 'lam_min', np.linalg.eigvalsh(A)[0])
