"""Check q_delta against -d^2_s G(l,l) = -1/2[K''(0) - K''(2l)] computed in the kernel domain (NUMERICAL)."""
import numpy as np, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
def triG2(y):
    a = abs(y); return -2+3*a if a <= 1 else ((2-a) if a <= 2 else 0.0)
prime2 = lambda d: sum(np.log(p)/np.sqrt(n)*(triG2(np.log(n)+d)+triG2(-np.log(n)+d)) for n,p in NS)
archw = wt*Psi2*(special.digamma(0.25+0.5j*t).real-np.log(np.pi))/np.pi
arch2 = lambda d: -archw @ (t**2*np.cos(t*d))
K2 = lambda d: 2*bI2**2*np.cosh(d/2)/4 - prime2(d) + arch2(d)
print("K''(0) =", K2(0.0), " K''(0.5) =", K2(0.5), " -d2G(l,l) =", -0.5*(K2(0.0)-K2(0.5)))
