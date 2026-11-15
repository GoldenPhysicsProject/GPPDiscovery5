"""Independent cross-check (NUMERICAL, float64, different code path) of P1 block entries: prime and pole parts directly in the
kernel domain, 1/2∫∫[k(s-u)-k(s+u)] phi_j(s) phi_k(u), with the in-tree triG prime kernel (composite Gauss-Legendre, 2000^2 nodes);
arch part by scipy quad of the spectral integral with the explicit closed-form S_j. Compare with p1_oddblock_L05_N12_results.json."""
import numpy as np, json
from scipy import integrate, special
l = 0.25
NS = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
def triG(y):
    a = np.abs(y)
    return np.where(a <= 1, 2/3 - y*y + a**3/2, np.where(a <= 2, (2-a)**3/6, 0.0))
prime = lambda d: sum(np.log(p)/np.sqrt(n)*(triG(np.log(n)+d)+triG(-np.log(n)+d)) for n,p in NS)
bI2 = 16*np.sinh(l)**2
pole = lambda d: 2*bI2**2*np.cosh(d/2)
x0, w0 = np.polynomial.legendre.leggauss(20)
edges = np.linspace(0, l, 101)
X = np.concatenate([(a+b)/2 + (b-a)/2*x0 for a, b in zip(edges[:-1], edges[1:])])
Wt = np.concatenate([(b-a)/2*w0 for a, b in zip(edges[:-1], edges[1:])])
s, u = np.meshgrid(X, X, indexing='ij')
phi = lambda j, x: np.sqrt(2/l)*np.sin(j*np.pi/l*x)
Gp = (prime(s-u) - prime(s+u))/2; Gq = (pole(s-u) - pole(s+u))/2
def part(G, j, k): return phi(j, X) @ (G*np.outer(Wt, Wt)) @ phi(k, X)
def Sj(j, t):
    w = j*np.pi/l
    return w*l*np.sinc((t-w)*l/np.pi)/(t+w)
def arch(j, k):
    g = lambda t: (np.sinc(t/(2*np.pi)))**4*(special.digamma(0.25+0.5j*t).real - np.log(np.pi))*Sj(j,t)*Sj(k,t)
    tot = sum(integrate.quad(g, a, a+10, epsabs=1e-17, epsrel=1e-13, limit=200)[0] for a in range(0, 3000, 10))
    return (2/l)*tot/np.pi
R = json.load(open('p1_oddblock_L05_N12_results.json'))
out = []
for (j, k) in [(1,1), (1,2), (3,4), (6,7), (12,12)]:
    v = arch(j,k) - part(Gp,j,k) + part(Gq,j,k)
    ref = float(R['M'][j-1][k-1])
    out.append(dict(j=j, k=k, crosscheck=f"{v:.10e}", arb=f"{ref:.10e}", rel_diff=f"{abs(v-ref)/abs(ref):.2e}"))
    print(out[-1], flush=True)
json.dump(out, open('crosscheck_results.json', 'w'), indent=1)
