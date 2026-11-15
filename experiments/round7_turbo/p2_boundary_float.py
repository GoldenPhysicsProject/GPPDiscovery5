"""Round 6 (Turbo 6889647810): structure of the odd-face form on E_12^perp at L=0.5. Float64 spectral quadrature (NUMERICAL test of a
CHECKED identity). Same exact reduction as Round 5: B(f,g) = (1/pi)∫Psi^2 V F G dt - 2 Psi(i/2)^2 sigma_f sigma_g,
F(t) = ∫_0^l f(x) sin(tx) dx, sigma_f = ∫_0^l f sinh(x/2). Boundary functional delta_l: F = sin(t l), sigma = sinh(l/2)."""
import numpy as np, json
from scipy import special
l = 0.25; NS = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
bI2 = 16*np.sinh(l)**2
T = 8000.0; nper = 4000
x0, w0 = np.polynomial.legendre.leggauss(40)
ed = np.linspace(0, T, nper+1)
t = np.concatenate([(a+b)/2+(b-a)/2*x0 for a,b in zip(ed[:-1], ed[1:])]); wt = np.concatenate([(b-a)/2*w0 for a,b in zip(ed[:-1], ed[1:])])
V = special.digamma(0.25+0.5j*t).real - np.log(np.pi) - 2*sum(np.log(p)/np.sqrt(n)*np.cos(t*np.log(n)) for n,p in NS)
Psi2 = np.sinc(t/(2*np.pi))**4
meas = wt*Psi2*V/np.pi
om = lambda j: j*np.pi/l
def S(j):  # spectral transform of phi_j (L^2-normalised)
    w = om(j); return np.sqrt(2/l)*w*l*np.sinc((t-w)*l/np.pi)/(t+w)
def sig(j):
    w = om(j); return -np.sqrt(2/l)*(-1)**j*w*np.sinh(l/2)/(0.25+w*w)
def B(F1, s1, F2, s2): return meas @ (F1*F2) - 2*bI2**2*s1*s2
Fd, sd = np.sin(t*l), np.sinh(l/2)
g = B(Fd, sd, Fd, sd)
# kernel-domain value of g = G(l,l) = 1/2[K(0) - K(2l)]: arch by the same t-grid, prime by triG, pole closed form
def triG(y):
    a = abs(y); return 2/3-y*y+a**3/2 if a <= 1 else ((2-a)**3/6 if a <= 2 else 0.0)
prime = lambda d: sum(np.log(p)/np.sqrt(n)*(triG(np.log(n)+d)+triG(-np.log(n)+d)) for n,p in NS)
archw = wt*Psi2*(special.digamma(0.25+0.5j*t).real-np.log(np.pi))/np.pi
arch = lambda d: archw @ np.cos(t*d)
K = lambda d: 2*bI2**2*np.cosh(d/2) - prime(d) + arch(d)
g_kernel = 0.5*(K(0.0) - K(2*l))
w_ = lambda j: np.sqrt(2/l)*(-1)**(j+1)/om(j)
u = [B(S(i), sig(i), Fd, sd) for i in range(1, 13)]
out = dict(g_spectral=g, g_kernel=g_kernel, u=u)
# leading-order test M_jk ~ w_j w_k g (j,k >= 13), M_ij ~ u_i w_j
rows = []
for (j,k) in [(13,13),(13,14),(13,20),(20,20),(20,31),(30,30),(40,40),(40,41),(60,60),(60,90),(90,90)]:
    m = B(S(j), sig(j), S(k), sig(k)); lead = w_(j)*w_(k)*g
    rows.append(dict(j=j,k=k,M=m,lead=lead,ratio=m/lead))
out['high'] = rows
cr = []
for i in [1,6,12]:
    for j in [13,30,60,90]:
        m = B(S(i), sig(i), S(j), sig(j)); lead = u[i-1]*w_(j)
        cr.append(dict(i=i,j=j,M=m,lead=lead,ratio=m/lead))
out['cross'] = cr
# bordered Gram on {phi_1..phi_12, delta_l}
N = 12
M = np.array([[B(S(i),sig(i),S(k),sig(k)) for k in range(1,N+1)] for i in range(1,N+1)])
R5 = json.load(open('../round5_turbo/p1_oddblock_L05_N12_results.json'))
Mref = np.array([[float(x) for x in r] for r in R5['M']])
out['max_rel_vs_arb_block'] = float(np.max(np.abs(M-Mref))/np.max(np.abs(Mref)))
Bd = np.block([[M, np.array(u)[:,None]],[np.array(u)[None,:], np.array([[g]])]])
out['eig_M12'] = sorted(np.linalg.eigvalsh(M).tolist())[:3]
out['eig_bordered13'] = sorted(np.linalg.eigvalsh(Bd).tolist())[:3]
Sc = M - np.outer(u,u)/g
out['eig_schur_M-uu/g'] = sorted(np.linalg.eigvalsh(Sc).tolist())[:3]
print(json.dumps(out, indent=1, default=float))
json.dump(out, open('p2_boundary_float_results.json','w'), indent=1, default=float)
