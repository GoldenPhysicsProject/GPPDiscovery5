"""Round 10 item 3 (float NUMERICAL test of a CHECKED identity).
Identity: for x_a in X' and j>=1, C~_{a,j} := omega_j^2 B(x_a, psi_j) = <k_a, phi_j>_{L2[0,l]}, k_a(y) = -d^2/dy^2 B(x_a, delta_y)
(Green: h_a(0)=0 since G(x,0)=0). Hence CC^T = Gram(k) - H H^T exactly, H_{a,j} = C~_{a,j}, j=1..72 (Bessel/Parseval on the complete
sine basis). Test: (i) <k_a,phi_j> vs spectral C for j in 73..412; (ii) Gram(k)-HH^T vs spectral C C^T (K=412) + Parseval tail."""
import numpy as np, json
from scipy.integrate import quad
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
l = 0.25
h = lambda u: 2*np.exp(-u/2)/(1-np.exp(-2*u))
w0 = float(special.digamma(0.25) - np.log(np.pi))
def tG2(y):
    a = np.abs(y); return np.where(a <= 1, -2+3*a, np.where(a <= 2, 2-a, 0.0))
def arch2(d):
    f = lambda u: h(u)*(tG2(d) - 0.5*(tG2(d+u)+tG2(d-u)))
    pts = sorted(set([abs(d), 1-abs(d), 1+abs(d), 2-abs(d), 2+abs(d)]))
    s = 0.0; a = 0.0
    for p in pts + [3.0]:
        if p > a: s += quad(f, a, p, limit=200, epsabs=1e-14, epsrel=1e-13)[0]; a = p
    s += float(tG2(d))*quad(h, 3.0, np.inf, epsabs=1e-14)[0]
    return w0*float(tG2(d)) + s
def prime2(d): return sum(np.log(p)/np.sqrt(n)*(tG2(np.log(n)+d)+tG2(-np.log(n)+d)) for n,p in NS)
dd = np.linspace(0, 0.5, 5001)
A2 = np.array([arch2(d) for d in dd])
def K2(d):
    a = np.abs(d)
    return 0.5*bI2**2*np.cosh(d/2) - prime2(d) + np.interp(a, dd, A2)
n = 4000; xs = np.linspace(0, l, n+1); wx = np.full(n+1, l/n); wx[0] = wx[-1] = l/(2*n)
X, Yg = np.meshgrid(xs, xs)          # rows y, cols x
Gyy = 0.5*(K2(X - Yg) - K2(X + Yg))
phi = lambda j, x: np.sqrt(2/l)*np.sin(j*np.pi/l*x)
N = 12; JY = range(13, 73)
kphi = {i: -(Gyy*(wx*phi(i, xs))[None, :]).sum(1) for i in range(1, 73)}
kd = -0.5*(K2(l - xs) - K2(l + xs))
w_ = lambda j: np.sqrt(2/l)*(-1)**(j+1)/om(j)
rows = [kphi[i] for i in range(1, N+1)] + [kd] + [kphi[k] - w_(k)*kd for k in JY]
Kr = np.array(rows)                  # 73 x (n+1)
Gram = (Kr*wx) @ Kr.T
Phi = np.array([phi(j, xs) for j in range(1, 2001)])
Hall = (Kr*wx) @ Phi.T               # <k_a, phi_j>, j = 1..2000
H72 = Hall[:, :72]
CCt_id = Gram - H72 @ H72.T
# spectral C for j = 73..412
Fx = [S(i) for i in range(1, N+1)] + [Fd] + [S(k) - w_(k)*Fd for k in JY]
sx = [sig(i) for i in range(1, N+1)] + [sd] + [sig(k) - w_(k)*sd for k in JY]
Jz = range(73, 413)
Csp = np.array([[ (B(Fa, sa, S(j) - w_(j)*Fd, sig(j) - w_(j)*sd))*om(j)**2 for j in Jz] for Fa, sa in zip(Fx, sx)])
Cx = Hall[:, 72:412]
rel = np.abs(Cx - Csp).max()/np.abs(Csp).max()
CCt_sp = Csp @ Csp.T
CCt_x412 = Cx @ Cx.T
tailx = Hall[:, 412:] @ Hall[:, 412:].T
res = dict(max_rel_C_x_vs_spectral_73_412=float(rel),
  C_row_norm_delta_x=float(np.sqrt(CCt_id[N, N])), C_row_norm_delta_sp412=float(np.sqrt(CCt_sp[N, N])),
  diag_CCt_identity_first=[float(CCt_id[i, i]) for i in [0, 6, 11, 12, 13, 40, 72]],
  diag_CCt_x_412_plus_x_tail_to2000=[float(CCt_x412[i, i] + tailx[i, i]) for i in [0, 6, 11, 12, 13, 40, 72]],
  diag_CCt_spectral_412=[float(CCt_sp[i, i]) for i in [0, 6, 11, 12, 13, 40, 72]],
  max_abs_diff_identity_vs_x2000=float(np.abs(CCt_id - CCt_x412 - tailx).max()),
  max_abs_CCt=float(np.abs(CCt_id).max()),
  Parseval_check_head_vs_gram_delta=[float(Gram[N, N]), float((Hall[N]**2).sum())])
print(json.dumps(res, indent=1)); json.dump(res, open('kgram_float_results.json', 'w'), indent=1)
np.save('CCt_identity.npy', CCt_id); np.save('Csp.npy', Csp)
