"""Round 14B shared Galerkin for the Weil-reflected pairing <f,f>_Theta = W(f * Theta f), Theta f(x) = x^-1 conj f(1/x).
In u = log x with g(u) = f(e^u) e^{u/2} (support [-L/2, L/2]) the pairing is the Round-13 Weil form
  W(g) = C A(0) + 2 int_0^inf w(t)[A(0)-A(t)] dt + sum_a  2 int_R A(t) cosh(a t) dt  -  2 sum_n c(n) n^-1/2 A(log n),
  A(t) = int g(x) g(x+t) dx.
Basis: Dirichlet sines phi_j on [-L/2, L/2]; the autocorrelation matrix A_ij(t) = int phi_i(x) phi_j(x+t) dx is used in
CLOSED FORM (no quadrature in x), so high frequencies (needed for the DH zero at height 85.7) are exact.
Integer data only: c(n) (from the model's own Dirichlet coefficients), arch (C, w) from the model's own gamma factor,
pole exponents a from the model's own poles. No zeros enter."""
import numpy as np, math
from scipy.integrate import quad
XG, WG = np.polynomial.legendre.leggauss(64)

def Aclosed(L, M, t):
    """A_ij(t) = int_{-L/2}^{L/2-t} phi_i(x) phi_j(x+t) dx, phi_j = sqrt(2/L) sin(j pi (x+L/2)/L), 0 <= t."""
    if t >= L: return np.zeros((M, M))
    a = np.arange(1, M+1)*np.pi/L; ai = a[:, None]; aj = a[None, :]; Y = L - t
    d = ai - aj; s = ai + aj
    with np.errstate(divide='ignore', invalid='ignore'):
        T1 = np.where(np.abs(d) > 1e-14, (np.sin(d*Y - aj*t) - np.sin(-aj*t))/np.where(np.abs(d) > 1e-14, d, 1.0),
                      Y*np.cos(aj*t))
    T2 = (np.sin(s*Y + aj*t) - np.sin(aj*t))/s
    return (T1 - T2)/L

def basis(L, M, x):
    j = np.arange(1, M+1)[:, None]
    return np.where(np.abs(x) <= L/2 + 1e-15, math.sqrt(2/L)*np.sin(j*np.pi*(x+L/2)/L), 0.0)

def coshvec(L, M, a):
    """c_j = int phi_j cosh(a x), s_j = int phi_j sinh(a x) (closed form)."""
    k = np.arange(1, M+1)*np.pi/L; nrm = math.sqrt(2/L)
    # int_{-L/2}^{L/2} sin(k(x+L/2)) e^{a x} dx = e^{-aL/2} int_0^L sin(k y) e^{a y} dy
    def I(b):  # int_0^L sin(k y) e^{b y} dy
        return (np.exp(b*L)*(b*np.sin(k*L) - k*np.cos(k*L)) + k)/(b*b + k*k)
    ep = math.exp(-a*L/2)*I(a); em = math.exp(a*L/2)*I(-a)
    return nrm*0.5*(ep + em), nrm*0.5*(ep - em)

def operators(model, L, M, nt=64*60):
    I = np.eye(M)
    arch = np.zeros((M, M))
    if model.get('C') is not None:
        w = model['w']
        arch = model['C']*I + 2*I*quad(lambda t: float(w(t)), L, np.inf, limit=400)[0]
        edges = np.concatenate([[0], L*np.geomspace(1e-7, 1, 70)])
        for a, b in zip(edges[:-1], edges[1:]):
            for xg, wg in zip(0.5*(b-a)*XG+0.5*(b+a), 0.5*(b-a)*WG):
                A = Aclosed(L, M, xg); arch += 2*wg*float(w(xg))*(I - 0.5*(A+A.T))
    pole = np.zeros((M, M))
    for a, wt in model.get('poles', []):
        c, s = coshvec(L, M, a); pole += wt*2*(np.outer(c, c) - np.outer(s, s))
    prime = np.zeros((M, M))
    m = np.log(model['n']) < L
    for n, cc in zip(model['n'][m], model['c'][m]):
        A = Aclosed(L, M, math.log(n)); prime += cc/math.sqrt(n)*(A+A.T)
    return arch, pole, prime

def peak_freq(L, v, nf=4000, fmax=None):
    M = len(v); fmax = fmax or (M+4)*np.pi/L
    x = np.linspace(-L/2, L/2, 4001); g = v @ basis(L, M, x)
    fr = np.linspace(0, fmax, nf); G = np.abs(np.trapezoid(g[None, :]*np.exp(1j*fr[:, None]*x[None, :]), x, axis=1))
    return float(fr[np.argmax(G)])
