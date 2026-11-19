"""Round 17 core: BPY-type laws Phi_M(u) for each model, built from integer data only (theta series /
Dirichlet coefficients + own gamma factor). Zeros are never inputs.

Convention (all models): Xi_M(z) = xi_M(1/2 + i z) = int Phi_M(u) e^{i z u} du,  i.e.
    xi_M(s) = int Phi_M(u) e^{(s-1/2) u} du,      Phi_M even  <=>  xi_M(s) = xi_M(1-s).
BPY form: if Phi_M >= 0 and xi_M(0)=1/2 then U = log Y has density 2 Phi_M(u) e^{-u/2} and E[Y^s] = 2 xi_M(s).

ZETA : xi(s) = s(s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s);  t = e^{2u} in theta(t) = sum e^{-pi n^2 t}
       Phi_Z(u) = sum_{n>=1} (4 pi^2 n^4 e^{9u/2} - 6 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u})   (= 2 x Polya's Phi)
EPS  : xi_E(s) = s(s-1) (sqrt20/2pi)^s Gamma(s) Z_Q(s), Q = m^2+5n^2;  t = e^u in Theta(t) = sum e^{-pi t Q/sqrt5}
       Phi_E(u) = sum'_{Q} r_Q (a^2 e^{5u/2} - 2 a e^{3u/2}) exp(-a e^u),  a = pi Q / sqrt5
DH   : xi_DH(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) L_DH(s)   (entire, no pole);  a(n) periodic mod 5: 0,1,k,-k,-1
       Phi_DH(u) = 2 sum_n a(n) n e^{3u/2} exp(-pi n^2 e^{2u}/5)
FTH  : xi(s-th) xi(s+th):  Phi = (Phi_Z e^{-th u}) * (Phi_Z e^{th u})     (convolution of two positive functions)
FTAU : xi(s-i tau) xi(s+i tau): Phi(v) = int Phi_Z(u) Phi_Z(v-u) cos(tau (v-2u)) du
BEUR2: no theta series, no functional equation of its own -> no law (N/A).
"""
import math
import mpmath as mp

KAP = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)

def kappa():
    s5 = mp.sqrt(5); return (mp.sqrt(10 - 2*s5) - 2)/(s5 - 1)

def dh_a(n):
    k = kappa(); return [mp.mpf(0), mp.mpf(1), k, -k, mp.mpf(-1)][n % 5]

def eps_reps(Qmax):
    """r_Q = #{(m,n) in Z^2: m^2+5n^2 = Q}, Q=1..Qmax"""
    r = [0]*(Qmax + 1)
    M = int(math.isqrt(Qmax)) + 1
    for m in range(-M, M + 1):
        for n in range(-M, M + 1):
            Q = m*m + 5*n*n
            if 0 < Q <= Qmax: r[Q] += 1
    return r

class Model:
    def __init__(self, name, cutoff_log=None):
        self.name = name
        dps = mp.mp.dps
        self.cut = cutoff_log or (dps*2.4 + 40)      # stop when exponent > cut
        if name == 'EPS':
            a1 = mp.pi/mp.sqrt(5)
            Qmax = int(float(self.cut/a1)) + 2
            self.r = eps_reps(Qmax); self.Qmax = Qmax
        if name == 'DH':
            self.a = [dh_a(n) for n in range(5)]

    def phi(self, u):
        u = mp.mpf(u); cut = self.cut
        if self.name == 'ZETA':
            E2 = mp.exp(2*u); s = mp.mpf(0); n = 1
            while mp.pi*(n*n-1)*E2 < cut:
                s += (4*mp.pi**2*n**4*mp.exp(4.5*u) - 6*mp.pi*n*n*mp.exp(2.5*u))*mp.exp(-mp.pi*n*n*E2); n += 1
            return s
        if self.name == 'EPS':
            E = mp.exp(u); s = mp.mpf(0); c = mp.pi/mp.sqrt(5)
            for Q in range(1, self.Qmax + 1):
                if self.r[Q] == 0: continue
                a = c*Q
                if (a - c)*E > cut: break
                s += self.r[Q]*(a*a*mp.exp(2.5*u) - 2*a*mp.exp(1.5*u))*mp.exp(-a*E)
            return s
        if self.name == 'DH':
            E2 = mp.exp(2*u); s = mp.mpf(0); n = 1
            while mp.pi*(n*n-1)*E2/5 < cut:
                an = self.a[n % 5]
                if an != 0: s += an*n*mp.exp(-mp.pi*n*n*E2/5)
                n += 1
            return 2*mp.exp(1.5*u)*s
        raise ValueError(self.name)

U_MAX = dict(ZETA=3.6, EPS=7.6, DH=4.6)

def grid(model, h, U=None):
    """even function: values on u_j = j h, j = 0..J"""
    U = U or U_MAX[model.name]
    J = int(U/h) + 1
    return [mp.mpf(j)*h for j in range(J + 1)], [model.phi(mp.mpf(j)*h) for j in range(J + 1)]

def trap_even(f0, fs, h):
    """int_{-inf}^{inf} f for even f, f0=f(0), fs = [f(h), f(2h), ...]"""
    return h*(f0 + 2*mp.fsum(fs))

def moments_even(us, ph, h, K):
    """m_{2k} = int Phi u^{2k} du, k = 0..K (Phi even)"""
    out = []
    pw = [mp.mpf(1)]*len(us)
    u2 = [u*u for u in us]
    for k in range(K + 1):
        out.append(trap_even(ph[0]*(1 if k == 0 else 0), [pw[j]*ph[j] for j in range(1, len(us))], h))
        pw = [pw[j]*u2[j] for j in range(len(us))]
    return out

def tilted_moments(us, ph, h, N, lam):
    """M_j = int Phi(u) e^{lam u} u^j du, j=0..N, over the full line (Phi even, lam real or complex)"""
    out = []
    ep = [mp.exp(lam*u) for u in us]; em = [mp.exp(-lam*u) for u in us]
    for j in range(N + 1):
        sgn = (-1)**j
        f = [ph[i]*(us[i]**j)*(ep[i] + sgn*em[i]) for i in range(1, len(us))]
        out.append(h*((ph[0] if j == 0 else 0) + mp.fsum(f)))
    return out

def product_moments(Ma, Mb, K):
    """even moments of the convolution (Phi_a * Phi_b): m_{2k} = sum_j C(2k,j) Ma_j Mb_{2k-j}"""
    return [mp.fsum(mp.binomial(2*k, j)*Ma[j]*Mb[2*k - j] for j in range(2*k + 1)) for k in range(K + 1)]

def Xi_heat(us, ph, h, z, t=0, deriv=0):
    """H_t^{(deriv)}(z) = d^deriv/dz^deriv int e^{t u^2} Phi(u) cos(z u) du (Phi even)"""
    def f(u, p):
        w = p*mp.exp(t*u*u)
        if deriv == 0: return w*mp.cos(z*u)
        if deriv == 1: return -w*u*mp.sin(z*u)
        if deriv == 2: return -w*u*u*mp.cos(z*u)
    return trap_even(f(us[0], ph[0]), [f(us[j], ph[j]) for j in range(1, len(us))], h)
