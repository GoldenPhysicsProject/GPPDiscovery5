# Grok S60 H60a: Polya/Riemann representation => coherent kernel xi((conj(s)+t)/2) is PD on all of C x C.
# Conventions: xi(s)=1/2 s(s-1) pi^{-s/2} Gamma(s/2) zeta(s);  Phi(u)=sum_n (2pi^2 n^4 e^{9u}-3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})
# de Bruijn-Newman: H_0(z)=int_0^oo Phi(u)cos(zu)du=(1/8)xi(1/2+iz/2)  =>  xi(1/2+w)=4 int_R Phi(u) e^{2wu} du.
# Hence K(s,t):=xi((conj s+t)/2)=<e_s,e_t>_{L^2(R,4Phi du)}, e_s(u)=e^{(s-1/2)u}.
import mpmath as mp, json, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
def xi(s): return mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def Phin(n,u): return (2*mp.pi**2*n**4*mp.e**(9*u)-3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u))
def PhiRaw(u):
    u = mp.mpf(u); s = 0; n = 1
    while True:
        t = Phin(n,u); s += t
        if n > 3 and abs(t) < mp.mpf(10)**(-mp.mp.dps-10): return s
        n += 1
def Phi(u): return PhiRaw(abs(mp.mpf(u)))   # evenness (theta functional eq.) checked below on PhiRaw
out = {}
# (0) Phi even, positive; per-n pieces signed
out['Phi_even_check'] = [mp.nstr(PhiRaw(mp.mpf(u))-PhiRaw(-mp.mpf(u)),5) for u in ['0.05','0.2','0.5']]
out['Phi_min_on_grid'] = mp.nstr(min(Phi(mp.mpf(k)/40) for k in range(-60,61)),8)
out['Phi_1(-0.3)'] = mp.nstr(Phin(1,mp.mpf('-0.3')),8)   # negative: Dirichlet sector n=1 is not a positive weight
# (1) Fourier/Laplace identity at complex w (inside and outside the strip)
def lap(wp, weight=lambda u: 1):
    f = lambda u: 4*Phi(u)*weight(u)*mp.e**(2*wp*u)
    return mp.quad(f, [-3,-1,-0.5,0,0.5,1,3])
chk = {}
for wp in [mp.mpc(0,0), mp.mpc('0.3','2.0'), mp.mpc('-0.7','5'), mp.mpc('1.2','-3'), mp.mpc('0','7.06736')]:
    chk[str(wp)] = mp.nstr(abs(lap(wp)-xi(mp.mpf(1)/2+wp))/abs(xi(mp.mpf(1)/2+wp)),5)
out['laplace_rel_err'] = chk
# (2) Gram matrices
pts = [mp.mpc(a,b) for a,b in [('0.2','0'),('0.2','9'),('0.5','3'),('0.5','17'),('0.8','6'),('0.8','-11'),('1.5','1'),('1.5','25')]]
def gram(F): return mp.matrix([[F((mp.conj(s)+t)/2) for t in pts] for s in pts])
def eig(M):
    E = mp.eighe(M, eigvals_only=True) if hasattr(mp,'eighe') else mp.eigh(M, eigvals_only=True)
    return sorted([mp.re(e) for e in E])
def rep(name, F):
    E = eig(gram(F)); out[name] = {'min': mp.nstr(E[0],6), 'max': mp.nstr(E[-1],6), 'all': [mp.nstr(e,4) for e in E]}
rep('xi', xi)
d = mp.mpf('0.3'); rep('shift_xi_delta0.3', lambda w: xi(w+d))           # zeros at rho-0.3: off the line
rep('teeth_(w-1/2)^2_xi', lambda w: (w-mp.mpf(1)/2)**2*xi(w))            # signed measure -Phi''/4: must FAIL
rep('teeth_codex_halfLambda', lambda w: mp.pi**(-w/2)*mp.gamma(w/2)*mp.zeta(w)/2)  # Codex pre-polar kernel outside Re>1
rep('teeth_sector_n1', lambda w: w*(w-1)*mp.pi**(-w/2)*mp.gamma(w/2)/2)  # Dirichlet n=1 sector of xi kernel
# dBN control t<0 (Rodgers-Tao Lambda>=0: H_t has non-real zeros for every t<0)
tt = mp.mpf('-0.2'); w_t = lambda wp: lap(wp, lambda u: mp.e**(tt*u*u))
if mp.mp.dps <= 50:
    E = eig(mp.matrix([[w_t((mp.conj(s)+t)/2-mp.mpf(1)/2) for t in pts] for s in pts]))
    out['dBN_t-0.2'] = {'min': mp.nstr(E[0],6), 'max': mp.nstr(E[-1],6)}
# (3) diagonal norms ||e_s||^2 = xi(Re s)
out['norm_check'] = {str(s): mp.nstr(abs(lap((mp.conj(s)+s)/2-mp.mpf(1)/2)-xi(mp.re(s))),5) for s in pts[:4]}
# (4) orthogonality at a zero: s=0.2+iy, t=0.8+i(y+2*gamma1)
g1 = mp.im(mp.zetazero(1)); s = mp.mpc('0.2','1'); t = mp.mpc('0.8',1+2*g1)
out['orth_at_rho1'] = mp.nstr(abs(xi((mp.conj(s)+t)/2)),5)
# (5) Codex rigging measure: 8 e^{-u} Phi(u) = (D^2+2D) psi(e^{4u}),  psi(x)=sum exp(-pi n^2 x)
psi = lambda u: mp.nsum(lambda n: mp.e**(-mp.pi*n**2*mp.e**(4*u)), [1, mp.inf])
out['D2+2D_check'] = [mp.nstr(mp.diff(psi,u,2)+2*mp.diff(psi,u,1)-8*mp.e**(-u)*Phi(u),5) for u in [mp.mpf('-0.3'),mp.mpf('0.1')]]
print(json.dumps(out, indent=1))
