# Grok S60: (a) (log xi)''(sigma) = 4 Var_sigma(u) under the tilted Polya measure Phi(u)e^{(2sigma-1)u}du  (log-convexity, unconditional)
# (b) symmetric off-line control xi(w+d)xi(w-d) (satisfies the functional equation, zeros at rho+-d): Gram PD by Schur product
import mpmath as mp, json
mp.mp.dps = 40
def xi(s):
    if s == 0 or s == 1: return mp.mpf(1)/2
    return mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def Phi(u):
    u = abs(mp.mpf(u)); s = 0; n = 1
    while True:
        t = (2*mp.pi**2*n**4*mp.e**(9*u)-3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u)); s += t
        if n > 3 and abs(t) < mp.mpf(10)**(-50): return s
        n += 1
out = {}
nodes = [-3,-1,-0.5,0,0.5,1,3]
for sig in [mp.mpf(1), mp.mpf('0.5'), mp.mpf(2)]:
    m = [mp.quad(lambda u: u**k*Phi(u)*mp.e**((2*sig-1)*u), nodes) for k in range(3)]
    var = m[2]/m[0]-(m[1]/m[0])**2
    out[f'4Var at sigma={sig}'] = mp.nstr(4*var, 20); out[f'(logxi)\'\' at {sig}'] = mp.nstr(mp.diff(lambda z: mp.log(xi(z)), sig, 2), 20)
for d in ['0.05', '0.1', '0.25']:
    d = mp.mpf(d); f = lambda z: mp.log(xi(z+d)) + mp.log(xi(z-d))
    out[f'sym control d={d}: (log xi(s+d)xi(s-d))\'\'(1)'] = mp.nstr(mp.diff(f, mp.mpf(1), 2), 10)
pts = [mp.mpc(a,b) for a,b in [('0.2','0'),('0.2','9'),('0.5','3'),('0.5','17'),('0.8','6'),('0.8','-11'),('1.5','1'),('1.5','25')]]
for d in ['0.3']:
    d = mp.mpf(d); F = lambda w: xi(w+d)*xi(w-d)
    G = mp.matrix([[F((mp.conj(s)+t)/2) for t in pts] for s in pts])
    E = sorted(mp.re(e) for e in mp.eighe(G, eigvals_only=True))
    out[f'Gram sym control d={d}'] = [mp.nstr(e, 5) for e in E]
print(json.dumps(out, indent=1))
