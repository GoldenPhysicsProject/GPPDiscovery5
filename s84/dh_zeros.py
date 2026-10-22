"""S84: zeros of the Davenport-Heilbronn function f(s)=((1-ik)/2)L(s,chi)+((1+ik)/2)L(s,chibar), chi mod 5, chi(2)=i,
k=(sqrt(10-2sqrt5)-2)/(sqrt5-1) (same coefficients a(n)=Re chi(n)+k Im chi(n) as s75 dhritz). NUMERICAL (mpmath):
on-line zeros from sign changes of the real Hardy function Z(t)=e^{i theta(t)} f(1/2+it), refined by findroot;
off-line zeros refined from literature approximations (Spira / Balanzario-Sanchez-Ortiz)."""
import mpmath as mp, json, sys
mp.mp.dps = 30
k = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def F(s): return (1-1j*k)/2*mp.dirichlet(s, chi) + (1+1j*k)/2*mp.dirichlet(s, chib)
def theta(t):  # odd character: gamma factor (5/pi)^{s/2} Gamma((s+1)/2)
    s = mp.mpf(0.5)+1j*t
    return mp.im(mp.loggamma((s+1)/2)) + t/2*mp.log(5/mp.pi)
ph0 = None
def Z(t):
    return mp.exp(1j*theta(t))*F(mp.mpf(0.5)+1j*t)
if __name__ == '__main__':
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 250
    z0 = Z(mp.mpf(1)); print("phase check Z(1)", z0, "Z(7.3)", Z(mp.mpf(7.3)))
    rot = mp.conj(z0)/abs(z0) if abs(mp.im(z0)) > 1e-10*abs(z0) else 1
    Zr = lambda t: mp.re(rot*Z(t))
    print("imag residue", mp.im(rot*Z(mp.mpf(33.3))), mp.re(rot*Z(mp.mpf(33.3))))
    h = 0.1; ts = [i*h for i in range(1, int(T/h))]; vals = [Zr(mp.mpf(t)) for t in ts]
    on = []
    mp.mp.dps = 45
    for i in range(len(ts)-1):
        if vals[i]*vals[i+1] < 0:
            on.append(mp.findroot(lambda t: mp.re(rot*Z(t)), (mp.mpf(ts[i]), mp.mpf(ts[i+1])), solver='anderson'))
    off = []
    for g0 in [0.808517+85.699348j, 0.650830+114.163343j, 0.574356+166.479306j, 0.724258+176.702461j]:
        r = mp.findroot(F, mp.mpc(g0)); off.append(r)
    json.dump({"T": T, "online": [mp.nstr(x, 40) for x in on], "offline": [[mp.nstr(mp.re(x), 40), mp.nstr(mp.im(x), 40)] for x in off],
               "online_count": len(on)}, open("dh_zeros_T%d.json" % int(T), "w"), indent=0)
    print(len(on), [mp.nstr(x, 12) for x in on[:8]], [mp.nstr(x, 12) for x in off])
    Nexp = T/(2*mp.pi)*mp.log(5*T/(2*mp.pi*mp.e)) + mp.mpf(1)/4  # rough RvM for q=5 (zeros with 0<t<T)
    print("RvM rough N(T)", Nexp)
