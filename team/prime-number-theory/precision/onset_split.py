"""Three-way explicit-formula split of the Davenport-Heilbronn Weil form on a normalised even vector f (||f||=1 on
[-L,L], L = HALF-width), f(u) = sum_m c_m cos(a_m u), a_m=(2m+1)pi/(2L):
  Q_DH(f) = (a) sum over on-line zeros 1/2+i gamma, |gamma|<=T, of fhat(gamma)^2  [+ smooth tail beyond T]
          + (b) quartet of rho0 = 1/2+beta0+i gamma0 = 0.808517...+85.699348 i : 4(A^2 - B^2)
          + (c) everything else = Q - (a) - (b)  (other off-line quartets, of which the three next known ones are
                listed separately, plus the error of the T-truncation / smooth tail).
fhat(w) = int f(u) cos(wu) du = sum_m c_m (-1)^m cos(wL) 2a_m/(a_m^2-w^2) (exact). A = Re fhat(gamma0+i beta0) =
<f, cos(gamma0 x)cosh(beta0 x)>, B = -Im fhat = <f, sin(gamma0 x) sinh(beta0 x)>. Queen's normalisation
2|<f,chi_L>|^2 - 2|<f,phi_L>|^2 = 2(A^2-B^2) = (b)/2, chi_L = e^{-i gamma0 x}cosh(beta0 x), phi_L = e^{-i gamma0 x} sinh(beta0 x).
Tail: int_T^inf fhat(t)^2 (1/pi) log(5t/2pi) dt (both signs of gamma; Riemann-von Mangoldt density for conductor 5).
Off-line zeros refined by mpmath findroot from the S-team values. NUMERICAL throughout.
usage: onset_split.py vec.json zeros.json[,zeros2.json..] T out.json   (zeros files from dh_zeros.py)
       onset_split.py --basis Wdump.json zeros.json T   (validation: W_00/L of basis mode 0 vs the zero sum)
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, mpmath as mp
mp.mp.dps = 60
s5 = mp.sqrt(5); kap = (mp.sqrt(10-2*s5)-2)/(s5-1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def fdh(s): return (1-1j*kap)/2*mp.dirichlet(s, chi)+(1+1j*kap)/2*mp.dirichlet(s, chib)
OFF0 = ['0.808517182457+85.6993484854j', '0.65083008061+114.163342731j', '0.574356050451+166.479305913j', '0.724257694627+176.702461243j']
def offline():
    mp.mp.dps = 30; r = [mp.findroot(fdh, mp.mpc(z)) for z in OFF0]; mp.mp.dps = 60; return r
def loadzeros(spec, T):
    zs = []
    for fn in spec.split(','): zs += json.load(open(fn))['zeros']
    zs = sorted(zs); u = [z for i, z in enumerate(zs) if i == 0 or z-zs[i-1] > 1e-6]   # dedupe chunk overlaps; NO rounding
    return [z for z in u if z <= T]
def make_fhat(c, L):
    a = [(2*m+1)*mp.pi/(2*L) for m in range(len(c))]
    def fh(w): return mp.cos(w*L)*mp.fsum(c[m]*(-1)**m*2*a[m]/(a[m]**2-w**2) for m in range(len(c)))
    return fh
def online(fh, zs, T, L):
    S = 2*mp.fsum(fh(mp.mpf(g))**2 for g in zs)
    mp.mp.dps = 30
    tail = mp.quad(lambda t: fh(t)**2*mp.log(5*t/(2*mp.pi))/mp.pi, mp.linspace(T, T+40, 41)+[mp.inf])
    mp.mp.dps = 60
    return S, tail
def Lval(s):
    p, q = (s.split('/')+['1'])[:2]; return mp.mpf(p)/mp.mpf(q)
if sys.argv[1] == '--basis':
    D = json.load(open(sys.argv[2])); L = Lval(D['L']); T = float(sys.argv[4]); zs = loadzeros(sys.argv[3], T)
    c = [mp.mpf(1)] + [0]*0; fh = make_fhat(c, L)
    S, tail = online(fh, zs, T, L); q = [4*mp.re(fh(mp.mpc(z.imag, 0.5-z.real) if False else (z.imag - 1j*(z.real-0.5))))**0 for z in []]
    off = offline(); quart = [4*mp.re(fh(mp.mpc(z.imag, -(z.real-mp.mpf(1)/2)))**2) for z in off if z.imag <= T]
    W00 = mp.mpf(D['W'][0][0][0])
    print(json.dumps({"L": D['L'], "W00": mp.nstr(W00, 15), "n_online": len(zs), "online_sum": mp.nstr(S, 15), "tail": mp.nstr(tail, 6),
                      "offline_quartets_T": [mp.nstr(x, 10) for x in quart], "zero_side": mp.nstr(S+tail+sum(quart), 15),
                      "diff": mp.nstr(W00-S-tail-sum(quart), 6)}))
    sys.exit()
V = json.load(open(sys.argv[1])); L = Lval(V['L']); T = float(sys.argv[3])
c = [mp.mpf(x) for x in V['coef']]; fh = make_fhat(c, L); Q = mp.mpf(V['Q_over_norm2'])
zs = loadzeros(sys.argv[2], T); S, tail = online(fh, zs, T, L)
off = offline(); rows = []
for z in off:
    beta = z.real-mp.mpf(1)/2; w = mp.mpc(z.imag, beta)      # fhat(gamma0 + i beta0), conj gives the other pair
    v = fh(w); A = mp.re(v); B = -mp.im(v)
    # overlaps with the normalised probes cos(g x)cosh(b x), sin(g x)sinh(b x) on [-L,L]
    nchi = mp.sqrt(mp.quad(lambda x: (mp.cos(z.imag*x)*mp.cosh(beta*x))**2, [-L, 0, L]))
    nphi = mp.sqrt(mp.quad(lambda x: (mp.sin(z.imag*x)*mp.sinh(beta*x))**2, [-L, 0, L]))
    rows.append({"rho": mp.nstr(z, 15), "A": mp.nstr(A, 25), "B": mp.nstr(B, 25), "quartet_4(A2-B2)": mp.nstr(4*(A*A-B*B), 25),
                 "queen_2(A2-B2)": mp.nstr(2*(A*A-B*B), 25), "overlap_chi": mp.nstr(abs(A)/nchi, 8), "overlap_phi": mp.nstr(abs(B)/nphi, 8)})
q0 = 4*(mp.mpf(rows[0]["A"])**2-mp.mpf(rows[0]["B"])**2)
qrest = mp.fsum(4*(mp.mpf(r["A"])**2-mp.mpf(r["B"])**2) for r in rows[1:])
res = {"vec": sys.argv[1], "form": V['form'], "L": V['L'], "K": V['K'], "Q": V['Q_over_norm2'], "fp_L": V['fp_L'],
       "T": T, "n_online_zeros": len(zs), "a_online_sum_T": mp.nstr(S, 25), "a_tail_beyond_T": mp.nstr(tail, 8),
       "tail_model_2fp2logT_over_3piT3": mp.nstr(2*mp.mpf(V['fp_L'])**2*mp.log(5*T/(2*mp.pi))/(3*mp.pi*T**3), 6),
       "b_quartet_rho0": mp.nstr(q0, 25), "other_known_quartets": mp.nstr(qrest, 12),
       "c_remainder": mp.nstr(Q-S-tail-q0, 12), "c_minus_known_quartets": mp.nstr(Q-S-tail-q0-qrest, 12), "offline": rows}
print(json.dumps(res, indent=1)); json.dump(res, open(sys.argv[4], 'w'), indent=1)
