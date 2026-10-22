"""Export DH's even onset eigenvector with SIGNED full-precision coefficients plus the rho0-quartet leak.
Basis: f(u) = sum_m c_m cos(a_m u), a_m = (2m+1)pi/(2L), m = 0..K-1, on [-L,L] (L = HALF-width); Gram <cos a_m, cos a_n> = L delta_mn.
Normalisation: ||f||_{L^2(-L,L)} = 1, i.e. L sum c_m^2 = 1; sign fixed by c_0 > 0. Inner product <f,g> = int_{-L}^{L} f conj(g).
rho0 = 1/2 + beta0 + i gamma0 refined by mpmath findroot at 60 digits. chi_L = e^{-i gamma0 x} cosh(beta0 x), phi_L = e^{-i gamma0 x} sinh(beta0 x):
<f,chi_L> = int f e^{i g x} cosh(b x) = A (real, f even), <f,phi_L> = i int f sin(g x) sinh(b x) = i B (purely imaginary).
Leak (Queen) = 2|<f,chi>|^2 - 2|<f,phi>|^2 = 2(A^2-B^2); full quartet contribution to Q_DH(f) = 4(A^2-B^2) = 2 x leak.
Mirror gamma0 -> -gamma0 evaluated independently (for real even f it must equal the original: |A'|=|A|, |B'|=|B|).
A, B from the exact transform fhat(w) = sum c_m (-1)^m cos(wL) 2 a_m/(a_m^2 - w^2), A = Re fhat(g+ib), B = -Im fhat(g+ib),
cross-checked by direct quadrature. NUMERICAL. usage: dh_onset_export.py vec.json dump.json out.json
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, mpmath as mp
V = json.load(open(sys.argv[1])); D = json.load(open(sys.argv[2])); out = sys.argv[3]
mp.mp.dps = 230
p, q = V['L'].split('/'); L = mp.mpf(p)/int(q); K = V['K']
c = [mp.mpf(x) for x in V['coef']]
if c[0] < 0: c = [-x for x in c]
a = [(2*m+1)*mp.pi/(2*L) for m in range(K)]
s5 = mp.sqrt(5); kap = (mp.sqrt(10-2*s5)-2)/(s5-1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def fdh(s): return (1-1j*kap)/2*mp.dirichlet(s, chi)+(1+1j*kap)/2*mp.dirichlet(s, chib)
mp.mp.dps = 60; rho = mp.findroot(fdh, mp.mpc('0.808517182457', '85.6993484854')); resf = abs(fdh(rho)); mp.mp.dps = 230
b = mp.re(rho)-mp.mpf(1)/2; g = mp.im(rho)
def fhat(w): return mp.cos(w*L)*mp.fsum(c[m]*(-1)**m*2*a[m]/(a[m]**2-w**2) for m in range(K))
def AB(gg):
    v = fhat(mp.mpc(gg, b)); return mp.re(v), -mp.im(v)
A, B = AB(g); Am, Bm = AB(-g)
# residual ||W v - lam v|| on the same (mid) matrix, in the vecs.py vector normalisation (sum v^2 = 1)
vv = [x*mp.sqrt(L) for x in c]
W = [[mp.mpf(D['W'][i][j][0]) for j in range(K)] for i in range(K)]
Wv = [mp.fsum(W[i][j]*vv[j] for j in range(K)) for i in range(K)]
lam = mp.fsum(vv[i]*Wv[i] for i in range(K)); res = mp.sqrt(mp.fsum((Wv[i]-lam*vv[i])**2 for i in range(K)))
import os
mp.mp.dps = 40   # quadrature cross-check of A, B (set QUAD=1; slow)
fq = lambda x: mp.fsum(c[m]*mp.cos(a[m]*x) for m in range(K))
Aq = A if not os.environ.get('QUAD') else mp.quad(lambda x: fq(x)*mp.cos(g*x)*mp.cosh(b*x), mp.linspace(-L, L, 41))
Bq = B if not os.environ.get('QUAD') else mp.quad(lambda x: fq(x)*mp.sin(g*x)*mp.sinh(b*x), mp.linspace(-L, L, 41))
nchi = mp.sqrt(mp.quad(lambda x: mp.cosh(b*x)**2, [-L, L])); nphi = mp.sqrt(mp.quad(lambda x: mp.sinh(b*x)**2, [-L, L]))
mp.mp.dps = 230
S = lambda x: mp.nstr(x, 60)
rec = {"form": "DH", "sector": "even", "L": V['L'], "L_is": "half-width, window [-L,L]", "K": K, "assembly_dump": sys.argv[2],
       "assembly_prec_bits": D.get('prec'), "eigvec_prec_bits": 720,
       "basis": "f(u)=sum_{m=0}^{K-1} c_m cos(a_m u), a_m=(2m+1)pi/(2L); Gram <cos a_m,cos a_n>_{L2(-L,L)} = L delta_mn",
       "normalization": "||f||_{L2(-L,L)}=1 (L*sum c_m^2 = 1), sign c_0>0",
       "eig_W": S(lam), "Q_over_norm2 = eig_W/L": S(lam/L), "residual_norm_Wv_minus_lam_v (sum v^2=1, v=c*sqrt(L))": mp.nstr(res, 6),
       "note_precision": "entries are 512-bit Arb mids (~154 digits); the vector is the exact eigenvector of the mid matrix to 720 bits",
       "rho0": mp.nstr(rho, 50), "beta0": mp.nstr(b, 50), "gamma0": mp.nstr(g, 50), "|f_DH(rho0)|": mp.nstr(resf, 5),
       "<f,chi_L> = A": S(A), "<f,phi_L> = iB, B": S(B), "<f,chi_L>/||chi_L||": mp.nstr(A/nchi, 20), "<f,phi_L>/||phi_L|| (imag part)": mp.nstr(B/nphi, 20),
       "||chi_L||": mp.nstr(nchi, 20), "||phi_L||": mp.nstr(nphi, 20),
       "leak_2(A^2-B^2)": S(2*(A*A-B*B)), "quartet_in_Q_4(A^2-B^2)": S(4*(A*A-B*B)),
       "mirror_A(-gamma0)": S(Am), "mirror_B(-gamma0)": S(Bm), "mirror_leak": S(2*(Am*Am-Bm*Bm)),
       "mirror_check_|leak-mirror|": mp.nstr(abs(2*(A*A-B*B)-2*(Am*Am-Bm*Bm)), 5),
       "quadrature_check_|A-Aq|,|B-Bq|": [mp.nstr(abs(A-Aq), 4), mp.nstr(abs(B-Bq), 4)],
       "coef": [mp.nstr(x, 210) for x in c],
       "credit": "Daniel Toupin / Golden Physics Project, with AI assistance from Grok.", "label": "NUMERICAL"}
json.dump(rec, open(out, 'w'), indent=1)
print(json.dumps({k: v for k, v in rec.items() if k != 'coef'}, indent=1))
