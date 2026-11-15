# NUMERICAL test of derivation steps (Heart, Round 2 reply to Flint 6888778018).
# Zeros used ONLY to test identities/limits, never as proof inputs.
# Unit-peak Abel kernel Ahat^q(w) = eps^2/(eps^2-(w-q)^2); tau(w)=1-conj(w).
# W(f,h) = sum_w m_w fhat(w) conj(hhat(tau w))   (sesquilinear spectral face)
import mpmath as mp, json, time
mp.mp.dps = 30
N = 400
t0=time.time()
g = [mp.im(mp.zetazero(k)) for k in range(1,N+1)]
half = mp.mpf(1)/2
Zz = [half+1j*x for x in g] + [half-1j*x for x in g]           # zeta zeros (m=1)
Zt = [z+mp.mpf(1)/4 for z in Zz] + [z-mp.mpf(1)/4 for z in Zz]   # F_theta zeros, theta=1/4
tau = lambda w: 1-mp.conj(w)
A_ = lambda q,w,e: e**2/(e**2-(w-q)**2)
out = {"label":"NUMERICAL","N_zeros":N,"rows":[]}
for name,Z,s in [("zeta_rho1",Zz,Zz[0]),("Ftheta_s_theta",Zt,Zz[0]+mp.mpf(1)/4)]:
    for e in [0.30,0.20,0.12,0.06]:
        e = mp.mpf(e)
        A = mp.fsum(A_(s,w,e)*mp.conj(A_(s,tau(w),e)) for w in Z)       # W(g,g)
        B = mp.fsum(A_(s,w,e)**2 for w in Z)                             # W(g,Jg)
        sel = e**2/(e**2-(2*mp.re(s)-1)**2)                               # Flint C3 factor
        fh = lambda w: A_(s,w,e)-A_(tau(s),w,e)                            # f = g - Jg
        Wexact = mp.fsum(fh(w)*mp.conj(fh(tau(w))) for w in Z)
        Wminus_sq = -mp.fsum(fh(w)**2 for w in Z)
        Wminus_abs = -mp.fsum(abs(fh(w))**2 for w in Z)
        D_direct = 2*mp.re(A)-2*mp.re(B)
        out["rows"].append({k:mp.nstr(v,12) for k,v in dict(case=0,eps=e,A=A,B=B,A_over_sel=A/sel if sel!=0 else 0,
            W_fJf_direct=Wexact, minus_sum_fhat_sq=Wminus_sq, minus_sum_abs_sq=Wminus_abs,
            abs_gap=Wexact-Wminus_abs, D_from_Gram=D_direct, det=A*A-abs(B)**2).items()} | {"case":name})
out["seconds"]=round(time.time()-t0,1)
json.dump(out,open("jgram_results.json","w"),indent=1)
for r in out["rows"]: print(r)
