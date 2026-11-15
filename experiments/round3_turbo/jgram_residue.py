# NUMERICAL test (Heart, Round 3, Turbo 6889016436 item ii).
# Zeta-side evaluation of the J-Gram with NO zero data: for rational R decaying like |w|^-4,
#   sum_w m_w R(w) = - sum_{poles p of R} Res_p( R * L ),  L = F'/F  (F = xi, or F_theta(w)=xi(w-th)xi(w+th)).
# Abel kernel Ahat^q(w) = eps^2/(eps^2-(w-q)^2); conj(Ahat^q(tau w)) = Ahat^{tau q}(w), tau q = 1-conj(q).
#   B = W(g,Jg) = sum m Ahat^q(w)^2
#     = -[ (eps^2/4)(L'(q+eps)+L'(q-eps)) + (eps/4)(L(q-eps)-L(q+eps)) ]
#   A = W(g,g)  = sum m Ahat^q(w) Ahat^{tau q}(w)
#     = -[ -(eps/2)Ahat^{tq}(q+eps)L(q+eps) + (eps/2)Ahat^{tq}(q-eps)L(q-eps)
#          -(eps/2)Ahat^{q}(tq+eps)L(tq+eps) + (eps/2)Ahat^{q}(tq-eps)L(tq-eps) ]   (tq != q)
#   A = B when tq = q (on-line label).
# L is evaluated from zeta's standard continuation (mpmath zeta, derivative=1) plus Gamma/pole terms.
# The 400-zero sums in ../round2_flint/jgram_results.json are used ONLY as the comparison target.
import mpmath as mp, json, time
mp.mp.dps = 40
t0 = time.time()
def Lxi(w):   # xi'/xi
    return mp.zeta(w, derivative=1)/mp.zeta(w) + 1/w + 1/(w-1) - mp.log(mp.pi)/2 + mp.digamma(w/2)/2
def dLxi(w):
    return mp.diff(Lxi, w)
th = mp.mpf(1)/4
models = {
  "zeta_rho1": (Lxi, dLxi, mp.mpf(1)/2 + 1j*mp.im(mp.zetazero(1))),   # label q = rho_1 (zero used only as LABEL)
  "Ftheta_s_theta": (lambda w: Lxi(w-th)+Lxi(w+th), lambda w: dLxi(w-th)+dLxi(w+th),
                     mp.mpf(1)/2 + 1j*mp.im(mp.zetazero(1)) + th),
}
Ah = lambda q,w,e: e**2/(e**2-(w-q)**2)
filed = {(r["case"], str(r["eps"])): r for r in json.load(open("../round2_flint/jgram_results.json"))["rows"]}
out = {"label": "NUMERICAL", "method": "zeta-side residue identity, no zero data", "rows": []}
for name,(L,dL,q) in models.items():
    tq = 1 - mp.conj(q)
    for e in ["0.3","0.2","0.12","0.06"]:
        ee = mp.mpf(e)
        B = -((ee**2/4)*(dL(q+ee)+dL(q-ee)) + (ee/4)*(L(q-ee)-L(q+ee)))
        if abs(tq-q) < mp.mpf(10)**-30:
            A = B
        else:
            A = -(-(ee/2)*Ah(tq,q+ee,ee)*L(q+ee) + (ee/2)*Ah(tq,q-ee,ee)*L(q-ee)
                  -(ee/2)*Ah(q,tq+ee,ee)*L(tq+ee) + (ee/2)*Ah(q,tq-ee,ee)*L(tq-ee))
        D = 2*mp.re(A) - 2*mp.re(B)
        f = filed[(name, e)]
        Af = mp.mpc(complex(f["A"].strip("()").replace(" ",""))); Bf = mp.mpc(complex(f["B"].strip("()").replace(" ","")))
        row = dict(case=name, eps=e, A_residue=mp.nstr(A,13), B_residue=mp.nstr(B,13), D_residue=mp.nstr(D,13),
                   det_residue=mp.nstr(A*A-abs(B)**2,13),
                   A_zero_sum_400=f["A"], B_zero_sum_400=f["B"],
                   absdiff_A=mp.nstr(abs(A-Af),3), absdiff_B=mp.nstr(abs(B-Bf),3))
        out["rows"].append(row); print(row)
out["seconds"] = round(time.time()-t0,1)
json.dump(out, open("jgram_residue_results.json","w"), indent=1)
