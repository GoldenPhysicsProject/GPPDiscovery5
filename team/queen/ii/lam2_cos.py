"""Queen S1d (ii): lambda_2(R_e) and lambda_1(R_o) of zeta's window form in Goddess of Prime Number Theory's Arb cosine
assembly (team/prime-number-theory/precision/cosritz_dump.py = Grok S62 cosritz.py, unchanged). Credit: Goddess of Prime
Number Theory / team/prime-number-theory, assembly from Grok S62.
Basis even c_m = cos((2m+1) pi u/(2L)), odd s_m = sin(m pi u/L), m>=1; Gram = L*I on L^2[-L,L] -> Rayleigh quotient = eig/L.
Pole (as in the dump): even +2 cp cp^T, cp_m = 2 a_m (-1)^m cosh(L/2)/(a_m^2+1/4); odd -2 sp sp^T,
sp_m = 2(-b_m (-1)^(m+1) sinh(L/2))/(b_m^2+1/4). R_e = W_e - 2 cp cp^T, R_o = W_o + 2 sp sp^T (= Q - pole, our convention).
Eigenvalue nearest 0 by inverse iteration on midpoints (arb_mat solve at prec P), Rayleigh quotient in Arb;
Weyl perturbation bound from radii: |delta lambda| <= ||Rad||_2/L <= k*max|rad|/L (float underflow-safe). Nested k-sections = leading k x k blocks. NUMERICAL."""
import sys, json, time
from flint import arb, arb_mat, ctx
def run(path, ks, P, iters=12):
    D = json.load(open(path)); K = D["K"]; ctx.prec = P
    Ls = D["L"]; L = arb(int(Ls.split('/')[0]))/int(Ls.split('/')[1]) if '/' in Ls else arb(Ls)
    PI = arb.pi(); sec = D["sector"]
    if sec == "even":
        a = [(2*m+1)*PI/(2*L) for m in range(K)]
        pv = [2*a[m]*(-1)**m*(L/2).cosh()/(a[m]**2+arb(1)/4) for m in range(K)]; sg = -2
    else:
        b = [m*PI/L for m in range(1, K+1)]
        pv = [2*(-b[m]*(-1)**(m+1)*(L/2).sinh())/(b[m]**2+arb(1)/4) for m in range(K)]; sg = +2
    out = []
    for k in ks:
        if k > K: continue
        t0 = time.time()
        R = arb_mat(k, k); rmax = 0.0
        for i in range(k):
            for j in range(k):
                R[i, j] = arb(D["W"][i][j][0]) + sg*pv[i]*pv[j]
                rmax = max(rmax, D["W"][i][j][1])
        Rm = arb_mat(k, k)
        for i in range(k):
            for j in range(k): Rm[i, j] = R[i, j].mid()
        x = arb_mat(k, 1)
        for i in range(k): x[i, 0] = arb(1)/(i+1)
        lam_prev = None
        for it in range(iters):
            y = Rm.solve(x)
            nrm = sum((y[i, 0].mid()**2 for i in range(k)), arb(0)).sqrt()
            for i in range(k): x[i, 0] = (y[i, 0]/nrm).mid()
        Rx = R*x
        num = sum((x[i, 0]*Rx[i, 0] for i in range(k)), arb(0)); den = sum((x[i, 0]**2 for i in range(k)), arb(0))
        lam = num/den/L
        # residual ||R x - lam L x|| / L  bounds distance to spectrum (for symmetric R)
        res = sum(((Rx[i, 0]-lam*L*x[i, 0]).mid()**2 for i in range(k)), arb(0)).sqrt()/L/den.sqrt()
        out.append(dict(L=Ls, sector=sec, k=k, prec=P, lam=lam.mid().str(20, radius=False), lam_rad=float(lam.rad()),
                        residual=float(res.mid()), weyl=k*rmax/float(L.mid()),
                        vec=[x[i, 0].mid().str(40, radius=False) for i in range(k)], secs=round(time.time()-t0, 1)))
        print({kk: v for kk, v in out[-1].items() if kk != "vec"}, flush=True)
    return out
if __name__ == "__main__":
    res = run(sys.argv[1], [int(x) for x in sys.argv[2].split(",")], int(sys.argv[3]))
    json.dump(res, open(sys.argv[4], "w"))
