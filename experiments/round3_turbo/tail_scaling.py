# NUMERICAL scaling test of the unit-peak TailLemma (Heart, Round 3). 400 zeta zeros, test-only.
# Far set: zeros with |w-q|>=delta and |w-tau q|>=delta. Checks pointwise bound and eps^2 / eps^4 scaling.
import mpmath as mp, json
mp.mp.dps = 25
g=[mp.im(mp.zetazero(k)) for k in range(1,401)]
Z=[mp.mpf(1)/2+1j*x for x in g]+[mp.mpf(1)/2-1j*x for x in g]
tau=lambda w:1-mp.conj(w)
Ah=lambda q,w,e:e**2/(e**2-(w-q)**2)
out={"label":"NUMERICAL","rows":[]}
for name,Zs,q in [("zeta_rho1",Z,Z[0]),("Ftheta_s_theta",[z+0.25 for z in Z]+[z-0.25 for z in Z],Z[0]+0.25)]:
  for delta in [mp.mpf("0.5"),mp.mpf("0.25")]:
    far=[w for w in Zs if abs(w-q)>=delta and abs(w-tau(q))>=delta]
    for e in ["0.12","0.06","0.03"]:
      e=mp.mpf(e)
      if 2*e>delta: continue
      lin=mp.fsum(abs(Ah(q,w,e)) for w in far)
      gram=mp.fsum(abs(Ah(q,w,e))*abs(Ah(q,tau(w),e)) for w in far)
      ok=all(abs(Ah(q,w,e))<=mp.mpf(4)/3*e**2/abs(w-q)**2 for w in far)
      r=dict(case=name,delta=str(delta),eps=str(e),lin_tail=mp.nstr(lin,6),lin_over_eps2=mp.nstr(lin/e**2,6),
             gram_tail=mp.nstr(gram,6),gram_over_eps4=mp.nstr(gram/e**4,6),pointwise_bound_holds=ok)
      out["rows"].append(r);print(r)
json.dump(out,open("tail_scaling_results.json","w"),indent=1)
