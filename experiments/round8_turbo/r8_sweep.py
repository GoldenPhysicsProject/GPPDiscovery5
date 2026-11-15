"""Round 8 (T2) sweep: omega^4 rho_j for j = 13..1012 step 3 (T = 2e5 + mean tail), envelope vs log omega (NUMERICAL)."""
import numpy as np, json
exec(open('r8_t2.py').read().split("out = dict")[0].replace("JJ = [73, 123, 173, 273, 472]", "JJ = list(range(13, 1013, 3))"))
q, r4, v = run(2e5, chunk=20000)
oms = np.array([om(j) for j in JJ]); mx = np.abs(r4).max(0)
out = dict(J=JJ, omega4_rho_max=mx.tolist(), omega4_rho_delta=r4[N].tolist(), omega4_rho_phirows_max=np.abs(r4[:N]).max(0).tolist(),
           omega3_v_max=(np.abs(v).max(0)*oms**3).tolist(), q=q.tolist())
for a,b in [(13,72),(73,172),(173,372),(373,672),(673,1012)]:
    sel = [(k) for k,j in enumerate(JJ) if a <= j <= b]
    print(a,b,'max w4|rho|', round(float(mx[sel].max()),3), 'max/log(om)', round(float((mx[sel]/np.log(oms[sel])).max()),3),
          'delta-row max', round(float(np.abs(r4[N,sel]).max()),3), 'phi-rows max', round(float(np.abs(r4[:N][:,sel]).max()),3), 'max w3|v|', round(float((np.abs(v[:,sel]).max(0)*oms[sel]**3).max()),4))
json.dump(out, open('r8_sweep_results.json','w'))
