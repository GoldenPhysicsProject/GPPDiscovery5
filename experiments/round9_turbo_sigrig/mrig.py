"""(M-rig) Rigorous upper bound, J=73, of the mass majorant (CHECKED, Cauchy-Schwarz):
M(T) = (2/l) int_0^T sin^2(t l) sum_{j>=J} w_j^2/(w_j^2-t^2)^2 dt,  w_j = j pi/l, l=1/4, T < w_J.
Terms j=J..N by acb.integral; tail j>N bounded by w_j^2/(w_j^2-t^2)^2 <= w_j^-2/(1-T^2/w_{N+1}^2)^2, sum_{j>N} j^-2 < 1/N."""
import json
from flint import arb, acb, ctx
ctx.prec = 80
l = arb(1)/4; pi = arb.pi(); J = 73; N = 600
def M(T):
    T = arb(T); tot = arb(0)
    for j in range(J, N+1):
        w2 = (j*pi/l)**2
        tot += acb.integral(lambda t,an: (t*l).sin()**2*w2/(w2-t*t)**2, 0, T).real
    wN1 = (N+1)*pi/l
    tail = (l/pi)**2/N/(1-T**2/wN1**2)**2 * (T/2 - (2*T*l).sin()/(4*l))
    return (2/l)*(tot + tail)
out = {}
for name, T in [('31.05', arb(3105)/100), ('omega13', 13*pi/l)]:
    m = M(T); out['M_upper_'+name] = float(m.upper()); out['M_ball_'+name] = m.str(12)
s = json.load(open('sigrig_results.json'))
m0 = arb(out['M_upper_31.05']); m13 = arb(out['M_upper_omega13'])
lam = (arb(s['min_lower_t_ge_omega13'])*(pi/2 - m13) + arb(s['min_lower_on_0_T0'])*m0)/pi
out['lamZ_lower'] = float(lam.lower()); out['N'] = N
print(json.dumps(out, indent=1)); json.dump(out, open('mrig_results.json','w'), indent=1)
