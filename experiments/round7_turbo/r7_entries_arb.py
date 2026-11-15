"""Round 7 (Turbo 6889791941) GAP (a): arb entries of the bordered Gram on {phi_1..phi_72, delta_l} at L = 0.5 (l = 1/4).
Same exact reduction as Round 5 / Round 6 (CHECKED there):
  B(phi_j,phi_k) = (2/l){(1/pi)∫_0^∞ Psi^2 V S_j S_k dt - 2 Psi(i/2)^2 sigma_j sigma_k},  S_j = w_j l sinc((t-w_j)l)/(t+w_j),
  B(phi_i,delta) = (1/pi)∫ Psi^2 V c2 S_i sin(t l) dt - 2 Psi(i/2)^2 c2 sigma_i sinh(l/2),   g = B(delta,delta),  c2 = sqrt(2/l).
Tail lemmas (CHECKED by hand in Rounds 5/6): Psi^2 <= 16/t^4, |V| <= log t + C (t >= 100), |S_j| <= (4/3) w_j/t^2 for t >= 2 w_j.
  phi-phi on [T,inf): (16*16/9/pi) w_j w_k T^-7 [(log T + C)/7 + 1/49]  (needs T >= 2 max(w_j,w_k));
  phi-delta: (16/pi)(4/3) c2 w_i T^-5 [(log T + C)/5 + 1/25];  delta-delta: (16/pi) T^-3 [(log T + C)/3 + 1/9].
Usage: python3 r7_entries_arb.py WORKER NWORKERS NTOT   -> shard_WORKER.json (entries as arb strings). No zero data."""
from flint import acb, arb, ctx
import json, time, sys, math
ctx.prec = 160
wid, nw, NT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
l = arb(1)/4; pi = arb.pi(); lp = pi.log(); c2 = (2/l).sqrt()
NS = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
W = [j*pi/l for j in range(1, NT+1)]
bI2 = 16*l.sinh()**2; sd = (l/2).sinh()
def S(j, t): w = W[j]; return w*l*((t-w)*l).sinc()/(t+w)
def V(t):
    i = acb(0,1)
    return ((acb(l)+i*t/2).digamma()+(acb(l)-i*t/2).digamma())/2 - lp - 2*sum(arb(p).log()/arb(n).sqrt()*(t*arb(n).log()).cos() for n,p in NS)
P0 = sum(arb(p).log()/arb(n).sqrt() for n,p in NS); C = 2.5 + 2*P0
sigma = [-(-1)**(j+1)*W[j]*sd/(arb(1)/4 + W[j]**2) for j in range(NT)]
def integ(f, T, step=50):
    cuts = list(range(0, T+1, step))
    I = sum((acb.integral(f, cuts[m], cuts[m+1]) for m in range(len(cuts)-1)), acb(0))
    assert abs(float(I.imag.mid())) < 1e-30
    return I.real
jobs = [(j,k) for j in range(NT) for k in range(j, NT)] + [(i,'d') for i in range(NT)] + [('d','d')]
mine = jobs[wid::nw]
out = {}; t0 = time.time()
for (a,b) in mine:
    if b == 'd' and a == 'd':
        T = 40000; TT = arb(T)
        tg = 16/pi*TT**-3*((TT.log()+C)/3+arb(1)/9)
        val = integ(lambda t, _: (t/2).sinc()**4*V(t)*(t*l).sin()**2/pi, T, 100) + arb(0, float(tg.upper())) - 2*bI2**2*sd**2
    elif b == 'd':
        T = 8000; TT = arb(T)
        tu = 16/pi*arb(4)/3*W[a]*c2*TT**-5*((TT.log()+C)/5+arb(1)/25)
        val = integ(lambda t, _: (t/2).sinc()**4*V(t)*c2*S(a,t)*(t*l).sin()/pi, T) + arb(0, float(tu.upper())) - 2*bI2**2*c2*sigma[a]*sd
    else:
        wmax = max(float(W[a].mid()), float(W[b].mid()))
        T = max(800, 50*math.ceil((2*wmax+1)/50)); TT = arb(T)
        tc = 16*arb(16)/9/pi*TT**-7*((TT.log()+C)/7+arb(1)/49)*W[a]*W[b]
        I = integ(lambda t, _: (t/2).sinc()**4*V(t)*S(a,t)*S(b,t)/pi, T)
        val = (2/l)*(I + arb(0, float(tc.upper())) - 2*bI2**2*sigma[a]*sigma[b])
    out[f"{a},{b}"] = val.str(30)
    if len(out) % 20 == 0:
        json.dump(out, open(f"shard_{wid}.json","w")); print(wid, len(out), len(mine), round(time.time()-t0,1), flush=True)
json.dump(out, open(f"shard_{wid}.json","w")); print("done", wid, round(time.time()-t0,1), flush=True)
