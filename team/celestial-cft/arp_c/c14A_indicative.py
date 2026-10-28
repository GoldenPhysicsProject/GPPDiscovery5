"""C.14.A indicative check (NUMERICAL, mpmath): LHS(t) = sigma_S(t) + 2eps'(1+)(1-kappa_odd)/(1/4+t^2) vs c2,
with 1-kappa_odd = 0.58 (Golden Eye sinc-kernel plausibility value) and 2eps'(1+) = 45.993 (CC NUMERICAL).
sigma_S unknown (GAP); lower bounds used: sigma_S >= 0 (Bochner, T3) and sigma_S >= 2theta' - B/(1/4+t^2),
B = 2eps'(1+) + ||q||_1 (S&M SECTOR_COERCIVITY_p2 sec.5: B = 46 if ||q||_1 ~ 0, B = 66). Also the q-hat = 0 model.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import mpmath as mp
c2=mp.log(2)/mp.sqrt(2); E=mp.mpf('45.993'); g=E*mp.mpf('0.58')
O=lambda t: mp.re(mp.digamma(0.25+0.5j*t))-mp.log(mp.pi)
ts=[mp.mpf(i)/100 for i in range(0,4001)]
def scan(f,name):
    v=[(f(t)-c2,t) for t in ts]; m=min(v)
    bad=[t for d,t in v if d<=0]
    print(f"{name}: min(LHS-c2) = {mp.nstr(m[0],5)} at t = {mp.nstr(m[1],5)}; fails on",
          (f"[{mp.nstr(bad[0],4)}, {mp.nstr(bad[-1],4)}] ({len(bad)} grid pts)" if bad else "nothing (t in [0,40])"))
scan(lambda t: g/(0.25+t*t), "sigma_S>=0 only")
for B in (46,66):
    scan(lambda t,B=B: g/(0.25+t*t)+max(0,O(t)-B/(0.25+t*t)), f"sigma_S>=max(0,2theta'-B/(1/4+t^2)), B={B}")
scan(lambda t: O(t)-(E-g)/(0.25+t*t), "q-hat = 0 model (sigma_S = 2theta' - 2eps'/(1/4+t^2))")
