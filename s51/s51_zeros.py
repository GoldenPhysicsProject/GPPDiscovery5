# Grok S51: Arb-certified zeta zeros (python-flint acb.zeta_zeros, Turing-method isolation) -> text file
import flint, sys, time
from flint import acb, ctx
ctx.prec = 1400
K = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
t = time.time(); out = []
for s in range(1, K + 1, 500):
    zz = acb.zeta_zeros(s, min(500, K - s + 1))
    for z in zz:
        assert z.real == 0.5  # exact 1/2 (on the line, certified)
        out.append(z.imag)
with open('zeros_arb_%d.txt' % K, 'w') as f:
    for g in out: f.write(g.mid().str(130, radius=False) + ' ' + g.rad().str(5, radius=False) + '\n')
print(K, 'zeros', time.time() - t, 's; last', out[-1].mid().str(20), 'max rad', max(float(g.rad()) for g in out))
