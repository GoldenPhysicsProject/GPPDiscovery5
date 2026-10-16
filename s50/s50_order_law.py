# S50: lattice-model order law. eps_+(n) upper band edge (closed form); order n (union of orders 1..n) covers all
# mid-gap delta <= 1/2 up to height gamma_n = 2 pi exp(4 pi eps_+(n)) (Delta = 2pi/log(gamma/2pi)). HEURISTIC (lattice model).
import math, json
def ep(n):
    if n == 1: return math.sqrt(8)/(2*math.pi)
    M = 2*n; return math.sqrt(M*(M - 1) + math.sqrt(M*(M - 1)*(4*M - 6)))/(2*math.pi)
def em(n):
    if n == 1: return 0.0
    M = 2*n; return math.sqrt(M*(M - 1) - math.sqrt(M*(M - 1)*(4*M - 6)))/(2*math.pi)
rows = []
for n in range(1, 13):
    rows.append(dict(n=n, eps_minus=em(n), eps_plus=ep(n), asym_plus=(2*n + math.sqrt(2*n) - 0.5)/(2*math.pi),
                     asym_minus=(2*n - math.sqrt(2*n) - 0.5)/(2*math.pi), log10_gamma_n=math.log10(2*math.pi) + 4*math.pi*ep(n)/math.log(10)))
    print(rows[-1])
H = 3e12; eH = math.log(H/(2*math.pi))/(4*math.pi)
nH = min(n for n in range(1, 50) if ep(n) >= eH)
print('at H=3e12: delta=1/2 is eps =', eH, ' order needed (lattice) =', nH)
json.dump(dict(rows=rows, eps_half_at_H=eH, order_needed_at_H=nH), open('s50_order_law.json', 'w'), indent=1)
