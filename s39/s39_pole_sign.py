#!/usr/bin/env python3
"""H39a.1: Independent audit of Codex polar sign.
Pole(f) = 2 Re(ell_+ conj(ell_-)) = ∫∫ f(x)conj(f(y)) 2 cosh((x-y)/2) dx dy
For real even: +2|ell_+|^2; real odd: -2|ell_+|^2.
Reproduce: odd f(x)=x on [-1,1] Pole≈-0.233574; even ≈8.68929 (identify test).
Label: CERTIFIED (mpmath dps>=50 analytic) / NUMERICAL.
"""
from mpmath import mp, mpf, quad, cosh, exp, fabs, nstr
mp.dps = 60

def ell_pm(f, L, sign):
    # ell_± = ∫_{-L}^L f(x) e^{±x/2} dx
    return quad(lambda x: f(x) * exp(sign * x / 2), [-L, L])

def pole_via_ells(f, L):
    ep = ell_pm(f, L, +1)
    em = ell_pm(f, L, -1)
    # real-valued f => 2 Re(ep * conj(em)) = 2*ep*em if real
    return 2 * (ep * em).real, ep, em

def pole_via_kernel(f, L):
    # double integral ∫∫ f(x)f(y) 2 cosh((x-y)/2)
    return quad(lambda x: f(x) * quad(lambda y: f(y) * 2 * cosh((x - y) / 2), [-L, L]), [-L, L])

results = {}

# --- Odd test: f(x)=x on [-1,1] ---
L = mpf(1)
f_odd = lambda x: x
P_ells, ep, em = pole_via_ells(f_odd, L)
P_ker = pole_via_kernel(f_odd, L)
# Analytic: ell_+ = ∫_{-1}^1 x e^{x/2} dx, ell_- = ∫ x e^{-x/2}
# Integration by parts: ∫ x e^{ax} = e^{ax}(ax-1)/a^2
# For odd f, ell_+ = - ell_- (real), Pole = -2 |ell_+|^2
two_abs = -2 * (ep ** 2)  # since ep real and em=-ep for odd real
results['odd_x'] = {
    'Pole_ells': str(P_ells),
    'Pole_kernel': str(P_ker),
    'ell_plus': str(ep),
    'ell_minus': str(em),
    'minus_2_abs_ell_plus_sq': str(two_abs),
    'Codex_claim': '-0.233574',
    'match_ells_vs_minus2': float(fabs(P_ells - two_abs)) < 1e-40,
    'match_ells_vs_kernel': float(fabs(P_ells - P_ker)) < 1e-25,
    'approx_float': float(P_ells),
}

# Analytic closed form for odd f(x)=x on [-1,1]:
# ell_+ = ∫_{-1}^1 x e^{x/2} dx = [e^{x/2}(x/2 - 1)/(1/4)*? wait a=1/2:
# ∫ x e^{ax} dx = e^{ax}(a x - 1)/a^2
a = mpf('0.5')
def I_xa(a, L):
    # ∫_{-L}^L x e^{a x} dx = [e^{a x}(a x - 1)/a^2 ]_{-L}^L
    def F(x):
        return exp(a * x) * (a * x - 1) / (a * a)
    return F(L) - F(-L)
ep_an = I_xa(a, L)
em_an = I_xa(-a, L)
Pole_an = 2 * ep_an * em_an
results['odd_x_analytic'] = {
    'ell_plus': str(ep_an),
    'ell_minus': str(em_an),
    'Pole': str(Pole_an),
    'Pole_float': float(Pole_an),
    'minus_2_ep_sq': float(-2 * ep_an ** 2),
}

# --- Even tests to identify Codex's 8.68929 ---
# Try f=1 on [-1,1], f=1-x^2, f=cos, etc.
candidates = []
def try_even(name, f, L):
    P, ep, em = pole_via_ells(f, L)
    two = 2 * (ep ** 2)
    candidates.append({
        'name': name, 'L': float(L),
        'Pole': float(P), 'plus_2_abs': float(two),
        'ell_plus': float(ep),
        'diff_to_8.68929': abs(float(P) - 8.68929),
        'parity_ok': abs(float(P - two)) < 1e-20,
    })

try_even('f=1 on [-1,1]', lambda x: mpf(1), mpf(1))
try_even('f=1-x^2 on [-1,1]', lambda x: 1 - x*x, mpf(1))
try_even('f=cos(x) on [-1,1]', lambda x: mp.cos(x), mpf(1))
try_even('f=1 on [-2,2]', lambda x: mpf(1), mpf(2))
try_even('f=1-x^2/4 on [-2,2]', lambda x: 1 - (x*x)/4, mpf(2))
try_even('f=exp(-x^2) on [-1,1]', lambda x: exp(-x*x), mpf(1))
try_even('f=sech(x) on [-1,1]', lambda x: 1/mp.cosh(x), mpf(1))
# maybe f=cosh(x/2) related
try_even('f=cosh(x/2) on [-1,1]', lambda x: cosh(x/2), mpf(1))
try_even('f=e^{-|x|/2} on [-1,1] approx smooth', lambda x: exp(-fabs(x)/2), mpf(1))

# Analytic for f=1 on [-L,L]: ell_+ = ∫ e^{x/2} = 2(e^{L/2}-e^{-L/2}) = 4 sinh(L/2)
# Pole = 2 (2 sinh(L/2)*2 / (1/2)? wait:
# ∫_{-L}^L e^{a x} dx = (e^{aL}-e^{-aL})/a = 2 sinh(aL)/a
# a=1/2: ell_+ = 2 sinh(L/2)/(1/2) = 4 sinh(L/2)
# Pole = 2*(4 sinh(L/2))^2 = 32 sinh^2(L/2)
L1 = mpf(1)
Pole_const = 32 * mp.sinh(L1/2)**2
candidates.append({
    'name': 'analytic f=1 [-1,1]',
    'Pole': float(Pole_const),
    'diff_to_8.68929': abs(float(Pole_const) - 8.68929),
})

# Sort by closeness to 8.68929
candidates.sort(key=lambda c: c.get('diff_to_8.68929', 1e9))
results['even_candidates'] = candidates[:8]
results['best_even_match'] = candidates[0]

# Also check: maybe they used ||f|| normalization or L=pi?
for Ltry in [mpf('1.5'), mpf('2'), mp.pi/2, mp.pi]:
    P, ep, em = pole_via_ells(lambda x: mpf(1), Ltry)
    candidates.append({'name': f'f=1 L={float(Ltry)}', 'Pole': float(P),
                       'diff_to_8.68929': abs(float(P)-8.68929)})

# f(x)=1-|x| on [-1,1]
P, ep, em = pole_via_ells(lambda x: 1 - fabs(x), mpf(1))
results['even_tent'] = {'Pole': float(P), 'diff': abs(float(P)-8.68929), 'plus2': float(2*ep**2)}

# Print summary
import json
print(json.dumps({
    'odd_x_float': results['odd_x']['approx_float'],
    'odd_analytic_float': results['odd_x_analytic']['Pole_float'],
    'odd_minus2': results['odd_x_analytic']['minus_2_ep_sq'],
    'odd_match_codex': abs(results['odd_x']['approx_float'] + 0.233574) < 1e-5,
    'parity_odd_ok': results['odd_x']['match_ells_vs_minus2'],
    'best_even': results['best_even_match'],
    'even_top3': candidates[:3],
}, indent=2, default=str))

# Save full
with open('/workspace/GPPDiscovery5/s39/h39a_pole_sign.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print('wrote h39a_pole_sign.json')
