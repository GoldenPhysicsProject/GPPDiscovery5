"""Grok S61: markdown table from budget/budget_M56.jsonl (signed decomposition of the Weil minimizer, NUMERICAL)."""
import json
rows = [json.loads(l) for l in open('budget/budget_M56.jsonl')]
print('| L | par | λ_min | pole | arch ω0N+E_γ | primes (signed) | E_prime/(2S_L N) | E_γ share of fwd | fwd = rev |')
print('|---|---|---|---|---|---|---|---|---|')
for d in rows:
    if 'L' not in d: continue
    for p in ('even', 'odd'):
        e = d[p]
        print(f"| {d['L']} | {p} | {float(e['lam_min']):.3e} | {float(e['pole']):+.4f} | {float(e['arch_total']):+.4f} | {float(e['prime_total_signed']):+.4f} | {float(e['Eprime_over_2S']):.4f} | {float(e['share_forward']['E_gamma']):.3f} | {float(e['forward']):.4f} |")
