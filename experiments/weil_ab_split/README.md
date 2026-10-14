# weil_ab_split

NUMERICAL (float). RH OPEN.

Split the Weil quadratic form of `weil_ritz.py` into
- **A**: archimedean/gamma digamma integral
- **B**: pole terms − von Mangoldt prime terms

and report signed Rayleigh quotients of A and B on the Ritz minimizer of A+B,
for ζ and F_θ (esp. θ=1/4) on the known failure window h∈[0.40,0.52].

## Run

```bash
source /workspace/venv/bin/activate
python weil_ab_split.py 18 0.40 0.52 0.01
```

## Outputs

- `result_M18_0.400_0.520.json` — own-minimizer A/B and cross-eval
- `result_cross_subsplit_M18.json` — cross-eval with B split into pole vs primes
- `out_M18_0.40_0.52.txt` — console log

See GPP-bridge `work/grok/2026-10-14_s7_weil_ab_split.md`.
