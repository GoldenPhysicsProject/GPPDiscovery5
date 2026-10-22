# Grok S90 (2026-10-22): prime anatomy of the in-band cancellation (H89a). RH OPEN.
`prime_anatomy.py` (S87 builder, mpmath dps 90, N=24; sign conventions from s89/calib.py). For ζ at L = 0.6, 0.8, 1.0, 1.2:
(1) splits vᵀQv for the near-null v into arch in band, arch out of band, pole, and each prime power −Λ(n)vᵀP_nv;
(2) gives the unit-atom profile g(x) on (1, e^{2L});
(3) runs exact perturbations: remove or flip one prime (all of its powers), primes only, Λ → (1±1%)Λ, jitter of every prime power x → x(1+δ), 6/10/12 inserted as "primes", and the smooth PNT density dψ = dx.
Log: `anatomy.log`; JSON: `prime_anatomy.json` (written when the run ends). HIGH-PRECISION NUMERICAL, not interval.
Prior art: Connes–Consani, "Spectral triples and ζ-cycles" (arXiv 2106.01715) §2.3–2.5 and §3. Moving p=2 off 2 by more than about 1e−3 breaks positivity. Dropping a prime power breaks it. The extremely small eigenvalues are explained by the near-radical E(f)(x) = x^{1/2} Σ_n f(nx) applied to prolate functions.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
