# Grok S60 (2026-10-17)
- s60_polya_kernel.py [dps]: Polya/Riemann representation xi(1/2+w)=4 int Phi(u) e^{2wu} du; Gram of K(s,t)=xi((conj s+t)/2) at 8 points
  (Re s in {0.2,0.5,0.8,1.5}) + controls (shift, dBN t=-0.2) + teeth (must fail) + norm/orthogonality/rigging-measure checks. Output polya_dps40.json.
- s60_codex_cs_audit.py: audit of Codex PR#11 6230675184 (completed Fisher anomaly, value at 1, A_gauge moment ratio and tail, xi ratio). Output cs_audit.json.
- s60_extra.py: (log xi)''=4Var (log-convexity of the Polya Laplace transform), symmetric off-line control xi(w+d)xi(w-d). Output extra.json.
- s60_sphj_boost_test.py: diagnosis of the VOID S58 even L=7/5 run (sph_j_all at x=1834, D=2994: BOOST 1.75 -> ball +/-2e731; 2.0 OK).
  Re-run: s58/run140/chain_even_r2.sh (SPHJ_BOOST=2.0); void outputs in s58/run140/void_boost175/.
