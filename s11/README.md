# S11 — Gibbs-state / Bochner obstruction into the strip (PROVED LEMMA, classical)
Lemma S11.1. For sigma>1, phi_sigma(t)=zeta(sigma+it)/zeta(sigma) is the characteristic function of the prime-gas Gibbs law on -log n (PD, infinitely divisible, Levy measure sum Lambda(n)/log n n^-sigma delta_{-log n}); each Euler factor is a PD (geometric) factor iff the prime mode is unitary (S9/S10).
For every real sigma<=1, sigma!=1 pole aside, NO constant c makes t->c*zeta(sigma+it) positive definite (and at sigma=1 the boundary is singular).
Proof: continuous PD f has |f(t)|<=f(0). zeta(sigma+it) is unbounded in t for every sigma<=1 (Bohr-Landau / Titchmarsh ch. 8 for 1/2<=sigma<=1; functional equation growth for sigma<1/2), while c*zeta(sigma) is finite. QED.
Witnesses (float, mpmath 30 dps): output.txt (sigma=0.25,0.5,0.75 beat |zeta(sigma)| by t<=200; sigma near 1 needs larger t, theorem covers it).
Consequence (obstruction): the Gibbs state cannot be continued into the strip as a state (probability measure on spectrum). Any strip positivity must live on the test-function side (Weil/reflection form), not as Bochner positivity in t. Prime-gas "probabilities sum to 1" ends exactly at sigma=1 (= Hagedorn gap, #347).
Does NOT claim: anything about zeros or RH; novelty (classical: Khinchine; Lin-Hu 2001 zeta distribution).
