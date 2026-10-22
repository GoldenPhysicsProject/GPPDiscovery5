# team/principal-series/epstein (New Bot S2; RH OPEN)
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. L = HALF-width (window [-L,L]).
- lkritz.py: even-sector cos-basis Galerkin (same basis as s75/dhritz.py) of the Queen/Celestial Levy-Khintchine window form,
  for EPS (x^2+5y^2: Celestial's Omega_E(0), nu_E = 1/sinh(u/2), her rQ/logderiv recursion re-run in mpmath), DH (Queen's data), ZETA.
  Needs copies of team/celestial-cft/queen_probe/epstein_sinh_lib.py (as epstein_sinh_lib_celestial_copy.py) and team/queen/dh_sinh.py.
  Validation: DH at L=1.7105 K=60 gives ritz3 = 6.35207e-16 vs dhritz Arb 6.35218e-16; the two bottom (~1e-29, 6e-22) eigenvalues are below
  this code's quadrature floor (~1e-21) and are NOT resolved here. Epstein eigenvalues quoted are O(1e-4) or larger.
- eps_overlap.py: nodes of the lowest even vector and overlaps with Im f_Q at Celestial's 13 CERTIFIED x^2+5y^2 zeros (certified_zero_boxes.json copy).
NUMERICAL throughout.
