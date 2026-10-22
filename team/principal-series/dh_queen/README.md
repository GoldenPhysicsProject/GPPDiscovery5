# team/principal-series/dh_queen (New Bot, session 2; RH OPEN)
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
L = HALF-width throughout (window [-L,L]); s75/dhritz.py and team/queen dh_sinh.py both use it.
- overlap.py: node counts of DH lowest even Ritz vectors (s75 cosvec witnesses + moderate K) and overlap with Queen's f_Q = e^{-i g0 x} sinh(b0 x). -> overlap_s75.jsonl, overlap_moderateK.jsonl
- ovcheck.py: dps-40 closed-form recheck of the overlap.
- parity_split.py: Q_DH(Re f_Q) (odd) and Q_DH(Im f_Q) (even) via Queen's Levy-Khintchine formula; needs a copy of team/queen/dh_sinh.py as dh_sinh_queen_copy.py (git show origin/team/queen:team/queen/dh_sinh.py). -> parity_split.json; qcross.py -> qcross.log
- run_moderateK.sh: unmodified s75/dhritz.py at L=3421/2000, K=60,120,180 (python-flint venv). -> dh_moderateK.jsonl, cosvec_*.json
All NUMERICAL except the S75 K300/K360 negativity at 1.7105 (CERTIFIED, S75).
Note: bridge work/grok/team/principal-series.md, session 2.
