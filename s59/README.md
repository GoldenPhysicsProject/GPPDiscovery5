# Grok S59 (2026-10-17). RH OPEN. GRH OPEN. φ absent.
See GPP-bridge work/grok/2026-10-17_s59_crossgroup_lift_knockout.md for results and labels.

- `lift/`: H59a test of the S58 H58c conjecture (cross-group lift). `s59_joint_lp.py DC budget tag [init]` = S55 proximal cutting-plane
  LP on the 6-torus {2,3,5,7,11,13}, L = 7/5, buffer EPS = 0.0056, started from the S58 CERTIFIED separable lift (A <= 5.189) plus all
  cross characters with |k|_1 <= DC. `s59_supcheck.py` = independent heavier 6-D sup estimate (products of group maxima on 80^3 grids +
  400k random, 8000 refined). Validation: it reproduces the separable sup 5.18343 (= sum of the two S58 B&B best samples).
  Results (NUMERICAL, sup estimates, not B&B): separable 5.1834; DC = 4 (287 cross chars) 5.1753; DC = 5 (993 cross chars) 5.1684.
  Logs dc4.log, dc5.log (dc5 stopped by hand after it 8), sup_*.jsonl. Lift jsons joint_dc4.json, joint_dc5.json committed.
- `ko/`: H59b knockout margins (`s59_ko.py` = s58_ko.py), EDGEK = 0. `chain.sh` runs detached; results in bis.jsonl / lam.jsonl,
  done-file ko.done.
