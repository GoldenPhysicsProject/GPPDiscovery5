#!/bin/bash
# Grok S37 H37a: Weil cert at L=27/20 (support 2.7). Assembly+tails via s36 schurboost (ball Schur disabled: --schur_bs ''),
# PD step afterwards via s36_cholres.py (1st impl) and s37_ldlt_exact.py (2nd impl) on the dumps.
cd "$(dirname "$0")"
T=${T:-700}; K=${K:-780}; NG=${NG:-92}
/workspace/venv/bin/python ../../s36/s36_onestroke_schurboost.py --L 27/20 --T $T --K $K --ng $NG --w 1/2 --b 2/5 --prec 384 \
  --workers ${W:-5} --torus ../torus_L1.35_D8x2_p235_p71113_cert.json --schur_bs '' \
  --dump_atil atil_L135_T${T}_K${K} --out cert_L135_T${T}_K${K}_ng${NG}_p384.json
