#!/bin/bash
# Grok S56 H56a: L = 69/50 (support 2.76) one parity, buffered lift (closest lift char log16: eps=0.01259).
cd "$(dirname "$0")"
T=${T:-1020}; K=${K:-1148}; NG=${NG:-92}
/workspace/venv/bin/python ../../s53/s53_onestroke_par.py --parities ${PARS:-0} --L 69/50 --T $T --K $K --ng $NG --w 1/2 --b 2/5 --prec 384 \
  --workers ${W:-5} --torus ${TORUS:-torus_L1.38_buf_PROVISIONAL.json} --schur_bs '' \
  --dump_atil atil_L138_T${T}_K${K} --out cert_L138_T${T}_K${K}_ng${NG}_p384_par${PARS:-0}.json
