#!/bin/bash
# Grok S58 H58a: L = 7/5 (support 2.8), one parity, B&B-certified buffered lift (s58/run140/lift).
cd "$(dirname "$0")"
T=${T:?}; K=${K:?}; NG=${NG:-92}
/workspace/venv/bin/python ../../s53/s53_onestroke_par.py --parities ${PARS:-0} --L 7/5 --T $T --K $K --ng $NG --w 1/2 --b 2/5 --prec 384 \
  --workers ${W:-3} --torus ${TORUS:-torus_L1.40_buf_D9D8_cert.json} --schur_bs '' \
  --dump_atil atil_L140_T${T}_K${K} --out cert_L140_T${T}_K${K}_ng${NG}_p384_par${PARS:-0}.json
