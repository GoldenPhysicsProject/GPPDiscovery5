#!/bin/bash
# Grok S53 H53a: ODD sector only at L=27/20 (support 2.7), same certified D8x2 lift as S37, larger T (beta 0.119 -> 0.253).
cd "$(dirname "$0")"
T=${T:-800}; K=${K:-880}; NG=${NG:-92}
/workspace/venv/bin/python ../s53_onestroke_par.py --parities ${PARS:-1} --L 27/20 --T $T --K $K --ng $NG --w 1/2 --b 2/5 --prec 384 \
  --workers ${W:-5} --torus ../../s37/torus_L1.35_D8x2_p235_p71113_cert.json --schur_bs '' \
  --dump_atil atil_L135_T${T}_K${K} --out cert_L135_T${T}_K${K}_ng${NG}_p384_par${PARS:-1}.json
