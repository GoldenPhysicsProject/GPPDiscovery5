#!/bin/bash
# Grok S55 H55a frontier: full-complex cert at L = 2701/2000 = 1.3505 (support 2.701) with the CERTIFIED S37 D8x2 lift
# (A <= 4.5928). Closest lift character (0,1,1) at log 15 = 2.7080502: gap eps = 7.05e-3, eps*T = 5.64 at T=800
# (S54 buffer rule eps*T >~ 5; S37 had 5.6 at T=700 and passed; S54 L=1.354 had 0.04 and failed). Tests the rule.
cd "$(dirname "$0")"
T=${T:-800}; K=${K:-881}; NG=${NG:-92}
/workspace/venv/bin/python ../../s53/s53_onestroke_par.py --parities ${PARS:-0} --L 2701/2000 --T $T --K $K --ng $NG --w 1/2 --b 2/5 --prec 384 \
  --workers ${W:-4} --torus ../../s37/torus_L1.35_D8x2_p235_p71113_cert.json --schur_bs '' \
  --dump_atil atil_L13505_T${T}_K${K} --out cert_L13505_T${T}_K${K}_ng${NG}_p384_par${PARS:-0}.json
