#!/bin/bash
# Grok S54 H54a': full-complex cert at L = 677/500 = 1.354 (support 2.708), just below log(15)/2 = 1.3540251.
# The S37 certified D8x2 lift (A_cert <= 4.5928) stays admissible up to L < log(15)/2: its closest character
# (0,1,1) has frequency log 15 = 2.7080502 > 2L = 2.708 (Arb-checked by the assembly's visibility assert).
# The comb is unchanged (n <= 13). Parities run sequentially (memory). Usage: PARS=0 ./run.sh ; PARS=1 ./run.sh
cd "$(dirname "$0")"
T=${T:-800}; K=${K:-884}; NG=${NG:-92}
/workspace/venv/bin/python ../../s53/s53_onestroke_par.py --parities ${PARS:-0} --L 677/500 --T $T --K $K --ng $NG --w 1/2 --b 2/5 --prec 384 \
  --workers ${W:-5} --torus ../../s37/torus_L1.35_D8x2_p235_p71113_cert.json --schur_bs '' \
  --dump_atil atil_L1354_T${T}_K${K} --out cert_L1354_T${T}_K${K}_ng${NG}_p384_par${PARS:-0}.json
