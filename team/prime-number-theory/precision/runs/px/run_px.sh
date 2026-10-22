#!/bin/bash
# P_X positivity tests. Half-width L. X=23 so (log X)/2 = 1.5678.
P=/workspace/goddess/venv/bin/python; S=../../cosritz_px.py; C=../../ldl_certify.py; B=../../ldl_scan.py
base="2,3,5,7,11,13,17,19,23"
run(){ tag=$1; gp=$2; L=$3; K=$4
  GPRIMES=$gp TAG=_$tag $P $S $L $K 640 200 even >/dev/null
  f=Wdump_ZETA_${tag}_even_L${L/\//_}_K${K}_p640.json
  echo "{\"tag\":\"$tag\",\"gprimes\":\"$gp\",\"L\":\"$L\"}"
  $P $C $f $K 576 ; $P $B $f $K 320 ball | tail -1; }
run zeta29 "$base,29" 3421/2000 200
run m29p "$base,29.001" 3421/2000 200
run m29m "$base,28.999" 3421/2000 200
run del29 "$base" 3421/2000 200
run syn2 "$base,24.7,28.3" 3421/2000 200
run syn2_L162 "$base,24.7,28.3" 81/50 160
run zeta_L162 "$base" 81/50 160
