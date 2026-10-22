#!/bin/bash
# Just above the cutoff: first synthetic g-prime 24.7 (log/2 = 1.60335) barely inside the window. X=23, half-width L.
P=/workspace/goddess/venv/bin/python; S=../../cosritz_px.py; C=../../ldl_certify.py; B=../../ldl_scan.py
gp="2,3,5,7,11,13,17,19,23,24.7"
for L in 16035/10000 16043/10000 16100/10000; do
  tag=edge; GPRIMES="$gp" TAG=_$tag $P $S $L 160 640 200 even >/dev/null
  f=Wdump_ZETA_${tag}_even_L${L/\//_}_K160_p640.json
  echo "{\"tag\":\"edge\",\"gprimes\":\"$gp\",\"L\":\"$L\"}"; $P $C $f 160 576; $P $B $f 160 576 ball | tail -1
  GPRIMES="2,3,5,7,11,13,17,19,23" TAG=_edgeref $P $S $L 160 640 200 even >/dev/null
  echo "{\"tag\":\"edgeref(zeta)\",\"L\":\"$L\"}"; $P $C Wdump_ZETA_edgeref_even_L${L/\//_}_K160_p640.json 160 576
done
