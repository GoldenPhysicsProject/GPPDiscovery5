#!/bin/bash
# Pinning test: move the rational prime 29 by +-10^-e inside the L=1.7105 window (half-width), X=23. 
P=/workspace/goddess/venv/bin/python; S=../../cosritz_px.py; C=../../ldl_certify.py; B=../../ldl_scan.py
base="2,3,5,7,11,13,17,19,23"
for e in 30 40 50; do for sg in p m; do
  z=$(printf '%0*d' $((e-1)) 0); if [ $sg = p ]; then q="29.${z}1"; else q="28.$(printf '9%.0s' $(seq 1 $e))"; fi
  tag=d${sg}${e}; GPRIMES="$base,$q" TAG=_$tag $P $S 3421/2000 200 640 200 even >/dev/null
  f=Wdump_ZETA_${tag}_even_L3421_2000_K200_p640.json
  echo "{\"tag\":\"$tag\",\"delta\":\"${sg}1e-$e\"}"; $P $C $f 200 576; $P $B $f 200 576 ball | tail -1
done; done
