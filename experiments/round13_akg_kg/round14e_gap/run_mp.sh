#!/bin/sh
# Round 14E2 driver: wait for each W matrix, then run mpgap for both whitening cuts.
cd "$(dirname "$0")"
for job in "ZETA 1.5 45 30" "ZETA 1.5 90 51" "ZETA 1.75 45 34" "ZETA 1.75 90 59" "ZETA 2.0 45 37" "ZETA 2.0 90 66" "FTH 2.0 45 37" "FTH 2.0 90 66"; do
  set -- $job
  for tau in 1e-13 1e-40; do
    echo "while [ ! -f mp80/W_$1_$2_$4.pkl ]; do sleep 10; done; python3 mpgap.py $1 $2 $3 $tau > mp80/gap_$1_$2_$3_$tau.log 2>&1"
  done
done | xargs -P 6 -I{} sh -c "{}"
