#!/bin/bash
# Grok S55 H55b knockout refinement queue (M=56, dps130, 1024 bits). NUMERICAL.
cd "$(dirname "$0")"; PY=/workspace/venv/bin/python; export OMP_NUM_THREADS=1
while pgrep -f "s55_knockout.py 56 1.205" >/dev/null; do sleep 5; done
KO=17,19 $PY ../s55_knockout.py 56 1.42,1.425,1.43,1.435,1.44,1.45,1.475,1.48,1.485,1.49,1.5 > ko_p17_19.jsonl 2> ko_p17_19.err
KO=5 $PY ../s55_knockout.py 56 0.8075,0.81,0.8125,0.815,0.8175 > ko_p5.jsonl 2> ko_p5.err
echo done > chain.done
