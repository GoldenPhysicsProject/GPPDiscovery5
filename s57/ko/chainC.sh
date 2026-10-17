#!/bin/bash
cd "$(dirname "$0")"; PY="nice taskset -c 6 /workspace/venv/bin/python"
for f in ../wit/w_p5_odd_k0_M56_L8073_10000.json ../wit/w_p5_odd_k0_M56_L101_125.json; do $PY ../s57_witness_check.py $f 70 >> ../wit/check.jsonl 2>>ko.err; done
EDGEK=0 $PY ../s57_ko.py bis odd 7 80 0.001 0.0016 7 >> bis_odd.jsonl 2>>ko.err
EDGEK=0 $PY ../s57_ko.py bis odd 11 40 0 0.008 8 >> bis_odd.jsonl 2>>ko.err
EDGEK=0 $PY ../s57_ko.py lam odd 40 "log(11)/2" >> lam_odd.jsonl 2>>ko.err
echo done > C.done
