#!/bin/bash
cd "$(dirname "$0")"; PY="nice taskset -c 6 /workspace/venv/bin/python ../s57_ko.py"
export EDGEK=0
for M in 40 56 80; do $PY bis odd 5 $M 0 0.006 9 >> bis_odd.jsonl 2>>ko.err; done
for M in 40 56; do $PY bis odd 7 $M 0 0.012 9 >> bis_odd.jsonl 2>>ko.err; done
$PY lam odd 80 "log(2)/2,log(3)/2,log(5)/2,log(7)/2" >> lam_odd.jsonl 2>>ko.err
echo done > A.done
