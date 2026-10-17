#!/bin/bash
cd "$(dirname "$0")"; PY="nice taskset -c 7 /workspace/venv/bin/python ../s57_ko.py"
export EDGEK=0
for M in 40 56 80; do $PY bis odd 3 $M 0 0.015 9 >> bis_odd.jsonl 2>>ko.err; done
for M in 40 56; do $PY bis odd 2 $M 0 0.04 9 >> bis_odd.jsonl 2>>ko.err; done
for M in 40 56; do $PY lam odd $M "log(2)/2,log(3)/2,log(5)/2,log(7)/2" >> lam_odd.jsonl 2>>ko.err; done
echo done > B.done
