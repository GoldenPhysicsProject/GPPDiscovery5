#!/bin/bash
cd "$(dirname "$0")"; PY="nice taskset -c 7 /workspace/venv/bin/python $PWD/s58_ko.py"; export EDGEK=0
for M in 40 56 80; do $PY bis even 7 $M 0 0.015 10 >> bis_even.jsonl 2>>ko.err; done
$PY lam even 80 "log(2)/2,log(3)/2,log(5)/2,log(7)/2" >> lam_even.jsonl 2>>ko.err
echo done > D.done
