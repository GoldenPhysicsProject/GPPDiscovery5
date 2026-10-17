#!/bin/bash
cd "$(dirname "$0")"; PY="nice taskset -c 4 /workspace/venv/bin/python $PWD/s58_ko.py"; export EDGEK=0
for M in 40 56; do $PY bis even 2 $M 0 0.12 10 >> bis_even.jsonl 2>>ko.err; done
$PY lam even 56 "log(2)/2,log(3)/2,log(5)/2,log(7)/2" >> lam_even.jsonl 2>>ko.err
$PY bis even 2 80 0 0.12 10 >> bis_even.jsonl 2>>ko.err
echo done > A.done
