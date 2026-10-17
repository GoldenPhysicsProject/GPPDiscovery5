#!/bin/bash
cd "$(dirname "$0")"; PY="nice taskset -c 5 /workspace/venv/bin/python $PWD/s58_ko.py"; export EDGEK=0
for M in 40 56 80; do $PY bis even 3 $M 0 0.05 10 >> bis_even.jsonl 2>>ko.err; done
echo done > B.done
