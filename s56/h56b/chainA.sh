#!/bin/bash
cd "$(dirname "$0")"
for M in 40 56 72 96; do /workspace/venv/bin/python ../s56_ko_bisect.py odd 5 $M 0.03 6 >> ko_odd_p5.jsonl 2>>ko.err; done
for K in 2 6; do EDGEK=$K /workspace/venv/bin/python ../s56_ko_bisect.py odd 5 56 0.03 6 >> ko_odd_p5_edgek.jsonl 2>>ko.err; done
echo done > A.done
