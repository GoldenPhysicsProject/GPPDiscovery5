#!/bin/bash
cd "$(dirname "$0")"
for M in 40 56 72 96; do /workspace/venv/bin/python ../s56_ko_bisect.py odd 3 $M 0.03 6 >> ko_odd_p3.jsonl 2>>ko.err; done
for M in 40 72; do /workspace/venv/bin/python ../s56_ko_bisect.py even 5 $M 0.03 6 >> ko_even_p5.jsonl 2>>ko.err; done
echo done > B.done
