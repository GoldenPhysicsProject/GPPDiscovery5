#!/bin/bash
cd "$(dirname "$0")"; PY="nice taskset -c 7 /workspace/venv/bin/python"
EDGEK=0 $PY ../s57_ko.py bis odd 7 104 0.0009 0.0013 6 >> bis_odd.jsonl 2>>ko.err
EDGEK=0 $PY ../s57_ko.py lam odd 104 "log(7)/2" >> lam_odd.jsonl 2>>ko.err
echo done > D.done
