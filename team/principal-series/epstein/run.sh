#!/bin/bash
cd /workspace/newbot/epstein
for K in 40 60 90; do for L in 1.5 2 2.5 3; do nice /workspace/sm/venv/bin/python lkritz.py EPS $L $K >> eps_ritz.jsonl 2>>eps_ritz.err; done; done
