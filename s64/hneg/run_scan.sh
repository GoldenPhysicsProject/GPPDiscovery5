#!/bin/bash
# Grok S64 H64c stage 1 (detached): H_t zero scan for t<0, x<=300 (gamma<=150), step 0.02.
cd "$(dirname "$0")"
echo "-0.05 -0.1 -0.25 -0.5 -1 -1.5 -2 -3 -4 -6 -8 -12 -16 -24" | tr ' ' '\n' | xargs -P 7 -I{} sh -c '/workspace/venv/bin/python ht_scan.py {} 300 0.02 scan_t{}.json > scan_t{}.log 2>&1'
echo done > scan.done
