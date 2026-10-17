#!/bin/bash
# Grok S59 H59b: knockout margins (EDGEK=0 free-edge Legendre Ritz), sequential, one core, niced.
cd "$(dirname "$0")"; PY="nice -n 15 taskset -c 7 /workspace/venv/bin/python $PWD/s59_ko.py"; export EDGEK=0
$PY bis odd 11 56 0 0.0025 10 >> bis.jsonl 2>>ko.err
$PY bis even 11 56 0 0.006 10 >> bis.jsonl 2>>ko.err
$PY bis even 7 104 0 0.0045 10 >> bis.jsonl 2>>ko.err
$PY bis odd 11 80 0 0.0025 10 >> bis.jsonl 2>>ko.err
$PY bis even 11 80 0 0.006 10 >> bis.jsonl 2>>ko.err
$PY lam odd 80 "log(11)/2" >> lam.jsonl 2>>ko.err
$PY lam even 80 "log(11)/2" >> lam.jsonl 2>>ko.err
echo done > ko.done
