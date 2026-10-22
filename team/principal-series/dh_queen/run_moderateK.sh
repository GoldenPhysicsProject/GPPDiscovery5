#!/bin/bash
cd /workspace/newbot/dh
for K in 60 120 180; do FORM=DH nice /workspace/sm/venv/bin/python dhritz.py 3421/2000 $K 320 60 even >> dh_moderateK.jsonl 2>>dh_moderateK.err; done
