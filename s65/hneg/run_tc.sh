#!/bin/bash
cd "$(dirname "$0")"
for t in -0.25 -0.3 -0.35 -0.4 -0.45 -0.5; do echo $t; done | xargs -P 3 -I{} sh -c '/workspace/venv/bin/python tc_track.py {} >> tc.jsonl 2>>tc.err'
echo done > tc.done
