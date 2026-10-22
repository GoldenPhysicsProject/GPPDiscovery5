#!/bin/bash
# Independent replay chain (resumable: finished cells/panels are cached on disk). Spectrum-measure seat, 2026-10-22.
PY=/workspace/sm/venv/bin/python; R=/workspace/sm/GPPDiscovery5/team/spectrum-measure/replay/mori_replay.py
cd /workspace/sm/runs
run() { tag=$1; shift; echo "$tag start $(date +%H:%M:%S)"; $PY $R all "$@" --out $tag > $tag.log 2>&1; echo "$tag exit $? $(date +%H:%M:%S)"; }
run L7150_O1000_p0 --L 71/50 --Omega 1000 --s 763/200 --N 2376 --u0 607 --u1 200 --delta_log2 4200 --parity 0 --workers 4 --chunk 48 --ells 1.0711e-81,1.07645555e-81
run L7150_O1000_p1 --L 71/50 --Omega 1000 --s 763/200 --N 2376 --u0 607 --u1 200 --delta_log2 4200 --parity 1 --workers 4 --chunk 48 --ells 0,-4.0e-78,-4.2e-78
run L7150_O1150_p1 --L 71/50 --Omega 1150 --s 763/200 --N 2732 --u0 698 --u1 230 --delta_log2 4900 --parity 1 --workers 4 --chunk 40 --ells 2.4134e-77,2.4254670e-77
echo CHAIN_DONE $(date +%H:%M:%S)
