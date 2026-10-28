#!/bin/bash
# Chain 3 (S103): re-certify the cached T matrices with the rigorous budget (eps_J + eps_inv) folded into the shift,
# exact decimal ell, and both directions (PD claim at ell+2*extra, fail claim at ell-2*extra).  T.pkl from chain 2 (post-b9a215c).
PY=/workspace/sm/venv/bin/python; R=/workspace/sm/GPPDiscovery5/team/spectrum-measure/replay/mori_replay.py
cd /workspace/sm/runs
A12="--L log12/2 --Omega 416 --s 57/16 --N 832 --u0 207 --u1 68 --delta_log2 1515 --parity 0 --ells 1.2092639998e-54,1.2153e-54"
A0="--L 71/50 --Omega 1000 --s 763/200 --N 2376 --u0 607 --u1 200 --delta_log2 4200 --parity 0 --ells 1.0711e-81,1.07645555e-81"
A1="--L 71/50 --Omega 1000 --s 763/200 --N 2376 --u0 607 --u1 200 --delta_log2 4200 --parity 1 --ells 0,-4.0e-78,-4.2e-78"
A2="--L 71/50 --Omega 1150 --s 763/200 --N 2732 --u0 698 --u1 230 --delta_log2 4900 --parity 1 --ells 2.4134e-77,2.4254670e-77"
st() { tag=$1; shift; echo "$tag certify3 start $(date +%H:%M:%S)"; s=$(date +%s); $PY $R certify "$@" --out $tag > $tag.cert3.log 2>&1; echo "$tag certify3 exit $? $(date +%H:%M:%S) wall $(( $(date +%s)-s ))s"; }
st log12_p0 $A12
st L7150_O1000_p0 $A0 &
st L7150_O1000_p1 $A1 &
st L7150_O1150_p1 $A2 &
wait
echo CHAIN3_DONE $(date +%H:%M:%S)
