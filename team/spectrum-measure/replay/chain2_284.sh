#!/bin/bash
# Chain 2 (after the Arb bessel_j turning-point fix in sph_j_rows). Resumable: cells/panels cached on disk.
PY=/workspace/sm/venv/bin/python; R=/workspace/sm/GPPDiscovery5/team/spectrum-measure/replay/mori_replay.py
cd /workspace/sm/runs
A0="--L 71/50 --Omega 1000 --s 763/200 --N 2376 --u0 607 --u1 200 --delta_log2 4200 --parity 0 --chunk 48 --ells 1.0711e-81,1.07645555e-81"
A1="--L 71/50 --Omega 1000 --s 763/200 --N 2376 --u0 607 --u1 200 --delta_log2 4200 --parity 1 --chunk 48 --ells 0,-4.0e-78,-4.2e-78"
A2="--L 71/50 --Omega 1150 --s 763/200 --N 2732 --u0 698 --u1 230 --delta_log2 4900 --parity 1 --chunk 40 --ells 2.4134e-77,2.4254670e-77"
st() { tag=$1; stage=$2; shift 2; echo "$tag $stage start $(date +%H:%M:%S)"; $PY $R $stage "$@" --out $tag >> $tag.log 2>&1; echo "$tag $stage exit $? $(date +%H:%M:%S)"; }
st L7150_O1000_p0 panels $A0 --workers 2
# wait for old chain's odd Omega=1000 cells to finish, then stop it (its panels use the pre-fix code)
while [ $(ls L7150_O1000_p1/cells/*.pkl 2>/dev/null | wc -l) -lt 12 ]; do sleep 30; done
sleep 20
pkill -f chain_284.sh; pkill -f "mori_replay.py all --L 71/50 --Omega 1000"; sleep 5
rm -rf L7150_O1000_p1/panels; mkdir -p L7150_O1000_p1/panels; echo "old chain stopped $(date +%H:%M:%S)"
st L7150_O1000_p0 assemble $A0 --workers 4
st L7150_O1000_p0 certify $A0 --workers 4
st L7150_O1000_p1 panels $A1 --workers 4
st L7150_O1000_p1 assemble $A1 --workers 4
st L7150_O1000_p1 certify $A1 --workers 4
st L7150_O1150_p1 all $A2 --workers 4
echo CHAIN2_DONE $(date +%H:%M:%S)
