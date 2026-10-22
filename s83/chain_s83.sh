#!/bin/bash
# Grok S83 H80c step 3 (2026-10-22): rebuild the even support-2.84 dump (S71 ref sha b609d6f5...) with batched Arb cells
# (--chunk 1024: K = 1024 nodes per cell job, 36 jobs instead of 456; unchanged s71/mori_cert.py), compare within radius,
# then PSD: s36 cholres 0.999 pass / 1.005 control fail, s37 exact LDLt 0.999 pass / 1.005 control fail. RSS guard: kill if MemAvailable < 1.2 GB.
cd /workspace/GPPDiscovery5/s71; PY=/workspace/.venv-os/bin/python; R=../s83/run; tag=L7150_m12_O1000_p0_b1024
echo "start $(date +%H:%M:%S) $(date -u +%H:%M:%SZ)"
( while sleep 10; do a=$(awk '/MemAvailable/{print $2}' /proc/meminfo); if [ "$a" -lt 1200000 ]; then echo "RSS GUARD $(date +%H:%M:%S) avail=$a kB"; pkill -f "mori_cert.py --L 71/50"; fi; done ) & GUARD=$!
$PY mori_cert.py --L 71/50 --N 2376 --Omega 1000 --s 763/200 --m 1/2 --u0 607 --u1 200 --prec 512 --parity 0 --workers 4 --chunk 1024 --png 96 --delta_log2 4200 --out $R/cert_s83.jsonl --dump $R/T_$tag.json > $R/cert_$tag.log 2>&1 || { echo "ASSEMBLY FAILED"; tail -5 $R/cert_$tag.log; kill $GUARD; exit 1; }
kill $GUARD; f=$R/T_$tag.json; sha=$(sha256sum $f | cut -d' ' -f1); echo "assembled sha256 $sha $(date +%H:%M:%S)"
rn=$(tail -1 $R/cert_s83.jsonl | $PY -c "import json,sys;print(json.loads(sys.stdin.read())['radT_rowsum'])")
$PY ../s83/cmp_dumps.py $f run/T_L7150_m12_O1000_p0.json $rn 3.826567968227845e-145 | tee $R/cmp_$tag.json
$PY lam_from_dump.py $f > $R/lam_$tag.json; est=$($PY -c "import json;print(json.load(open('$R/lam_$tag.json'))['lam_min_est'])"); echo "lam_est $est"
$PY s36_cholres.py $f 640 0.999,1.005 $est > $R/cholres_$tag.log 2>&1; grep "fac=" $R/cholres_$tag.log
pos=$($PY -c "print('%.4e'%(0.999*$est))"); ctl=$($PY -c "print('%.4e'%(1.005*$est))")
$PY s37_ldlt_exact.py $f $sha 640 $pos,$ctl $R/ldlt_$tag.json > $R/ldlt_$tag.log 2>&1; tail -3 $R/ldlt_$tag.log
echo "CHAIN_S83_DONE $(date +%H:%M:%S)"
