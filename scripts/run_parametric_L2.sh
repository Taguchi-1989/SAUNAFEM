#!/bin/bash
# Run Case L-2: 18kW buoyantPimpleFoam transient 3600s
set -e
source /usr/lib/openfoam/openfoam2312/etc/bashrc 2>/dev/null || true
export FOAM_SIGFPE=false

BASE="/mnt/d/dev/SaunaFEM"
SRC="$BASE/results/L2"
CASE="$HOME/saunaflow_run/L2"

echo "============================================"
echo "=== L-2: 18kW transient (3600s) ==="
echo "============================================"

rm -rf "$CASE"
mkdir -p "$(dirname "$CASE")"
cp -r "$SRC" "$CASE"
cd "$CASE"

find . -type f ! -path '*/polyMesh/*' | xargs -r sed -i 's/\r$//'

echo "=== blockMesh ===" && blockMesh > log.blockMesh 2>&1 && echo "done"

if [ -f system/topoSetDict ]; then
    echo "=== topoSet ===" && topoSet > log.topoSet 2>&1 && echo "done"
fi

if grep -q "radiationModel.*viewFactor" constant/radiationProperties 2>/dev/null; then
    echo "=== faceAgglomerate ===" && faceAgglomerate -dict constant/viewFactorsDict > log.faceAgglomerate 2>&1 && echo "done"
    echo "=== viewFactorsGen ===" && viewFactorsGen > log.viewFactorsGen 2>&1 && echo "done"
fi

echo "=== buoyantPimpleFoam (3600s) ==="
START=$(date +%s)
buoyantPimpleFoam > log.solver 2>&1 || true
END=$(date +%s)
echo "Elapsed: $((END - START)) sec"

echo ""
echo "--- Errors ---"
grep -i "fatal\|abort" log.solver | head -3 || echo "No fatal errors"

echo ""
echo "--- deltaT ---"
grep "deltaT = " log.solver | tail -3

echo ""
echo "--- Time ---"
grep "^Time = " log.solver | tail -5

echo ""
echo "--- Probes ---"
PROBE_DIR=$(ls -d postProcessing/probes/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
if [ -n "$PROBE_DIR" ]; then
    echo "# upper_bench, lower_bench, floor_level"
    tail -20 "$PROBE_DIR/T" 2>/dev/null
fi

echo ""
echo "--- wallHeatFlux ---"
if [ -d postProcessing/wallHeatFlux ]; then
    WHFLUX_DIR=$(ls -d postProcessing/wallHeatFlux/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
    [ -n "$WHFLUX_DIR" ] && tail -7 "$WHFLUX_DIR/wallHeatFlux.dat" 2>/dev/null
fi

echo ""
echo "--- volAverageT ---"
if [ -d postProcessing/volAverageT ]; then
    VAT_DIR=$(ls -d postProcessing/volAverageT/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
    [ -n "$VAT_DIR" ] && tail -10 "$VAT_DIR/volFieldValue.dat" 2>/dev/null
fi

echo ""
echo "--- Copy back ---"
for d in [0-9]*; do cp -r "$d" "$SRC/" 2>/dev/null || true; done
cp -r postProcessing "$SRC/" 2>/dev/null || true
cp log.* "$SRC/" 2>/dev/null || true
echo "=== L-2 DONE ==="
