#!/bin/bash
# Run Case L-1: buoyantPimpleFoam transient 3600s with adjustTimeStep
set -e
source /usr/lib/openfoam/openfoam2312/etc/bashrc 2>/dev/null || true
export FOAM_SIGFPE=false

BASE="/mnt/d/dev/SaunaFEM"
SRC="$BASE/results/L1"
CASE_NAME="L1"
CASE="$HOME/saunaflow_run/$CASE_NAME"

echo "============================================"
echo "=== L-1: 13kW transient (3600s) ==="
echo "============================================"

rm -rf "$CASE"
mkdir -p "$(dirname "$CASE")"
cp -r "$SRC" "$CASE"
cd "$CASE"

find . -type f ! -path '*/polyMesh/*' | xargs -r sed -i 's/\r$//'

echo "=== blockMesh ==="
blockMesh > log.blockMesh 2>&1
echo "blockMesh done, exit=$?"

if [ -f system/topoSetDict ]; then
    echo "=== topoSet ==="
    topoSet > log.topoSet 2>&1
    echo "topoSet done, exit=$?"
fi

if grep -q "radiationModel.*viewFactor" constant/radiationProperties 2>/dev/null; then
    echo "=== faceAgglomerate ==="
    faceAgglomerate -dict constant/viewFactorsDict > log.faceAgglomerate 2>&1
    echo "faceAgglomerate done, exit=$?"
    echo "=== viewFactorsGen ==="
    viewFactorsGen > log.viewFactorsGen 2>&1
    echo "viewFactorsGen done, exit=$?"
fi

echo "=== buoyantPimpleFoam (3600s, adjustTimeStep) ==="
START=$(date +%s)
buoyantPimpleFoam > log.solver 2>&1 || true
END=$(date +%s)
echo "Elapsed: $((END - START)) sec"

echo ""
echo "--- Errors ---"
grep -i "fatal\|abort" log.solver | head -3 || echo "No fatal errors"

echo ""
echo "--- deltaT progression ---"
grep "deltaT = " log.solver | head -5
echo "..."
grep "deltaT = " log.solver | tail -5

echo ""
echo "--- Time progression ---"
grep "^Time = " log.solver | tail -10

echo ""
echo "--- Probes ---"
PROBE_DIR=$(ls -d postProcessing/probes/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
if [ -n "$PROBE_DIR" ]; then
    echo "# upper_bench, lower_bench, floor_level"
    tail -20 "$PROBE_DIR/T" 2>/dev/null
else
    echo "No probe data"
fi

echo ""
echo "--- wallHeatFlux ---"
if [ -d postProcessing/wallHeatFlux ]; then
    WHFLUX_DIR=$(ls -d postProcessing/wallHeatFlux/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
    if [ -n "$WHFLUX_DIR" ]; then
        tail -7 "$WHFLUX_DIR/wallHeatFlux.dat" 2>/dev/null
    fi
fi

echo ""
echo "--- volAverageT ---"
if [ -d postProcessing/volAverageT ]; then
    VAT_DIR=$(ls -d postProcessing/volAverageT/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
    if [ -n "$VAT_DIR" ]; then
        tail -10 "$VAT_DIR/volFieldValue.dat" 2>/dev/null
    fi
fi

echo ""
echo "--- Copy results back ---"
for d in [0-9]*; do cp -r "$d" "$SRC/" 2>/dev/null || true; done
cp -r postProcessing "$SRC/" 2>/dev/null || true
cp log.* "$SRC/" 2>/dev/null || true
echo "=== L-1 DONE ==="
