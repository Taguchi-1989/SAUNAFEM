#!/bin/bash
# Run Case F: transient + ventilation + viewFactor radiation
set -e
export FOAM_SIGFPE=false

SRC="/mnt/d/dev/SaunaFEM/results/dry_sauna_transient_rad"
CASE_NAME="dry_sauna_transient_rad"
CASE="$HOME/saunaflow_run/$CASE_NAME"

echo "============================================"
echo "=== Case F (radiation): $CASE_NAME ==="
echo "============================================"

rm -rf "$CASE"
mkdir -p "$(dirname "$CASE")"
cp -r "$SRC" "$CASE"
cd "$CASE"

# Fix Windows CRLF only in OpenFOAM case files
for f in $(find 0 constant system -type f 2>/dev/null); do
    sed -i 's/\r$//' "$f"
done

echo "=== blockMesh ==="
blockMesh > log.blockMesh 2>&1
echo "blockMesh done, exit=$?"

if [ -f system/topoSetDict ]; then
    echo "=== topoSet ==="
    topoSet > log.topoSet 2>&1
    echo "topoSet done, exit=$?"
fi

# viewFactor generation: faceAgglomerate -> viewFactorsGen
if grep -q "radiationModel.*viewFactor" constant/radiationProperties 2>/dev/null; then
    echo "=== faceAgglomerate ==="
    faceAgglomerate > log.faceAgglomerate 2>&1 || true
    echo "faceAgglomerate exit=$?"
    grep -i "error\|fatal" log.faceAgglomerate | head -3 || echo "No faceAgglomerate errors"

    echo "=== viewFactorsGen ==="
    viewFactorsGen > log.viewFactorsGen 2>&1 || true
    echo "viewFactorsGen exit=$?"
    grep -i "error\|fatal" log.viewFactorsGen | head -3 || echo "No viewFactor errors"
    grep "coarse faces" log.viewFactorsGen || true
fi

SOLVER=$(grep "^application" system/controlDict | awk '{print $2}' | tr -d ';\r\n')
echo "=== $SOLVER (600s physical time) ==="
$SOLVER > log.solver 2>&1 || true

echo ""
echo "--- Errors ---"
grep -i "fatal\|abort" log.solver | head -5 || echo "No fatal errors"

echo ""
echo "--- Progress ---"
grep "^Time = " log.solver | tail -5

echo ""
echo "--- rho ---"
grep "rho min/max" log.solver | tail -3

echo ""
echo "--- Probes ---"
PROBE_DIR=$(ls -d postProcessing/probes/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
if [ -n "$PROBE_DIR" ]; then
    cat "$PROBE_DIR/T" 2>/dev/null
else
    echo "No probe data"
fi

echo ""
echo "--- volAverageT ---"
if [ -d postProcessing/volAverageT ]; then
    VAT_DIR=$(ls -d postProcessing/volAverageT/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
    if [ -n "$VAT_DIR" ]; then
        cat "$VAT_DIR/volFieldValue.dat" 2>/dev/null
    fi
else
    echo "No volAverageT data"
fi

echo ""
echo "--- Copy back ---"
for d in [0-9]*; do cp -r "$d" "$SRC/" 2>/dev/null || true; done
cp -r postProcessing "$SRC/" 2>/dev/null || true
cp log.* "$SRC/" 2>/dev/null || true
echo "=== Case F Done ==="
