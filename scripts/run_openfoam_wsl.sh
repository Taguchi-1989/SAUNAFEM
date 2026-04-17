#!/bin/bash
# Generic OpenFOAM runner for SaunaFlow cases via WSL
# Usage: run_openfoam_wsl.sh <case_dir> [--toposet]
#
# Copies the case to Linux filesystem for performance, runs solver,
# then copies results back to the Windows source directory.
set -e
export FOAM_SIGFPE=false

SRC="${1:-/mnt/d/dev/SaunaFEM/results/openfoam_dry}"
RUN_TOPOSET=false

for arg in "$@"; do
    if [ "$arg" = "--toposet" ]; then
        RUN_TOPOSET=true
    fi
done

CASE_NAME=$(basename "$SRC")
WORK_DIR="$HOME/saunaflow_run/$CASE_NAME"

echo "=== Copying case to Linux FS ==="
rm -rf "$WORK_DIR"
mkdir -p "$(dirname "$WORK_DIR")"
cp -r "$SRC" "$WORK_DIR"
cd "$WORK_DIR"

# Detect solver from controlDict
SOLVER=$(grep "^application" system/controlDict | awk '{print $2}' | tr -d ';')
echo "=== Case: $CASE_NAME | Solver: $SOLVER ==="

echo "=== blockMesh ==="
blockMesh > log.blockMesh 2>&1

if [ "$RUN_TOPOSET" = true ] && [ -f system/topoSetDict ]; then
    echo "=== topoSet (heaterZone) ==="
    topoSet > log.topoSet 2>&1
fi

echo "=== $SOLVER ==="
$SOLVER > log.solver 2>&1 || true

echo ""
echo "=== Errors? ==="
grep -i "fatal\|abort\|Negative" log.solver | head -3 || echo "No fatal errors"

echo ""
echo "=== Progress ==="
grep "^Time = " log.solver | tail -5

echo ""
echo "=== Courant ==="
grep "Courant Number" log.solver | tail -3

echo ""
echo "=== rho range ==="
grep "rho min/max" log.solver | tail -3

echo ""
echo "=== Probe temps ==="
PROBE_DIR=$(ls -d postProcessing/probes/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
if [ -n "$PROBE_DIR" ]; then
    echo "Probe dir: $PROBE_DIR"
    tail -10 "$PROBE_DIR/T" 2>/dev/null
fi

echo ""
echo "=== Heat balance ==="
if [ -d postProcessing/wallHeatFlux ]; then
    WHFLUX_DIR=$(ls -d postProcessing/wallHeatFlux/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
    if [ -n "$WHFLUX_DIR" ]; then
        echo "wallHeatFlux dir: $WHFLUX_DIR"
        for f in "$WHFLUX_DIR"/*.dat; do
            echo "--- $(basename "$f") ---"
            tail -3 "$f" 2>/dev/null
        done
    fi
fi

if [ -d postProcessing/volAverageT ]; then
    VAT_DIR=$(ls -d postProcessing/volAverageT/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
    if [ -n "$VAT_DIR" ]; then
        echo "volAverageT:"
        tail -3 "$VAT_DIR/volFieldValue.dat" 2>/dev/null
    fi
fi

echo ""
echo "=== Copying results back ==="
# Copy time directories, postProcessing, and logs back to source
for d in [0-9]*; do
    cp -r "$d" "$SRC/" 2>/dev/null || true
done
cp -r postProcessing "$SRC/" 2>/dev/null || true
cp log.* "$SRC/" 2>/dev/null || true
echo "=== Done ==="
