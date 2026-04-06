#!/bin/bash
# Run OpenFOAM case in WSL2 on Linux filesystem
# Usage: wsl -d Ubuntu -- /usr/bin/openfoam2312 bash /mnt/d/dev/SaunaFEM/scripts/run_openfoam_wsl.sh
set -e

export FOAM_SIGFPE=false

SRC="/mnt/d/dev/SaunaFEM/results/openfoam_dry"
CASE="$HOME/openfoam_dry"

echo "=== Copying case to Linux filesystem ==="
rm -rf "$CASE"
cp -r "$SRC" "$CASE"

cd "$CASE"

echo ""
echo "=== blockMesh ==="
blockMesh 2>&1 | tail -5

echo ""
echo "=== checkMesh ==="
checkMesh 2>&1 | grep -E "OK|Mesh" | tail -3

echo ""
echo "=== Configure: 600s run (P0 validation) ==="
foamDictionary system/controlDict -entry endTime -set 600
foamDictionary system/controlDict -entry writeInterval -set 100
foamDictionary system/controlDict -entry deltaT -set 0.0001
foamDictionary system/controlDict -entry maxCo -set 0.3
foamDictionary system/controlDict -entry maxDeltaT -set 0.5
echo "deltaT=0.0001, maxCo=0.3, endTime=600, nOuterCorrectors=4"

echo ""
echo "=== buoyantPimpleFoam (600s with G_b enabled) ==="
buoyantPimpleFoam 2>&1 | tail -60

echo ""
echo "=== Copying results back ==="
cp -r "$CASE"/[0-9]* "$SRC/" 2>/dev/null || true
cp -r "$CASE/postProcessing" "$SRC/" 2>/dev/null || true
# Copy compiled dynamicCode so G_b doesn't need recompile next time
cp -r "$CASE/dynamicCode" "$SRC/" 2>/dev/null || true

echo ""
echo "=== Done ==="
ls -la "$SRC"
