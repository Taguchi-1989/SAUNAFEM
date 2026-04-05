#!/bin/bash
# Run OpenFOAM case in WSL2
# Usage: wsl -d Ubuntu -- /usr/bin/openfoam2312 bash /mnt/d/dev/SaunaFEM/scripts/run_openfoam_wsl.sh
set -e

export FOAM_SIGFPE=false

# Copy case to Linux filesystem (NTFS /mnt/d can't compile codedSource)
SRC="/mnt/d/dev/SaunaFEM/results/openfoam_dry"
CASE="$HOME/openfoam_dry"

echo "=== Copying case to Linux filesystem ==="
rm -rf "$CASE"
cp -r "$SRC" "$CASE"

# Disable codedSource (wmake fails on NTFS, G_b is Phase 2 refinement)
cat > "$CASE/constant/fvOptions" << 'ENDOPTIONS'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvOptions;
}
// G_b buoyancy production disabled (codedSource requires Linux FS for wmake)
ENDOPTIONS

cd "$CASE"

echo ""
echo "=== blockMesh ==="
blockMesh 2>&1 | tail -5

echo ""
echo "=== Configure: 30s run ==="
foamDictionary system/controlDict -entry endTime -set 30
foamDictionary system/controlDict -entry writeInterval -set 10
foamDictionary system/controlDict -entry deltaT -set 0.001
foamDictionary system/controlDict -entry maxCo -set 0.3
foamDictionary system/controlDict -entry maxDeltaT -set 0.5

echo ""
echo "=== buoyantPimpleFoam (30s) ==="
buoyantPimpleFoam 2>&1 | tail -40

echo ""
echo "=== Copying results back ==="
cp -r "$CASE"/[0-9]* "$SRC/" 2>/dev/null || true
cp -r "$CASE/postProcessing" "$SRC/" 2>/dev/null || true

echo ""
echo "=== Done ==="
ls -la "$SRC"
