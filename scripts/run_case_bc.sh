#!/bin/bash
# Run Case B (no vent) and Case C (with vent) in sequence
set -e
export FOAM_SIGFPE=false

run_case() {
    local SRC="$1"
    local LABEL="$2"
    local USE_TOPOSET="$3"
    local CASE_NAME=$(basename "$SRC")
    local CASE="$HOME/saunaflow_run/$CASE_NAME"

    echo "============================================"
    echo "=== $LABEL: $CASE_NAME ==="
    echo "============================================"

    rm -rf "$CASE"
    mkdir -p "$(dirname "$CASE")"
    cp -r "$SRC" "$CASE"
    cd "$CASE"

    # Fix Windows line endings
    find . -type f ! -path '*/polyMesh/*' | xargs -r sed -i 's/\r$//'

    echo "=== blockMesh ==="
    blockMesh > log.blockMesh 2>&1
    echo "blockMesh done, exit=$?"

    if [ "$USE_TOPOSET" = "yes" ] && [ -f system/topoSetDict ]; then
        echo "=== topoSet ==="
        topoSet > log.topoSet 2>&1
        echo "topoSet done, exit=$?"
    fi

    SOLVER=$(grep "^application" system/controlDict | awk '{print $2}' | tr -d ';\r\n')
    echo "=== $SOLVER ==="
    $SOLVER > log.solver 2>&1 || true

    echo ""
    echo "--- Errors ---"
    grep -i "fatal\|abort" log.solver | head -3 || echo "No fatal errors"

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
        tail -10 "$PROBE_DIR/T" 2>/dev/null
    else
        echo "No probe data"
    fi

    echo ""
    echo "--- wallHeatFlux ---"
    if [ -d postProcessing/wallHeatFlux ]; then
        WHFLUX_DIR=$(ls -d postProcessing/wallHeatFlux/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
        if [ -n "$WHFLUX_DIR" ]; then
            for f in "$WHFLUX_DIR"/*.dat; do
                echo "$(basename "$f"): $(tail -1 "$f" 2>/dev/null)"
            done
        fi
    else
        echo "No wallHeatFlux data"
    fi

    echo ""
    echo "--- volAverageT ---"
    if [ -d postProcessing/volAverageT ]; then
        VAT_DIR=$(ls -d postProcessing/volAverageT/[0-9]* 2>/dev/null | sort -t/ -k3 -n | tail -1)
        if [ -n "$VAT_DIR" ]; then
            tail -3 "$VAT_DIR/volFieldValue.dat" 2>/dev/null
        fi
    else
        echo "No volAverageT data"
    fi

    echo ""
    echo "--- Copy back ---"
    for d in [0-9]*; do cp -r "$d" "$SRC/" 2>/dev/null || true; done
    cp -r postProcessing "$SRC/" 2>/dev/null || true
    cp log.* "$SRC/" 2>/dev/null || true
    echo "=== $LABEL Done ==="
    echo ""
}

# Case B: volume source, no ventilation
run_case "/mnt/d/dev/SaunaFEM/results/dry_sauna_steady_volsource" "Case B (no vent)" "yes"

# Case C: volume source + ventilation
run_case "/mnt/d/dev/SaunaFEM/results/dry_sauna_steady_vent" "Case C (with vent)" "yes"

echo "============================================"
echo "=== ALL DONE ==="
echo "============================================"
