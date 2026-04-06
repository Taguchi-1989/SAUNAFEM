#!/bin/bash
# Setup OpenFOAM 2312 on WSL2 Ubuntu
# Usage: wsl -d Ubuntu -- bash /mnt/d/dev/SaunaFEM/scripts/setup_openfoam_wsl.sh

set -e

echo "=== Checking existing OpenFOAM installation ==="
if command -v blockMesh &>/dev/null; then
    echo "OpenFOAM already installed:"
    blockMesh -help 2>&1 | head -3
    exit 0
fi

echo "=== Adding OpenFOAM repository ==="
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash

echo "=== Installing OpenFOAM 2312 ==="
sudo apt-get update
sudo apt-get install -y openfoam2312

echo "=== Adding to bashrc ==="
echo 'source /usr/lib/openfoam/openfoam2312/etc/bashrc' >> ~/.bashrc

echo "=== Verifying ==="
source /usr/lib/openfoam/openfoam2312/etc/bashrc
which blockMesh && echo "SUCCESS: OpenFOAM installed"
blockMesh -help 2>&1 | head -3
