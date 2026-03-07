#!/bin/bash
# build_gui.sh — Rebuild the GUI frontend and the CoolSolve binary with embedded assets.
# Run from any directory; the script locates paths relative to itself.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GUI_DIR="$SCRIPT_DIR/gui"
BUILD_DIR="$SCRIPT_DIR/build"

echo "=== 1/3  Building GUI frontend ==="
cd "$GUI_DIR"
npm install --silent
npm run build
echo "    GUI built → $GUI_DIR/dist/"

echo ""
echo "=== 2/3  Configuring CMake (GUI=ON) ==="
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake -DCOOLSOLVE_BUILD_GUI=ON ..

echo ""
echo "=== 3/3  Compiling CoolSolve ==="
make -j"$(nproc)"

echo ""
echo "=== Done. Binary: $BUILD_DIR/coolsolve ==="
