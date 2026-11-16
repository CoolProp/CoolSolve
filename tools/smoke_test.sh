#!/bin/bash
#
# smoke_test.sh -- quick check of a packaged coolsolve binary (CI of the release packages, or by hand).
#
#   tools/smoke_test.sh <binary> [port]
#
# 1. `coolsolve --version` names the embedded model library (a GUI build);
# 2. an example with CoolProp properties is solved by the CLI (in a temporary copy: nothing written in examples/);
# 3. the GUI server starts with --no-browser, CoolProp warms up, and the library menu is served.
# Needs bash and curl (Git Bash on Windows). See .github/workflows/packages.yml.

set -euo pipefail

BINARY="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
PORT="${2:-18650}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
SERVER=""
cleanup() {
    [ -n "$SERVER" ] && kill "$SERVER" 2>/dev/null || true
    rm -rf "$WORK"
}
trap cleanup EXIT

echo "=== 1/3 Version ==="
"$BINARY" --version | tee "$WORK/version.txt"
grep -q "^CoolSolve Library" "$WORK/version.txt" || { echo "FAIL: no embedded model library"; exit 1; }

echo "=== 2/3 CLI solve ==="
cp "$PROJECT_DIR/examples/compressor_refrigeration_simple.eescode" "$WORK/"
(cd "$WORK" && "$BINARY" --no-sol compressor_refrigeration_simple.eescode) | tee "$WORK/solve.txt"
grep -q "Solver: SUCCESS" "$WORK/solve.txt" || { echo "FAIL: the example did not solve"; exit 1; }

echo "=== 3/3 GUI server ==="
(cd "$WORK" && COOLSOLVE_GUI_LOG="$WORK/gui.log" exec "$BINARY" --gui "$PORT" --no-browser) > "$WORK/server.txt" 2>&1 &
SERVER=$!
ready=no
for _ in $(seq 1 120); do
    if curl -sf -m 2 "http://127.0.0.1:$PORT/api/v1/health" | grep -q '"coolpropReady":true'; then
        ready=yes
        break
    fi
    sleep 1
done
[ "$ready" = yes ] || { echo "FAIL: the GUI server did not become ready"; cat "$WORK/server.txt"; exit 1; }
curl -sf -m 10 "http://127.0.0.1:$PORT/api/v1/library/menu" > "$WORK/menu.json"
[ "$(wc -c < "$WORK/menu.json")" -gt 1000 ] || { echo "FAIL: empty library menu"; exit 1; }
curl -sf -m 10 -o /dev/null "http://127.0.0.1:$PORT/" || { echo "FAIL: the GUI page is not served"; exit 1; }
echo "GUI server: ready, library menu $(wc -c < "$WORK/menu.json") bytes"
echo "=== Smoke test passed ==="
