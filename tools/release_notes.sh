#!/bin/bash
#
# release_notes.sh -- print the notes of a GitHub release from docs/versions.md.
#
#   tools/release_notes.sh vX.Y > notes.md
#   gh release create vX.Y --title "CoolSolve vX.Y" --notes-file notes.md <packages>
#
# The notes are the "## vX.Y — ..." section of docs/versions.md (up to the next "---"), its links to the other
# documents of docs/ made absolute (https://github.com/CoolProp/CoolSolve/blob/vX.Y/docs/...): relative links do not
# resolve on a release page. See docs/contributing.md §9.7.

set -euo pipefail

TAG="${1:?usage: tools/release_notes.sh vX.Y}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
REPO_URL="https://github.com/CoolProp/CoolSolve"

NOTES="$(awk -v tag="$TAG" '
    index($0, "## " tag " ") == 1 { found = 1 }
    found && /^---$/ { exit }
    found { print }
' "$PROJECT_DIR/docs/versions.md")"
[ -n "$NOTES" ] || { echo "Error: no section \"## $TAG ...\" in docs/versions.md" >&2; exit 1; }

# The heading is the release title already; links: a.md, a.md#anchor -> $REPO_URL/blob/$TAG/docs/a.md...
echo "$NOTES" | tail -n +2 \
    | sed -E "s#\]\(([A-Za-z0-9_.-]+\.md)#](${REPO_URL}/blob/${TAG}/docs/\1#g; s#\]\(\.\./README\.md#](${REPO_URL}/blob/${TAG}/README.md#g"
