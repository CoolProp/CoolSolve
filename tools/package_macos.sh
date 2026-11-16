#!/bin/bash
#
# package_macos.sh -- build the macOS application bundle (CoolSolve.app) and its disk image (.dmg).
#
#   tools/package_macos.sh [binary] [output-dir]
#
#   binary      a GUI build of coolsolve made on macOS (default: build/coolsolve)
#   output-dir  where the .dmg is written (default: the directory of the binary)
#
# Runs on macOS only (sips, iconutil, codesign, hdiutil): on a Mac, or on a GitHub Actions macOS
# runner (.github/workflows/packages.yml). The bundle:
#
#   CoolSolve.app/Contents/Info.plist
#   CoolSolve.app/Contents/MacOS/coolsolve-launcher       the bundle executable: opens the .command in Terminal
#   CoolSolve.app/Contents/MacOS/coolsolve                the binary, CLI and GUI
#   CoolSolve.app/Contents/Resources/coolsolve-gui.command   coolsolve --gui, usage log in ~/Library/Logs/CoolSolve
#   CoolSolve.app/Contents/Resources/coolsolve.icns
#
# As the console window of the Windows installer, the server runs in a Terminal window: closing it stops
# CoolSolve. The disk image is named CoolSolve_v<MAJOR.MINOR>_macOS_<arm64|x86_64|universal>.dmg.
#
# Signing: without CODESIGN_IDENTITY the bundle is signed ad hoc (enough to run on Apple Silicon, but
# Gatekeeper blocks it on the first opening: see docs/desktop_app_options.md). With
#   CODESIGN_IDENTITY="Developer ID Application: <name> (<team>)"
# it is signed with that certificate and the hardened runtime, and with NOTARY_PROFILE=<profile> (made by
# `xcrun notarytool store-credentials`) the disk image is also notarized and stapled.

set -euo pipefail
umask 022

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BINARY="${1:-$PROJECT_DIR/build/coolsolve}"
OUT_DIR="${2:-$(dirname "$BINARY")}"

[ "$(uname -s)" = "Darwin" ] || { echo "Error: run this script on macOS" >&2; exit 1; }
[ -x "$BINARY" ] || { echo "Error: $BINARY is not an executable" >&2; exit 1; }
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

VERSION_TEXT="$("$BINARY" --version)"
VERSION="$(echo "$VERSION_TEXT" | sed -n 's/^CoolSolve \([0-9][0-9.]*\).*/\1/p' | head -1)"
[ -n "$VERSION" ] || { echo "Error: no version in the output of $BINARY --version" >&2; exit 1; }
SHORT_VERSION="$(echo "$VERSION" | cut -d. -f1-2)"
if ! echo "$VERSION_TEXT" | grep -q "^CoolSolve Library"; then
    echo "Warning: $BINARY embeds no model library: is it a GUI build (COOLSOLVE_BUILD_GUI=ON, gui/dist built)?" >&2
fi
case "$(lipo -archs "$BINARY")" in
    arm64) ARCH=arm64 ;;
    x86_64) ARCH=x86_64 ;;
    *) ARCH=universal ;;
esac

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
APP="$WORK/dmg/CoolSolve.app"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"

# --- Binary and launchers ---
cp "$BINARY" "$APP/Contents/MacOS/coolsolve"
strip -x "$APP/Contents/MacOS/coolsolve"

# A file name distinct from coolsolve even on a case-insensitive file system
cat > "$APP/Contents/MacOS/coolsolve-launcher" <<'EOF'
#!/bin/sh
# Opens the CoolSolve server in a Terminal window (close the window to stop CoolSolve).
resources="$(cd "$(dirname "$0")/../Resources" && pwd)"
exec open -a Terminal "$resources/coolsolve-gui.command"
EOF

cat > "$APP/Contents/Resources/coolsolve-gui.command" <<'EOF'
#!/bin/sh
# Starts the CoolSolve GUI (coolsolve --gui) and opens it in the default browser.
# Close this window (or Ctrl+C) to stop CoolSolve.
macos="$(cd "$(dirname "$0")/../MacOS" && pwd)"
logdir="$HOME/Library/Logs/CoolSolve"
mkdir -p "$logdir" 2>/dev/null
if [ -z "$COOLSOLVE_GUI_LOG" ]; then
    COOLSOLVE_GUI_LOG="$logdir/coolsolve_gui.log"
    export COOLSOLVE_GUI_LOG
fi
cd "$HOME"
exec "$macos/coolsolve" --gui
EOF
chmod 755 "$APP/Contents/MacOS/coolsolve-launcher" "$APP/Contents/Resources/coolsolve-gui.command"

# --- Icon (from the 256 px PNG) ---
ICONSET="$WORK/coolsolve.iconset"
mkdir -p "$ICONSET"
ICON_PNG="$PROJECT_DIR/docs/icon_coolsolve_256x256.png"
for size in 16 32 128 256; do
    sips -z "$size" "$size" "$ICON_PNG" --out "$ICONSET/icon_${size}x${size}.png" >/dev/null
done
sips -z 32 32 "$ICON_PNG" --out "$ICONSET/icon_16x16@2x.png" >/dev/null
sips -z 64 64 "$ICON_PNG" --out "$ICONSET/icon_32x32@2x.png" >/dev/null
cp "$ICON_PNG" "$ICONSET/icon_128x128@2x.png"
iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/coolsolve.icns"

cat > "$APP/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key>                 <string>CoolSolve</string>
    <key>CFBundleDisplayName</key>          <string>CoolSolve</string>
    <key>CFBundleIdentifier</key>           <string>be.uliege.coolsolve</string>
    <key>CFBundleExecutable</key>           <string>coolsolve-launcher</string>
    <key>CFBundleIconFile</key>             <string>coolsolve</string>
    <key>CFBundlePackageType</key>          <string>APPL</string>
    <key>CFBundleShortVersionString</key>   <string>$VERSION</string>
    <key>CFBundleVersion</key>              <string>$VERSION</string>
    <key>LSMinimumSystemVersion</key>       <string>12.0</string>
    <key>NSHumanReadableCopyright</key>     <string>MIT License</string>
</dict>
</plist>
EOF
plutil -lint "$APP/Contents/Info.plist" >/dev/null

# --- Signing: the binary first, then the bundle ---
if [ -n "${CODESIGN_IDENTITY:-}" ]; then
    SIGN=(codesign --force --options runtime --timestamp --sign "$CODESIGN_IDENTITY")
else
    SIGN=(codesign --force --sign -)
fi
"${SIGN[@]}" "$APP/Contents/MacOS/coolsolve"
"${SIGN[@]}" "$APP"
codesign --verify --strict "$APP"

# --- Disk image: the app and a link to /Applications ---
ln -s /Applications "$WORK/dmg/Applications"
DMG="$OUT_DIR/CoolSolve_v${SHORT_VERSION}_macOS_${ARCH}.dmg"
rm -f "$DMG"
hdiutil create -volname "CoolSolve $SHORT_VERSION" -srcfolder "$WORK/dmg" -ov -format UDZO "$DMG" >/dev/null

if [ -n "${CODESIGN_IDENTITY:-}" ]; then
    codesign --force --timestamp --sign "$CODESIGN_IDENTITY" "$DMG"
    if [ -n "${NOTARY_PROFILE:-}" ]; then
        xcrun notarytool submit "$DMG" --keychain-profile "$NOTARY_PROFILE" --wait
        xcrun stapler staple "$DMG"
    fi
fi

echo "Disk image: $DMG"
