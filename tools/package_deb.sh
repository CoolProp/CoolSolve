#!/bin/bash
#
# package_deb.sh -- build a Debian package (.deb) of the CoolSolve desktop application.
#
#   tools/package_deb.sh [binary] [output-dir]
#
#   binary      a GUI build of coolsolve (default: build/coolsolve)
#   output-dir  where the .deb is written (default: the directory of the binary)
#
# The package installs:
#   /usr/bin/coolsolve                         the binary (stripped), CLI and GUI
#   /usr/bin/coolsolve-gui                     the desktop launcher (coolsolve --gui, usage log in ~/.local/state)
#   /usr/share/applications/coolsolve.desktop  the menu entry, opens .eescode files
#   /usr/share/mime/packages/coolsolve.xml     the MIME type of .eescode files (text/x-eescode)
#   /usr/share/icons/hicolor/256x256/apps/coolsolve.png
#
# The version is the one printed by `coolsolve --version`; the dependencies are computed by
# dpkg-shlibdeps from the binary, so the package requires the glibc of the machine that built it
# (a build on Ubuntu 24.04 does not install on Ubuntu 22.04: build on the oldest target, see
# docs/desktop_app_options.md). Needs dpkg-dev and fakeroot (sudo apt install dpkg-dev fakeroot).
# No root rights are needed. See docs/desktop_app_options.md, "Linux: a .deb package".

set -euo pipefail
umask 022

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BINARY="$(realpath "${1:-$PROJECT_DIR/build/coolsolve}")"
OUT_DIR="$(realpath "${2:-$(dirname "$BINARY")}")"

for tool in dpkg-deb dpkg-shlibdeps fakeroot strip; do
    command -v "$tool" >/dev/null || { echo "Error: $tool not found (sudo apt install dpkg-dev fakeroot binutils)" >&2; exit 1; }
done
[ -x "$BINARY" ] || { echo "Error: $BINARY is not an executable" >&2; exit 1; }

VERSION_TEXT="$("$BINARY" --version)"
VERSION="$(echo "$VERSION_TEXT" | sed -n 's/^CoolSolve \([0-9][0-9.]*\).*/\1/p' | head -1)"
[ -n "$VERSION" ] || { echo "Error: no version in the output of $BINARY --version" >&2; exit 1; }
if ! echo "$VERSION_TEXT" | grep -q "^CoolSolve Library"; then
    echo "Warning: $BINARY embeds no model library: is it a GUI build (COOLSOLVE_BUILD_GUI=ON, gui/dist built)?" >&2
fi
ARCH="$(dpkg --print-architecture)"
MAINTAINER="Sylvain Quoilin <squoilin@ulg.ac.be>"
PACKAGE="coolsolve_${VERSION}_${ARCH}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
ROOT="$WORK/$PACKAGE"

# --- Files ---
install -d "$ROOT/DEBIAN" "$ROOT/usr/bin" "$ROOT/usr/share/applications" "$ROOT/usr/share/mime/packages" \
    "$ROOT/usr/share/icons/hicolor/256x256/apps" "$ROOT/usr/share/doc/coolsolve"
install -m 755 "$BINARY" "$ROOT/usr/bin/coolsolve"
strip --strip-unneeded "$ROOT/usr/bin/coolsolve"

# The launcher of the menu entry: the usage log of the GUI goes to ~/.local/state/coolsolve instead of
# the working directory (the home folder when started from the menu).
cat > "$ROOT/usr/bin/coolsolve-gui" <<'EOF'
#!/bin/sh
# coolsolve-gui -- start the CoolSolve GUI (coolsolve --gui) and open it in the default browser.
# Arguments: an optional .eescode file to open. Close this window (or Ctrl+C) to stop CoolSolve.
state="${XDG_STATE_HOME:-$HOME/.local/state}/coolsolve"
mkdir -p "$state" 2>/dev/null
if [ -z "$COOLSOLVE_GUI_LOG" ]; then
    COOLSOLVE_GUI_LOG="$state/coolsolve_gui.log"
    export COOLSOLVE_GUI_LOG
fi
exec /usr/bin/coolsolve --gui "$@"
EOF
chmod 755 "$ROOT/usr/bin/coolsolve-gui"

# Terminal=true: the server runs in a terminal window, closed to stop CoolSolve (as the console
# window of the Windows installer).
cat > "$ROOT/usr/share/applications/coolsolve.desktop" <<'EOF'
[Desktop Entry]
Type=Application
Name=CoolSolve
GenericName=Equation Solver
Comment=Engineering equation solver compatible with EES, with CoolProp properties
Exec=coolsolve-gui %f
Icon=coolsolve
Terminal=true
Categories=Science;Engineering;
MimeType=text/x-eescode;
Keywords=EES;thermodynamics;CoolProp;equations;solver;
EOF

cat > "$ROOT/usr/share/mime/packages/coolsolve.xml" <<'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="text/x-eescode">
    <comment>CoolSolve model</comment>
    <sub-class-of type="text/plain"/>
    <glob pattern="*.eescode"/>
  </mime-type>
</mime-info>
EOF

install -m 644 "$PROJECT_DIR/docs/icon_coolsolve_256x256.png" "$ROOT/usr/share/icons/hicolor/256x256/apps/coolsolve.png"
{
    echo "CoolSolve is released under the MIT License:"
    echo
    cat "$PROJECT_DIR/LICENSE"
} > "$ROOT/usr/share/doc/coolsolve/copyright"

# --- Dependencies (from the shared libraries the binary links) ---
mkdir -p "$WORK/shlibs/debian"
echo "Source: coolsolve" > "$WORK/shlibs/debian/control"
DEPENDS="$(cd "$WORK/shlibs" && dpkg-shlibdeps -O "$ROOT/usr/bin/coolsolve" 2>/dev/null | sed -n 's/^shlibs:Depends=//p')"
[ -n "$DEPENDS" ] || { echo "Error: dpkg-shlibdeps found no dependencies" >&2; exit 1; }

INSTALLED_SIZE="$(du -sk --exclude=DEBIAN "$ROOT" | cut -f1)"
cat > "$ROOT/DEBIAN/control" <<EOF
Package: coolsolve
Version: $VERSION
Architecture: $ARCH
Maintainer: $MAINTAINER
Installed-Size: $INSTALLED_SIZE
Depends: $DEPENDS
Recommends: xdg-utils, fonts-dejavu-core
Section: science
Priority: optional
Homepage: https://github.com/CoolProp/CoolSolve
Description: Engineering equation solver compatible with EES
 CoolSolve parses, analyses and solves systems of equations written in a
 language compatible with EES (Engineering Equation Solver), with the
 thermodynamic properties of CoolProp. The GUI (coolsolve --gui) runs in the
 web browser and embeds the CoolSolve Library of models.
EOF

mkdir -p "$OUT_DIR"
fakeroot dpkg-deb --build -Zxz "$ROOT" "$OUT_DIR/$PACKAGE.deb" >/dev/null
echo "Package: $OUT_DIR/$PACKAGE.deb"
echo "Depends: $DEPENDS"
echo "Install: sudo apt install $OUT_DIR/$PACKAGE.deb"
