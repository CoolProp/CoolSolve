# CoolSolve Desktop Application — Installers and Options

[← Back to overview](../README.md)

This document covers how CoolSolve is distributed as a desktop application:

- §1: the packages of each platform;
- §2: building them for a release (the release process itself is [Contributing §9](contributing.md#9-release-process));
- §3: Windows (`.exe` installer);
- §4: macOS (`.dmg` disk image);
- §5: Linux (`.deb` package);
- §6: the issues common to the desktop packages;
- §7: the evaluation of the options for a native application window (not implemented).

---

## 1. Overview

The desktop application is the `coolsolve` binary of a GUI build (`COOLSOLVE_BUILD_GUI=ON`, `gui/dist` built
first): one self-contained executable that embeds the React GUI, the CoolSolve Library and CoolProp.
`coolsolve --gui` starts the embedded HTTP server on `127.0.0.1:8550` (the next free port when it is taken) and
opens the default browser. The installers only add a launcher, an icon and, where possible, the `.eescode` file
association around that binary.

| Platform | Package | Packaged by | Signed | Since |
|----------|---------|-------------|--------|-------|
| Windows 10/11 x64 | `CoolSolve_vX.Y_Installer.exe` (NSIS) | `coolsolve.nsi` | no | v0.1 |
| macOS 12+, Apple Silicon and Intel | `CoolSolve_vX.Y_macOS_universal.dmg` | `tools/package_macos.sh` | ad hoc | v0.5 |
| Ubuntu 22.04+ / Debian 12+, x64 | `coolsolve_X.Y.0_amd64.deb` | `tools/package_deb.sh` | no | v0.5 |

All three are built by the GitHub Actions workflow **`.github/workflows/packages.yml`** (*Release packages*), from one
commit, on GitHub-hosted machines: no Mac and no Windows machine is needed. The Windows installer can still be built
by hand with `build_installer.bat` (§3), and the `.deb` on any Debian/Ubuntu machine (§5.1).

The window is the user's browser in all three cases (the console or Terminal window of the server is the way to
stop it, §6). A real application window is the subject of §7.

---

## 2. Building the packages of a release

The steps of a release (version strings, library tag, tests, changelog, tag, GitHub release) are
[Contributing §9](contributing.md#9-release-process). The packages:

```bash
git push origin main                                    # the workflow builds a commit that is on GitHub
gh workflow run packages.yml --ref main -f library_tag=vX.Y
gh run list --workflow packages.yml --limit 1           # the id of the run
gh run watch <run-id>                                   # 30-90 min; or the Actions tab of the repository
gh run view <run-id> --log-failed                       # when a job fails
gh run download <run-id> --dir dist                     # dist/coolsolve-{windows,macos,linux}/<package>
```

| Job | Runner | Steps | Artifact |
|-----|--------|-------|----------|
| `windows` | `windows-latest` (Visual Studio, CMake, Git, Python, Node) | GUI, CMake (`-A x64`, no zlib), MSVC build with the 64-bit host compiler, smoke test, NSIS (`choco install nsis` if missing) | `coolsolve-windows` |
| `macos` | `macos-latest` (Apple Silicon) | GUI, CMake (`CMAKE_OSX_ARCHITECTURES`, macOS 12+), build, smoke test, `tools/package_macos.sh` | `coolsolve-macos` |
| `linux` | `ubuntu-22.04` (the oldest target: the `.deb` needs the glibc of its build machine) | GUI, CMake, build, smoke test, `tools/package_deb.sh` | `coolsolve-linux` |

Inputs: `library_tag` (the library tag to embed; `main` by default) and `mac_architectures` (`arm64;x86_64`, a
universal binary, by default; `arm64` builds faster). The smoke test (`tools/smoke_test.sh`, also usable by hand)
checks that `--version` names the embedded library, solves one example with the CLI, and starts the GUI server
(health, library menu, page). The workflow never runs on a push: a run costs nothing for a public repository but
takes an hour of three machines. Artifacts are kept 90 days and need a GitHub login to download.

Check, on the downloaded packages: the Linux binary (`dpkg-deb -x` the `.deb`, then `usr/bin/coolsolve --version`)
prints `CoolSolve X.Y.0` and the commit of the library tag; the run log shows the same for Windows and macOS.

---

## 3. Windows

The installer (`coolsolve.nsi`) installs `coolsolve.exe` in `Program Files\CoolSolve`, a Start Menu shortcut, the
`.eescode` association (`coolsolve.exe --gui "%1"`) and an uninstaller. It is not signed: Windows SmartScreen
shows *Windows protected your PC* on the first run (*More info* → *Run anyway*).

The `windows` job of §2 builds it. By hand, on a Windows machine with Visual Studio 2022, Git, Python 3, Node.js and
NSIS (README, *Building on Windows*):

```bat
set "LIBRARY_TAG=vX.Y"
build_installer.bat
```

writes `CoolSolve_vX.Y_Installer.exe` in the project root. Both pass `PreferredToolArchitecture=x64`: the embedded
library is about 50 MB of C arrays, on which the 32-bit host compiler runs out of heap (C1060). The CI job also
disables zlib (`CMAKE_DISABLE_FIND_PACKAGE_ZLIB`): a zlib found on the runner would make `coolsolve.exe` depend on
a DLL the installer does not ship (coolplot-core then writes uncompressed PDF streams, as in the hand build).

---

## 4. macOS

### 4.1 How to build without a Mac: the options

CoolSolve has no macOS-specific code beyond `open <url>` in `openBrowser()` (`src/server.cpp`), and CoolProp,
cpp-httplib and Eigen build with Apple Clang: the work is finding a Mac and packaging.

| Option | Cost | Effort | Verdict |
|--------|------|--------|---------|
| **GitHub Actions macOS runner** (`macos-latest`, Apple Silicon) | free for a public repository | the `macos` job of `.github/workflows/packages.yml` | **chosen**: reproducible, from the release commit, no hardware |
| A Mac (own, borrowed, a colleague's) | — | install Xcode command-line tools, Homebrew `cmake node`, then the same commands as the workflow | good for testing the app by hand; needed anyway to try the `.dmg` once |
| Rented cloud Mac (MacinCloud, AWS EC2 Mac, Scaleway Apple silicon) | hourly or monthly rental (AWS and Scaleway bill 24 h minimum) | as a Mac, through SSH/VNC | only if Actions is not an option |
| Cross-compiling from Linux (osxcross) | free | high: needs the macOS SDK extracted from Xcode (Apple licence: on Apple hardware only), no way to run the result | not recommended |

The runner is an Apple Silicon Mac. With `mac_architectures = arm64;x86_64` (the default) CMake builds a universal
binary, which also runs on Intel Macs (about twice the compile time; the file is `..._macOS_universal.dmg`, with
`arm64` alone `..._macOS_arm64.dmg`). The minimum system is macOS 12 (`CMAKE_OSX_DEPLOYMENT_TARGET`).

### 4.2 Building

The `macos` job of §2. On a Mac by hand: Xcode command-line tools, Homebrew `cmake node`, then the commands of the
job (`npm ci && npm run build` in `gui/`, the CMake configure and build, `tools/package_macos.sh build/coolsolve`).
Nobody has opened the app on a Mac yet: the smoke test of the job checks the binary (CLI and GUI server), not the
bundle, the launcher or Gatekeeper; report what a first user sees.

### 4.3 The application bundle

`tools/package_macos.sh [binary] [output-dir]` (on macOS only) builds:

```
CoolSolve.app/Contents/
├── Info.plist                       # be.uliege.coolsolve, version from coolsolve --version, macOS 12+
├── MacOS/coolsolve-launcher         # the bundle executable: open -a Terminal coolsolve-gui.command
├── MacOS/coolsolve                  # the binary (stripped), CLI and GUI
└── Resources/coolsolve-gui.command  # coolsolve --gui, usage log in ~/Library/Logs/CoolSolve/
    Resources/coolsolve.icns         # from docs/icon_coolsolve_256x256.png
```

and puts it in `CoolSolve_vX.Y_macOS_universal.dmg` with a link to `/Applications` (drag and drop install). Opening the
app opens a Terminal window running the server, and the browser on the GUI; closing the Terminal window stops
CoolSolve — the behaviour of the Windows console window. The CLI is the binary inside the bundle:

```bash
sudo ln -s /Applications/CoolSolve.app/Contents/MacOS/coolsolve /usr/local/bin/coolsolve
```

Not supported: double-clicking a `.eescode` file. macOS hands the file to an application through an Apple Event,
not on the command line, which a shell launcher cannot receive (the native window of §7 would); files are opened
with the GUI's *Open* button.

### 4.4 Signing, Gatekeeper and the first opening

Without an Apple Developer ID the bundle is signed **ad hoc** (`codesign --sign -`), which Apple Silicon requires
to run any binary but which Gatekeeper does not trust. A `.dmg` downloaded with a browser is quarantined, and the
first opening of the app is refused (*“CoolSolve” cannot be opened because Apple cannot check it for malicious
software* / *Apple could not verify…*). Since macOS 15 the right-click → *Open* bypass is gone. The note for the
release page:

> **macOS**: open the `.dmg` and drag CoolSolve to Applications. The app is not signed by Apple: on the first
> opening, macOS refuses it. Either go to *System Settings → Privacy & Security* and click *Open Anyway* under the
> CoolSolve message (repeat if macOS asks again), or run once in Terminal:
> `xattr -dr com.apple.quarantine /Applications/CoolSolve.app`

To remove the warning altogether:

1. **Apple Developer Program** (99 $/year, personal or through ULiège if the university has an organisation
   account): a *Developer ID Application* certificate.
2. Sign and notarize: `tools/package_macos.sh` already does it when the environment provides
   `CODESIGN_IDENTITY="Developer ID Application: <name> (<team id>)"` (hardened runtime, timestamp) and
   `NOTARY_PROFILE=<profile>` (created once with `xcrun notarytool store-credentials`): it signs the binary, the
   bundle and the `.dmg`, submits the `.dmg` to Apple's notary service and staples the ticket.
3. On GitHub Actions, the certificate (`.p12`) and the App Store Connect API key go in repository secrets, imported
   into a temporary keychain before the packaging step (a few lines with `security import`; not written yet).

Not tested yet; the shell launcher may need replacing by a small compiled launcher for notarization to pass.

### 4.5 The dependencies on macOS (and Windows)

| Dependency | How it is obtained | Windows (MSVC) | macOS (Apple Clang) |
|------------|-------------------|----------------|---------------------|
| CoolProp (pinned commit) | FetchContent `git clone` from GitHub, compiled from source; Python 3 at configure | works (v0.4 installer) | expected to work: CoolProp builds on macOS in its own CI; first CoolSolve build on the runner |
| CoolSolve Library | FetchContent `git clone` (public repository), embedded as C arrays by a pure CMake script (`cmake/embed_common.cmake`): no platform tool | works (64-bit compiler host needed, set by `build_installer.bat`) | expected to work: same script |
| coolplot-core, cpp-httplib, peglib, nlohmann/json, Eigen | FetchContent, header-only or plain C++17 | works | expected to work |
| Fonts of the PDF diagrams | found at run time by coolplot-core in `C:/Windows/Fonts`, `/Library/Fonts`, `/System/Library/Fonts`, `~/Library/Fonts` (TrueType) | Arial | Arial (`/System/Library/Fonts/Supplemental`), to check once |

Nothing has to be installed on the user's machine: the binary embeds everything, and the system libraries it links
(C++ runtime, zlib) are part of macOS.

---

## 5. Linux: a `.deb` package

### 5.1 Building

The release `.deb` is built by the `linux` job of §2 (Ubuntu 22.04). On any Debian/Ubuntu machine:


```bash
sudo apt install dpkg-dev fakeroot          # once (both were already installed here)
./build_gui.sh                              # or the usual cmake/make of a GUI build, with -DCOOLSOLVE_LIBRARY_TAG=vX.Y
tools/package_deb.sh                        # build/coolsolve → build/coolsolve_X.Y.0_amd64.deb
```

`tools/package_deb.sh [binary] [output-dir]` needs no root rights. It reads the version from
`coolsolve --version`, strips the binary, computes the dependencies with `dpkg-shlibdeps` and writes the package
with `fakeroot dpkg-deb`. It installs:

| Path | Content |
|------|---------|
| `/usr/bin/coolsolve` | the binary: the CLI everywhere in the `PATH`, and the GUI |
| `/usr/bin/coolsolve-gui` | the menu launcher: `coolsolve --gui [file]`, usage log in `~/.local/state/coolsolve/` |
| `/usr/share/applications/coolsolve.desktop` | the menu entry *CoolSolve* (Science), opens `.eescode` files; `Terminal=true` |
| `/usr/share/mime/packages/coolsolve.xml` | the MIME type `text/x-eescode` of `*.eescode` (a kind of text) |
| `/usr/share/icons/hicolor/256x256/apps/coolsolve.png` | the icon |
| `/usr/share/doc/coolsolve/copyright` | the MIT licence |

`Depends`, built on Pop!_OS 24.04: `libc6 (>= 2.38), libgcc-s1, libstdc++6 (>= 13.1), zlib1g`; `Recommends`:
`xdg-utils` (opens the browser) and `fonts-dejavu-core` (PDF diagrams). The menu entry and the MIME type are
registered by the triggers of `desktop-file-utils` and `shared-mime-info`: no maintainer scripts. The current build
of `main` gives a 10 MB package (35 MB installed).

### 5.2 Installing and removing

```bash
sudo apt install ./coolsolve_X.Y.0_amd64.deb      # apt, not dpkg -i: it installs the dependencies
sudo apt remove coolsolve
```

### 5.3 Which distributions it installs on

The binary requires the glibc and libstdc++ of the machine that built it. The release `.deb`, built on Ubuntu
22.04 (glibc 2.35), installs on **Ubuntu 22.04+, Debian 12+, Linux Mint 21+, Pop!_OS 22.04+**. A build on an Ubuntu
24.04 machine (glibc 2.39) installs on Ubuntu 24.04+ / Debian 13+ only, and `apt` refuses it on Ubuntu 22.04. For an
older target than the CI's, build in a container of that system, e.g. with Docker (`sudo apt install docker.io`):

```bash
docker run --rm -it -v "$PWD":/src -w /src ubuntu:22.04 bash
# inside: apt update && apt install -y build-essential cmake git python3 nodejs npm dpkg-dev fakeroot
#         (Node 18+ is needed: take it from NodeSource if the distribution's is older), then a GUI build
#         in a separate build directory (build-jammy) and tools/package_deb.sh build-jammy/coolsolve
```

or change the runner of the `linux` job. Other formats, if asked for: an **AppImage** (one file for every
distribution, from the same binary with `linuxdeploy`) or a plain `.tar.gz` of the binary; `.rpm` with CPack or
`alien` from the `.deb`.

---

## 6. Issues common to the desktop packages

These show in a desktop setting (the GUI was first written for the server):

1. **The server listens on `127.0.0.1` only** (since v0.5; `--host <address>`, `ServerOptions::host`). Up to
   v0.4 it listened on `0.0.0.0`: on a laptop the GUI — files, solves — was reachable from the local network, and
   macOS with its firewall on asked *Do you want the application “coolsolve” to accept incoming network
   connections?*. The online server is unaffected: Apache proxies to `127.0.0.1:8550` and the deploy script checks
   `http://127.0.0.1:8550/api/v1/health` (checked on labothap, 2026-11-16). `--host 0.0.0.0` restores the former
   behaviour (a container, a machine reached without reverse proxy). The Vite dev server proxies to `127.0.0.1`,
   not `localhost`, which Node 18 resolves to `::1` first.
2. **Stopping the server**: closing the browser tab leaves it running; the console (Windows), Terminal (macOS) or
   terminal (Linux, `Terminal=true`) window is the way to stop it. With `Terminal=false` the server would be
   invisible and could only be killed.
3. **Usage log**: by default `coolsolve_gui.log` in the working directory. The macOS and Linux launchers set
   `COOLSOLVE_GUI_LOG` to a per-user folder; on Windows the Start Menu shortcut starts in `Program Files\CoolSolve`,
   not writable, so the log is silently not written (harmless; a per-user default such as
   `%LOCALAPPDATA%\CoolSolve` would be cleaner).
4. **A second launch** starts a second server on the next port (8551, …) instead of reusing the first one.

---

## 7. Native window: options evaluation (not implemented)

The evaluation below (written before v0.4) compares the ways of replacing *open the system browser* by an
application window. Its recommendation, webview/webview, would also solve the `.eescode` association on macOS
(§4.3) and the stopping of the server (§6, item 2): the process ends with its window.

### Goal

Distribute CoolSolve as a **single executable** that:
- Feels like a native desktop application (window title, taskbar icon, no browser address bar)
- Requires **no external programs or libraries** (no Python, no Node, no browser dependency)
- Targets **Windows** (primary), with Linux and macOS as secondary targets
- Can be distributed as a single `.exe` (Windows), `.AppImage` or binary (Linux), `.app` bundle (macOS)

### Current Architecture

CoolSolve already embeds a complete GUI stack:

| Component | Role | Size |
|-----------|------|------|
| **C++ solver** (parser, AD, CoolProp…) | Core engine | ~13 MB binary (with CoolProp); ~37 MB with the embedded model library (v0.4) |
| **cpp-httplib** | Embedded HTTP server (header-only, MIT) | ~100 KB source |
| **React + Monaco + Plotly SPA** | Frontend (TypeScript, compiled to JS/CSS) | ~1.7 MB (dist/) |
| **Embedded assets** | GUI files compiled into the binary | Part of the binary |

The `--gui` flag starts an HTTP server on localhost and opens the system browser. All of the complexity of this evaluation revolves around replacing "open system browser" with "open an embedded native window containing a web view."

### Option 1: webview/webview (C/C++ library) — **Recommended**

**What it is**: A tiny (~300 lines of code per platform) C/C++ library that provides a cross-platform API to open a native web view window. It uses:
- **Windows**: Edge WebView2 (Chromium-based, ships with Windows 10 1803+ and all Windows 11)
- **macOS**: WKWebView (ships with macOS, zero download)
- **Linux**: WebKitGTK (requires `libwebkit2gtk-4.0`, typically pre-installed on GNOME/Ubuntu)

**How it would integrate**:
```
CoolSolve binary (single .exe)
├── C++ solver engine
├── cpp-httplib server (localhost:PORT)
├── Embedded React GUI assets
└── webview window → navigates to http://localhost:PORT
```

The `--gui` codepath would change from `openBrowser("http://localhost:PORT")` to `webview_create() + webview_navigate("http://localhost:PORT")`. The HTTP server runs in a background thread; the webview event loop runs on the main thread.

#### Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Bundle size** | ★★★★★ | Adds ~50 KB to binary. Total: ~13 MB (unchanged). No Chromium bundled. |
| **No external deps (Windows)** | ★★★★★ | WebView2 runtime ships with Windows 10/11. For Windows 7/8 (rare), Microsoft provides a free bootstrapper. |
| **No external deps (macOS)** | ★★★★★ | WKWebView is part of the OS. Zero additional dependencies. |
| **No external deps (Linux)** | ★★★☆☆ | Needs `libwebkit2gtk-4.0-dev`. Pre-installed on Ubuntu/GNOME, but not on minimal server installs. Can document as a prerequisite or statically link with effort. |
| **Native feel** | ★★★★☆ | Real OS window with title bar, taskbar icon, proper Alt+F4/Cmd+Q. No browser address bar. Cannot customize the window chrome (no custom title bar). |
| **Codebase changes** | ★★★★★ | ~50–100 lines of new C++ code. Replace `openBrowser()` with webview calls. Add `webview.h` via FetchContent. No changes to React GUI or solver. |
| **Robustness** | ★★★★☆ | Mature library (5+ years, 13k+ GitHub stars). Edge cases: WebView2 missing on old Windows builds (can detect and fall back to browser). |
| **Maintainability** | ★★★★★ | Header-only C library. No build system changes beyond adding one FetchContent. No new language (stays C++). |
| **Cross-platform build** | ★★★★☆ | CMake native. Windows: MSVC or MinGW. macOS: Xcode. Linux: GCC. CI matrix for 3 platforms easily scriptable. |
| **Single-file installer (Windows)** | ★★★★★ | The `.exe` is standalone. Can also wrap with NSIS or WiX for a proper installer with Start Menu shortcut. |

#### Key implementation sketch

```cpp
#include "webview.h"

// In main.cpp, replace the openBrowser() call:
if (guiMode) {
    std::thread serverThread([&]() {
        coolsolve::startServer(opts);  // blocks on httplib listen loop
    });
    serverThread.detach();

    webview::webview w(false, nullptr);
    w.set_title("CoolSolve");
    w.set_size(1280, 800, WEBVIEW_HINT_NONE);
    w.navigate("http://localhost:" + std::to_string(opts.port));
    w.run();  // blocks on native event loop until window closed
    // When window closes, process exits (server thread dies with it)
}
```

#### Risks
- **WebView2 on Windows 7/8**: Would need the Evergreen WebView2 Runtime bootstrapper (~1.5 MB download). Negligible risk since Windows 7/8 market share is < 1%.
- **Linux WebKitGTK version**: Different distros ship different WebKitGTK versions. Target `webkit2gtk-4.0` (Ubuntu 20.04+) or `webkit2gtk-4.1` (Ubuntu 22.04+).

---

### Option 2: Electron

**What it is**: A framework that bundles Chromium + Node.js to run web apps as desktop applications. Very popular (VS Code, Slack, Discord).

**How it would integrate**: Wrap the CoolSolve binary as a child process. Electron spawns it, communicates over localhost HTTP, and displays the React GUI in its own Chromium instance.

#### Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Bundle size** | ★☆☆☆☆ | Adds ~150–200 MB (full Chromium + Node.js). Total: **~170+ MB**. This is the main dealbreaker. |
| **No external deps** | ★★★★★ | Self-contained (bundles its own Chromium). |
| **Native feel** | ★★★★★ | Full control over window chrome, system tray, menus, notifications, auto-update. Best DX for desktop features. |
| **Codebase changes** | ★★☆☆☆ | Need an `electron/` project directory with `main.js`, `preload.js`, packaging config. Need to handle spawning the C++ binary as a child process and managing its lifecycle. Different build pipeline. |
| **Robustness** | ★★★★★ | Battle-tested at scale (VS Code, Slack). Huge community. |
| **Maintainability** | ★★☆☆☆ | Adds Node.js/npm to the build chain. Electron major version upgrades (~quarterly) can introduce breaking changes. Two separate build systems to maintain. |
| **Cross-platform build** | ★★★★☆ | `electron-builder` handles Windows (NSIS), macOS (DMG), Linux (AppImage). But CI must build the C++ binary per-platform AND the Electron wrapper. |
| **Single-file installer** | ★★★★☆ | `electron-builder` can produce single `.exe` installers. But they're ~170 MB. |

#### Risks
- **Doubling the tech stack**: CoolSolve is a C++ project. Electron adds Node.js, npm, JavaScript packaging — a separate ecosystem to maintain.
- **Update burden**: Chromium security patches ship ~every 6 weeks. If not updated, the bundled Chromium becomes a liability.
- **Overkill**: CoolSolve already has an embedded HTTP server and compiled SPA. Electron's value-add is just the window frame.

---

### Option 3: Tauri

**What it is**: A Rust-based framework that uses the OS native web view (same as webview/webview under the hood) but adds an application framework: auto-update, system tray, file dialogs, signing, etc.

**How it would integrate**: The Tauri Rust process would spawn CoolSolve as a sidecar binary, connect its webview to `localhost:PORT`, and provide native integration (file dialogs, menus).

#### Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Bundle size** | ★★★★☆ | Adds ~5–10 MB for the Tauri runtime. Total: ~20–25 MB. Uses OS webview (same as Option 1). |
| **No external deps** | ★★★★☆ | Same webview story as Option 1. Additionally needs Rust toolchain for building (not for end users). |
| **Native feel** | ★★★★★ | System tray, custom menus, file dialogs, auto-update, code signing — all built-in. Best native integration of all options. |
| **Codebase changes** | ★★☆☆☆ | Need a `src-tauri/` directory with Rust config, sidecar binary configuration, Tauri config file. Must configure sidecar for each platform (different binary names). |
| **Robustness** | ★★★★☆ | Tauri v2 is stable. Large and active community. But younger than Electron. |
| **Maintainability** | ★★☆☆☆ | **Adds Rust as a second language** to the project. Tauri config format changes between major versions. Requires Rust toolchain on CI. |
| **Cross-platform build** | ★★★☆☆ | `tauri-cli` handles per-platform packaging. But cross-compilation is harder than native builds. CI must build C++ sidecar + Rust Tauri wrapper per-platform. |
| **Single-file installer** | ★★★★★ | Produces `.msi`/`.exe` (Windows), `.dmg` (macOS), `.AppImage`/`.deb` (Linux). Excellent packaging. |

#### Risks
- **Rust dependency**: CoolSolve is pure C++. Adding Rust introduces a language and toolchain that must be maintained by the team. For a thin wrapper, this is high overhead.
- **Sidecar complexity**: Tauri's sidecar model requires signing the binary on macOS, and the sidecar must be configured per-architecture (x86_64, aarch64).
- **Overengineered**: For what is essentially "open a window pointing at localhost", the Tauri framework adds significant ceremony.

---

### Option 4: Qt WebEngine

**What it is**: Qt's built-in Chromium-based web view. Qt is a mature C++ GUI framework.

#### Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Bundle size** | ★★☆☆☆ | Qt WebEngine bundles Chromium: **~80–120 MB** on top of the base binary. |
| **No external deps** | ★★★★★ | Self-contained (ships its own Chromium). |
| **Native feel** | ★★★★★ | Full native menus, dialogs, system tray, multi-window. Qt is the gold standard for C++ desktop apps. |
| **Codebase changes** | ★★★☆☆ | Stay in C++, but need Qt's build system (qmake or CMake with Qt). The Qt `QWebEngineView` widget would host the SPA. ~100 lines of Qt code. |
| **Robustness** | ★★★★★ | Qt is decades old, rock-solid, commercially supported. |
| **Maintainability** | ★★★☆☆ | Qt is C++ (no new language), but it's a large dependency. Qt's licensing (LGPL/commercial) may constrain distribution if not careful. |
| **Cross-platform build** | ★★★☆☆ | Qt supports all 3 platforms, but the WebEngine module is heavy to build and link. Static linking is non-trivial. |
| **Single-file installer** | ★★★☆☆ | On Windows, must either bundle Qt DLLs (use `windeployqt`) or statically link (complex). Not a single `.exe` without extra work. |

#### Risks
- **Licensing**: Qt WebEngine is under LGPL. Must dynamically link Qt or obtain a commercial license if distributing a proprietary application.
- **Build complexity**: Qt WebEngine takes 30+ minutes to build from source. Using pre-built binaries adds ~500 MB to the development setup.
- **Size**: 80–120 MB is excessive for a window wrapper.

---

### Option 5: Status Quo (Browser-based, enhanced)

**What it is**: Keep the current architecture but improve the "launch browser" experience: auto-detect Chrome/Edge and launch in `--app` mode (kiosk-like, no address bar), with a system tray icon.

#### Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Bundle size** | ★★★★★ | Zero additional bytes. 13 MB binary as-is. |
| **No external deps** | ★★☆☆☆ | **Requires a browser**. Chrome/Edge `--app` mode gives near-native feel, but the user needs Chrome or Edge installed. |
| **Native feel** | ★★★☆☆ | `--app` mode removes address bar and tabs. Still shows a browser icon in the taskbar (Chrome icon, not CoolSolve icon). No custom window title without PWA manifest. |
| **Codebase changes** | ★★★★★ | ~20 lines: detect Chrome/Edge path, launch with `--app=http://localhost:PORT --window-size=1280,800`. |
| **Robustness** | ★★★☆☆ | Works only if Chrome or Edge is installed. Firefox has no equivalent `--app` mode. |
| **Maintainability** | ★★★★★ | Trivial. No new dependencies. |
| **Single-file installer** | ★★★★★ | Already a single binary. |

#### Enhancement details
```bash
# Chrome --app mode (near-native, no address bar):
chrome --app=http://localhost:8550 --window-size=1280,800

# Edge --app mode (Windows, almost always available):
msedge --app=http://localhost:8550 --window-size=1280,800
```

Can add a **PWA manifest** (`manifest.json`) so that when the user "installs" the app from the browser, it gets a proper icon and title bar. This is zero-cost since the manifest is just a JSON file embedded with the other GUI assets.

---

### Comparison Matrix

| Criterion | webview/webview | Electron | Tauri | Qt WebEngine | Browser --app |
|-----------|:-:|:-:|:-:|:-:|:-:|
| **Bundle size overhead** | ~0 MB | +150 MB | +10 MB | +100 MB | 0 MB |
| **Total estimated size** | ~13 MB | ~170 MB | ~25 MB | ~120 MB | ~13 MB |
| **No external deps (Win)** | ✅ | ✅ | ✅ | ✅ | ❌ (needs browser) |
| **No external deps (macOS)** | ✅ | ✅ | ✅ | ✅ | ❌ (needs browser) |
| **No external deps (Linux)** | ⚠️ WebKitGTK | ✅ | ⚠️ WebKitGTK | ✅ | ❌ (needs browser) |
| **Native window** | ✅ | ✅ | ✅ | ✅ | ⚠️ (browser --app) |
| **Custom icon/title** | ✅ | ✅ | ✅ | ✅ | ⚠️ (via PWA) |
| **System tray** | ❌ (add manually) | ✅ built-in | ✅ built-in | ✅ built-in | ❌ |
| **Auto-update** | ❌ (add manually) | ✅ built-in | ✅ built-in | ❌ (add manually) | ❌ |
| **File dialogs** | ❌ (use HTML5) | ✅ native | ✅ native | ✅ native | ✅ (HTML5) |
| **New language required** | No (C/C++) | JavaScript | **Rust** | No (C++) | No |
| **Build complexity delta** | Minimal | High | High | High | None |
| **Lines of new code** | ~50–100 | ~300+ | ~200+ | ~100–200 | ~20 |
| **Single .exe (Windows)** | ✅ | ✅ (large) | ✅ | ⚠️ (DLLs) | ✅ (already) |

---

### Recommendation

#### Primary recommendation: **webview/webview** (Option 1)

This is the best fit for CoolSolve because:

1. **Minimal disruption**: CoolSolve already does 95% of what a desktop app needs (embedded server, compiled SPA, single binary). The only missing piece is "open a window instead of a browser tab." webview/webview provides exactly that and nothing more.

2. **Single binary stays single**: The ~13 MB executable remains ~13 MB. No Chromium, no Node.js, no Rust runtime. This aligns perfectly with the "single executable" requirement.

3. **No new languages or toolchains**: Stays pure C++/CMake. No Rust, no npm for the wrapper, no separate build pipeline.

4. **Windows coverage is excellent**: WebView2 ships with every modern Windows installation. Edge (Chromium) has been the default browser since Windows 10 20H2 (October 2020).

5. **macOS works out of the box**: WKWebView is part of the OS since macOS 10.10.

6. **The fallback is free**: If webview creation fails (e.g., very old Windows without WebView2), fall back to `openBrowser()` — the current behavior. Users lose nothing.

#### Implementation plan

1. **Add webview/webview via FetchContent** in CMakeLists.txt (header-only, ~50 KB)
2. **Modify `main.cpp`**: When `--gui` is passed, start the HTTP server in a background thread, then create a webview window pointing to `localhost:PORT`
3. **Add a platform-specific icon**: `.ico` for Windows, `.icns` for macOS — set via webview API or resource file
4. **Windows installer**: Optionally wrap with NSIS or WiX to create a Start Menu shortcut and file association for `.eescode` files
5. **Keep `--gui --no-browser` mode**: For headless/server deployments, the current behavior is preserved unchanged

#### Secondary recommendation: **Browser --app** (Option 5) as an interim step

If immediate results are needed before implementing webview, enhance `openBrowser()` to detect Chrome/Edge and launch in `--app` mode. This gives 80% of the native feel with ~20 lines of code and zero new dependencies. Add a PWA `manifest.json` for proper icon and title. This can ship today and coexist with the webview approach later.
