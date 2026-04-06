# CoolSolve Standalone Desktop Application — Options Evaluation

## Goal

Distribute CoolSolve as a **single executable** that:
- Feels like a native desktop application (window title, taskbar icon, no browser address bar)
- Requires **no external programs or libraries** (no Python, no Node, no browser dependency)
- Targets **Windows** (primary), with Linux and macOS as secondary targets
- Can be distributed as a single `.exe` (Windows), `.AppImage` or binary (Linux), `.app` bundle (macOS)

## Current Architecture

CoolSolve already embeds a complete GUI stack:

| Component | Role | Size |
|-----------|------|------|
| **C++ solver** (parser, AD, CoolProp…) | Core engine | ~13 MB binary (with CoolProp) |
| **cpp-httplib** | Embedded HTTP server (header-only, MIT) | ~100 KB source |
| **React + Monaco + Plotly SPA** | Frontend (TypeScript, compiled to JS/CSS) | ~1.7 MB (dist/) |
| **Embedded assets** | GUI files compiled into the binary | Part of the 13 MB |

The `--gui` flag starts an HTTP server on localhost and opens the system browser. All of the complexity of this evaluation revolves around replacing "open system browser" with "open an embedded native window containing a web view."

---

## Option 1: webview/webview (C/C++ library) — **Recommended**

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

### Assessment

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

### Key implementation sketch

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

### Risks
- **WebView2 on Windows 7/8**: Would need the Evergreen WebView2 Runtime bootstrapper (~1.5 MB download). Negligible risk since Windows 7/8 market share is < 1%.
- **Linux WebKitGTK version**: Different distros ship different WebKitGTK versions. Target `webkit2gtk-4.0` (Ubuntu 20.04+) or `webkit2gtk-4.1` (Ubuntu 22.04+).

---

## Option 2: Electron

**What it is**: A framework that bundles Chromium + Node.js to run web apps as desktop applications. Very popular (VS Code, Slack, Discord).

**How it would integrate**: Wrap the CoolSolve binary as a child process. Electron spawns it, communicates over localhost HTTP, and displays the React GUI in its own Chromium instance.

### Assessment

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

### Risks
- **Doubling the tech stack**: CoolSolve is a C++ project. Electron adds Node.js, npm, JavaScript packaging — a separate ecosystem to maintain.
- **Update burden**: Chromium security patches ship ~every 6 weeks. If not updated, the bundled Chromium becomes a liability.
- **Overkill**: CoolSolve already has an embedded HTTP server and compiled SPA. Electron's value-add is just the window frame.

---

## Option 3: Tauri

**What it is**: A Rust-based framework that uses the OS native web view (same as webview/webview under the hood) but adds an application framework: auto-update, system tray, file dialogs, signing, etc.

**How it would integrate**: The Tauri Rust process would spawn CoolSolve as a sidecar binary, connect its webview to `localhost:PORT`, and provide native integration (file dialogs, menus).

### Assessment

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

### Risks
- **Rust dependency**: CoolSolve is pure C++. Adding Rust introduces a language and toolchain that must be maintained by the team. For a thin wrapper, this is high overhead.
- **Sidecar complexity**: Tauri's sidecar model requires signing the binary on macOS, and the sidecar must be configured per-architecture (x86_64, aarch64).
- **Overengineered**: For what is essentially "open a window pointing at localhost", the Tauri framework adds significant ceremony.

---

## Option 4: Qt WebEngine

**What it is**: Qt's built-in Chromium-based web view. Qt is a mature C++ GUI framework.

### Assessment

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

### Risks
- **Licensing**: Qt WebEngine is under LGPL. Must dynamically link Qt or obtain a commercial license if distributing a proprietary application.
- **Build complexity**: Qt WebEngine takes 30+ minutes to build from source. Using pre-built binaries adds ~500 MB to the development setup.
- **Size**: 80–120 MB is excessive for a window wrapper.

---

## Option 5: Status Quo (Browser-based, enhanced)

**What it is**: Keep the current architecture but improve the "launch browser" experience: auto-detect Chrome/Edge and launch in `--app` mode (kiosk-like, no address bar), with a system tray icon.

### Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Bundle size** | ★★★★★ | Zero additional bytes. 13 MB binary as-is. |
| **No external deps** | ★★☆☆☆ | **Requires a browser**. Chrome/Edge `--app` mode gives near-native feel, but the user needs Chrome or Edge installed. |
| **Native feel** | ★★★☆☆ | `--app` mode removes address bar and tabs. Still shows a browser icon in the taskbar (Chrome icon, not CoolSolve icon). No custom window title without PWA manifest. |
| **Codebase changes** | ★★★★★ | ~20 lines: detect Chrome/Edge path, launch with `--app=http://localhost:PORT --window-size=1280,800`. |
| **Robustness** | ★★★☆☆ | Works only if Chrome or Edge is installed. Firefox has no equivalent `--app` mode. |
| **Maintainability** | ★★★★★ | Trivial. No new dependencies. |
| **Single-file installer** | ★★★★★ | Already a single binary. |

### Enhancement details
```bash
# Chrome --app mode (near-native, no address bar):
chrome --app=http://localhost:8550 --window-size=1280,800

# Edge --app mode (Windows, almost always available):
msedge --app=http://localhost:8550 --window-size=1280,800
```

Can add a **PWA manifest** (`manifest.json`) so that when the user "installs" the app from the browser, it gets a proper icon and title bar. This is zero-cost since the manifest is just a JSON file embedded with the other GUI assets.

---

## Comparison Matrix

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

## Recommendation

### Primary recommendation: **webview/webview** (Option 1)

This is the best fit for CoolSolve because:

1. **Minimal disruption**: CoolSolve already does 95% of what a desktop app needs (embedded server, compiled SPA, single binary). The only missing piece is "open a window instead of a browser tab." webview/webview provides exactly that and nothing more.

2. **Single binary stays single**: The ~13 MB executable remains ~13 MB. No Chromium, no Node.js, no Rust runtime. This aligns perfectly with the "single executable" requirement.

3. **No new languages or toolchains**: Stays pure C++/CMake. No Rust, no npm for the wrapper, no separate build pipeline.

4. **Windows coverage is excellent**: WebView2 ships with every modern Windows installation. Edge (Chromium) has been the default browser since Windows 10 20H2 (October 2020).

5. **macOS works out of the box**: WKWebView is part of the OS since macOS 10.10.

6. **The fallback is free**: If webview creation fails (e.g., very old Windows without WebView2), fall back to `openBrowser()` — the current behavior. Users lose nothing.

### Implementation plan

1. **Add webview/webview via FetchContent** in CMakeLists.txt (header-only, ~50 KB)
2. **Modify `main.cpp`**: When `--gui` is passed, start the HTTP server in a background thread, then create a webview window pointing to `localhost:PORT`
3. **Add a platform-specific icon**: `.ico` for Windows, `.icns` for macOS — set via webview API or resource file
4. **Windows installer**: Optionally wrap with NSIS or WiX to create a Start Menu shortcut and file association for `.eescode` files
5. **Keep `--gui --no-browser` mode**: For headless/server deployments, the current behavior is preserved unchanged

### Secondary recommendation: **Browser --app** (Option 5) as an interim step

If immediate results are needed before implementing webview, enhance `openBrowser()` to detect Chrome/Edge and launch in `--app` mode. This gives 80% of the native feel with ~20 lines of code and zero new dependencies. Add a PWA `manifest.json` for proper icon and title. This can ship today and coexist with the webview approach later.
