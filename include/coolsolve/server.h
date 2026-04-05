#pragma once

#ifdef COOLSOLVE_GUI

#include <string>
#include <functional>

namespace coolsolve {

/// Options for the embedded HTTP server
struct ServerOptions {
    int port = 8550;                    // Default port
    bool openBrowser = true;            // Open default browser on start
    bool onlineMode = false;            // Enable session sandboxing (future)
    std::string guiDevDir;              // If set, serve GUI from filesystem (dev mode)
    std::string docsDir;                // If set, serve docs from filesystem (dev/fallback)
};

/// Start the embedded HTTP server (blocking).
/// Returns the exit code (0 = normal shutdown).
int startServer(const ServerOptions& options);

} // namespace coolsolve

#endif // COOLSOLVE_GUI
