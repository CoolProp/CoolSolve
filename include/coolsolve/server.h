#pragma once

#ifdef COOLSOLVE_GUI

#include <string>
#include <functional>
#include <vector>

namespace coolsolve {

/// Options for the embedded HTTP server
struct ServerOptions {
    int port = 8550;                    // Default port
    std::string host = "127.0.0.1";     // Address to listen on. The default accepts connections
                                        // from this computer only; "0.0.0.0" listens on every
                                        // network interface (the API has no authentication)
    std::vector<std::string> allowedHosts;  // Host names the server answers to, besides
                                            // "localhost" and IP addresses (DNS rebinding guard)
    bool openBrowser = true;            // Open default browser on start
    bool onlineMode = false;            // Enable session sandboxing (future)
    std::string guiDevDir;              // If set, serve GUI from filesystem (dev mode)
    std::string docsDir;                // If set, serve docs from filesystem (dev/fallback)
    std::string initialFile;            // If set, auto-open this file in the GUI on first load
    std::string usageLogFile;           // Solve-attempt log path (empty = COOLSOLVE_GUI_LOG env
                                        // var, else "coolsolve_gui.log" in the working directory)
};

/// Start the embedded HTTP server (blocking).
/// Returns the exit code (0 = normal shutdown).
int startServer(const ServerOptions& options);

} // namespace coolsolve

#endif // COOLSOLVE_GUI
