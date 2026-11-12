#pragma once

#include <string>

namespace coolsolve {

/** Versions of the running CoolSolve binary and of the CoolProp it is linked with. */
struct VersionInfo {
    std::string coolsolve;          ///< CoolSolve version (from CMakeLists.txt)
    std::string coolprop;           ///< CoolProp version string, e.g. "8.1.0dev"
    std::string coolpropRevision;   ///< CoolProp git revision (full SHA-1), empty if unknown
    std::string library;            ///< embedded model library ("<commit>, snapshot <hash>"), empty if none
    std::string coolplotCore;       ///< coolplot-core version (coolplot_version()), empty if unknown
    std::string coolplotCommit;     ///< CoolPlot commit coolplot-core is built from (COOLSOLVE_COOLPLOT_TAG)
};

/**
 * Read the versions. The CoolProp values are queried at run time from the linked
 * library (`get_global_param_string("version")` / `("gitrevision")`), so they
 * describe the property library that produced the results, not what CMake asked for.
 */
VersionInfo getVersionInfo();

/**
 * Text printed by `coolsolve --version` (newline-terminated lines):
 *
 *     CoolSolve 0.3.0
 *     CoolProp 8.1.0dev (git revision <40 hex digits>)
 *     CoolSolve Library <commit> (snapshot <hash>)        <- only when the library is embedded
 *     coolplot-core 1.5.0 (CoolPlot <commit>)             <- the diagrams of the GUI
 *
 * Tools (e.g. the model library's `tools/test_models.py`) parse the second line.
 */
std::string versionText(const VersionInfo& info);

}  // namespace coolsolve
