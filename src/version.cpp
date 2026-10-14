#include "coolsolve/version.h"
#include "CoolProp.h"
#ifdef COOLSOLVE_EMBEDDED_LIBRARY
#include "embedded_library.h"
#endif

namespace coolsolve {

VersionInfo getVersionInfo() {
    VersionInfo info;
    info.coolsolve = COOLSOLVE_VERSION;
    info.coolprop = CoolProp::get_global_param_string("version");
    info.coolpropRevision = CoolProp::get_global_param_string("gitrevision");
#ifdef COOLSOLVE_EMBEDDED_LIBRARY
    info.library = std::string(embeddedLibraryCommit()) + ", snapshot " + embeddedLibraryHash();
#endif
    return info;
}

std::string versionText(const VersionInfo& info) {
    std::string text = "CoolSolve " + info.coolsolve + "\n";
    text += "CoolProp " + info.coolprop;
    text += " (git revision " + (info.coolpropRevision.empty() ? std::string("unknown") : info.coolpropRevision) + ")\n";
    if (!info.library.empty()) text += "CoolSolve Library " + info.library + "\n";
    return text;
}

}  // namespace coolsolve
