#include "coolsolve/version.h"
#include "CoolProp.h"

namespace coolsolve {

VersionInfo getVersionInfo() {
    VersionInfo info;
    info.coolsolve = COOLSOLVE_VERSION;
    info.coolprop = CoolProp::get_global_param_string("version");
    info.coolpropRevision = CoolProp::get_global_param_string("gitrevision");
    return info;
}

std::string versionText(const VersionInfo& info) {
    std::string text = "CoolSolve " + info.coolsolve + "\n";
    text += "CoolProp " + info.coolprop;
    text += " (git revision " + (info.coolpropRevision.empty() ? std::string("unknown") : info.coolpropRevision) + ")\n";
    return text;
}

}  // namespace coolsolve
