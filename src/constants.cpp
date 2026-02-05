#include "coolsolve/constants.h"
#include <cmath>

namespace coolsolve {

// Define the registry of constants
const std::map<std::string, ConstantInfo> Constants::registry_ = {
    // PI variants (special case, no # suffix)
    {"PI", {3.14159265358979323846, "", "Pi"}},
    {"pi", {3.14159265358979323846, "", "Pi"}},
    {"Pi", {3.14159265358979323846, "", "Pi"}},
    
    // Physical constants (ending in #)
    {"g#", {9.80665, "m/s^2", "Standard gravity"}},
    {"R#", {8314.462618, "J/kmol-K", "Universal gas constant"}},
    {"k#", {1.380649e-23, "J/K", "Boltzmann constant"}},
    {"sigma#", {5.670374419e-8, "W/m^2-K^4", "Stefan-Boltzmann constant"}},
    {"Po#", {101325.0, "Pa", "Standard atmospheric pressure"}},
    {"To#", {298.15, "K", "Standard temperature (25 C)"}},
    {"c#", {299792458.0, "m/s", "Speed of light in vacuum"}},
    {"h#", {6.62607015e-34, "J-s", "Planck constant"}},
    {"e#", {1.602176634e-19, "C", "Elementary charge"}},
    {"me#", {9.10938356e-31, "kg", "Electron mass"}},
    {"mp#", {1.6726219e-27, "kg", "Proton mass"}},
    {"mn#", {1.674927471e-27, "kg", "Neutron mass"}},
    {"NA#", {6.02214076e23, "1/mol", "Avogadro constant"}},
    {"Avogadro#", {6.02214076e23, "1/mol", "Avogadro constant"}}
};

bool Constants::isConstant(const std::string& name) {
    return registry_.find(name) != registry_.end();
}

double Constants::getValue(const std::string& name) {
    auto it = registry_.find(name);
    if (it != registry_.end()) {
        return it->second.value;
    }
    return 0.0;
}

std::optional<ConstantInfo> Constants::getInfo(const std::string& name) {
    auto it = registry_.find(name);
    if (it != registry_.end()) {
        return it->second;
    }
    return std::nullopt;
}

const std::map<std::string, ConstantInfo>& Constants::getAll() {
    return registry_;
}

} // namespace coolsolve
