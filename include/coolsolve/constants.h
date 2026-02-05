#pragma once

#include <string>
#include <map>
#include <optional>

namespace coolsolve {

struct ConstantInfo {
    double value;
    std::string units;
    std::string description;
};

class Constants {
public:
    // Check if a name corresponds to a known constant
    static bool isConstant(const std::string& name);
    
    // Get the value of a constant (returns 0.0 if not found)
    static double getValue(const std::string& name);
    
    // Get full info for a constant
    static std::optional<ConstantInfo> getInfo(const std::string& name);
    
    // Get the full map of constants
    static const std::map<std::string, ConstantInfo>& getAll();

private:
    static const std::map<std::string, ConstantInfo> registry_;
};

} // namespace coolsolve
