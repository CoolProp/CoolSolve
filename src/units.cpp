#include "coolsolve/units.h"
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace coolsolve {

// Helper to normalize unit strings (lowercase, remove spaces)
static std::string normalize(const std::string& unit) {
    std::string s = unit;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    // Remove spaces
    s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
    return s;
}

double UnitConverter::toSI(double value, UnitType type, const std::string& unit) {
    std::string u = normalize(unit);
    
    switch (type) {
        case UnitType::Temperature:
            if (u == "c" || u == "celsius") return value + 273.15;
            if (u == "k" || u == "kelvin") return value;
            if (u == "f" || u == "fahrenheit") return (value - 32.0) * 5.0/9.0 + 273.15;
            if (u == "r" || u == "rankine") return value * 5.0/9.0;
            break;
            
        case UnitType::Pressure:
            if (u == "pa" || u == "pascal") return value;
            if (u == "kpa") return value * 1000.0;
            if (u == "mpa") return value * 1e6;
            if (u == "bar") return value * 1e5;
            if (u == "atm") return value * 101325.0;
            if (u == "psi" || u == "psia") return value * 6894.75729;
            break;
            
        case UnitType::Energy:
            if (u == "j" || u == "joule") return value;
            if (u == "kj") return value * 1000.0;
            if (u == "mj") return value * 1e6;
            if (u == "cal") return value * 4.184;
            if (u == "kcal") return value * 4184.0;
            if (u == "btu") return value * 1055.05585;
            break;
            
        case UnitType::Mass:
            if (u == "kg" || u == "kilogram") return value;
            if (u == "g" || u == "gram") return value / 1000.0;
            if (u == "lb" || u == "lbm") return value * 0.45359237;
            break;
            
        case UnitType::Length:
            if (u == "m" || u == "meter") return value;
            if (u == "cm") return value / 100.0;
            if (u == "mm") return value / 1000.0;
            if (u == "ft") return value * 0.3048;
            if (u == "in" || u == "inch") return value * 0.0254;
            break;
            
        case UnitType::Time:
            if (u == "s" || u == "sec" || u == "second") return value;
            if (u == "min" || u == "minute") return value * 60.0;
            if (u == "h" || u == "hr" || u == "hour") return value * 3600.0;
            break;
            
        case UnitType::Power:
            if (u == "w" || u == "watt") return value;
            if (u == "kw") return value * 1000.0;
            if (u == "mw") return value * 1e6;
            break;
            
        case UnitType::SpecificHeat:
        case UnitType::SpecificEntropy:
            if (u == "j/kg-k" || u == "j/kgk") return value;
            if (u == "kj/kg-k" || u == "kj/kgk") return value * 1000.0;
            break;
            
        case UnitType::SpecificEnergy:
            if (u == "j/kg") return value;
            if (u == "kj/kg") return value * 1000.0;
            break;
            
        case UnitType::Conductivity:
            if (u == "w/m-k" || u == "w/mk") return value;
            break;
            
        case UnitType::Viscosity:
            if (u == "pa-s" || u == "pas") return value;
            if (u == "poise") return value * 0.1;
            if (u == "cp" || u == "centipoise") return value * 0.001;
            break;
            
        case UnitType::Density:
            if (u == "kg/m^3" || u == "kg/m3") return value;
            break;
            
        case UnitType::Dimensionless:
            return value;
    }
    
    // Default: assume no conversion if unit not recognized
    // This allows passing through values if units match implicitly
    return value;
}

double UnitConverter::fromSI(double value, UnitType type, const std::string& unit) {
    std::string u = normalize(unit);
    
    switch (type) {
        case UnitType::Temperature:
            if (u == "c" || u == "celsius") return value - 273.15;
            if (u == "k" || u == "kelvin") return value;
            if (u == "f" || u == "fahrenheit") return (value - 273.15) * 9.0/5.0 + 32.0;
            if (u == "r" || u == "rankine") return value * 9.0/5.0;
            break;
            
        case UnitType::Pressure:
            if (u == "pa" || u == "pascal") return value;
            if (u == "kpa") return value / 1000.0;
            if (u == "mpa") return value / 1e6;
            if (u == "bar") return value / 1e5;
            if (u == "atm") return value / 101325.0;
            if (u == "psi" || u == "psia") return value / 6894.75729;
            break;
            
        case UnitType::Energy:
            if (u == "j" || u == "joule") return value;
            if (u == "kj") return value / 1000.0;
            if (u == "mj") return value / 1e6;
            if (u == "cal") return value / 4.184;
            if (u == "kcal") return value / 4184.0;
            if (u == "btu") return value / 1055.05585;
            break;
            
        case UnitType::Mass:
            if (u == "kg" || u == "kilogram") return value;
            if (u == "g" || u == "gram") return value * 1000.0;
            if (u == "lb" || u == "lbm") return value / 0.45359237;
            break;
            
        case UnitType::Length:
            if (u == "m" || u == "meter") return value;
            if (u == "cm") return value * 100.0;
            if (u == "mm") return value * 1000.0;
            if (u == "ft") return value / 0.3048;
            if (u == "in" || u == "inch") return value / 0.0254;
            break;
            
        case UnitType::Time:
            if (u == "s" || u == "sec" || u == "second") return value;
            if (u == "min" || u == "minute") return value / 60.0;
            if (u == "h" || u == "hr" || u == "hour") return value / 3600.0;
            break;
            
        case UnitType::Power:
            if (u == "w" || u == "watt") return value;
            if (u == "kw") return value / 1000.0;
            if (u == "mw") return value / 1e6;
            break;
            
        case UnitType::SpecificHeat:
        case UnitType::SpecificEntropy:
            if (u == "j/kg-k" || u == "j/kgk") return value;
            if (u == "kj/kg-k" || u == "kj/kgk") return value / 1000.0;
            break;
            
        case UnitType::SpecificEnergy:
            if (u == "j/kg") return value;
            if (u == "kj/kg") return value / 1000.0;
            break;
            
        case UnitType::Conductivity:
            if (u == "w/m-k" || u == "w/mk") return value;
            break;
            
        case UnitType::Viscosity:
            if (u == "pa-s" || u == "pas") return value;
            if (u == "poise") return value * 10.0;
            if (u == "cp" || u == "centipoise") return value * 1000.0;
            break;
            
        case UnitType::Density:
            if (u == "kg/m^3" || u == "kg/m3") return value;
            break;
            
        case UnitType::Dimensionless:
            return value;
    }
    
    return value;
}

} // namespace coolsolve
