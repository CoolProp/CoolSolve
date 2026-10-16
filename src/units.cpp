#include "coolsolve/units.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <stdexcept>
#include <unordered_map>
#include <vector>
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

namespace {

/// Exponents of the base dimensions of a unit: mass, length, time, temperature (a temperature difference),
/// amount of substance and electric current.
using Dimension = std::array<int, 6>;

/// A unit as a product of SI base units: 1 [unit] = factor * kg^e0 * m^e1 * s^e2 * K^e3 * mol^e4 * A^e5.
struct UnitValue {
    double factor = 1.0;
    Dimension dim{};
};

constexpr double kPi = 3.14159265358979323846;

constexpr Dimension kNone{0, 0, 0, 0, 0, 0};
constexpr Dimension kMass{1, 0, 0, 0, 0, 0};
constexpr Dimension kLength{0, 1, 0, 0, 0, 0};
constexpr Dimension kTime{0, 0, 1, 0, 0, 0};
constexpr Dimension kTemperature{0, 0, 0, 1, 0, 0};
constexpr Dimension kAmount{0, 0, 0, 0, 1, 0};
constexpr Dimension kCurrent{0, 0, 0, 0, 0, 1};
constexpr Dimension kArea{0, 2, 0, 0, 0, 0};
constexpr Dimension kVolume{0, 3, 0, 0, 0, 0};
constexpr Dimension kVelocity{0, 1, -1, 0, 0, 0};
constexpr Dimension kFrequency{0, 0, -1, 0, 0, 0};
constexpr Dimension kForce{1, 1, -2, 0, 0, 0};
constexpr Dimension kPressure{1, -1, -2, 0, 0, 0};
constexpr Dimension kEnergy{1, 2, -2, 0, 0, 0};
constexpr Dimension kPower{1, 2, -3, 0, 0, 0};
constexpr Dimension kVoltage{1, 2, -3, 0, 0, -1};
constexpr Dimension kViscosity{1, -1, -1, 0, 0, 0};

/// The units that `CONVERT` knows by name, with their SI factor and dimension. The factors of the units that
/// UnitConverter converts are taken from it; the others are written here, in SI.
const std::unordered_map<std::string, UnitValue>& simpleUnits() {
    static const std::unordered_map<std::string, UnitValue> table = [] {
        std::unordered_map<std::string, UnitValue> t;
        const auto converted = [&](UnitType type, const Dimension& dim, std::initializer_list<const char*> names) {
            for (const char* name : names) t[name] = {UnitConverter::toSI(1.0, type, name), dim};
        };
        const auto add = [&](const Dimension& dim, std::initializer_list<std::pair<const char*, double>> units) {
            for (const auto& [name, factor] : units) t[name] = {factor, dim};
        };

        converted(UnitType::Mass, kMass, {"kg", "g", "lb", "lbm", "kilogram", "gram"});
        converted(UnitType::Length, kLength, {"m", "cm", "mm", "ft", "in", "meter", "inch"});
        converted(UnitType::Time, kTime, {"s", "sec", "min", "h", "hr", "second", "minute", "hour"});
        converted(UnitType::Pressure, kPressure, {"Pa", "kPa", "MPa", "bar", "atm", "psi", "psia", "pascal"});
        converted(UnitType::Energy, kEnergy, {"J", "kJ", "MJ", "cal", "kcal", "Btu", "joule"});
        converted(UnitType::Power, kPower, {"W", "kW", "MW", "watt"});
        converted(UnitType::Viscosity, kViscosity, {"poise", "cP", "centipoise"});

        // a temperature unit in a unit is a difference of temperature: K = C, R = F = 5/9 K
        const double rankine = UnitConverter::toSI(1.0, UnitType::Temperature, "R");
        add(kTemperature, {{"K", 1.0}, {"kelvin", 1.0}, {"C", 1.0}, {"celsius", 1.0},
                           {"F", rankine}, {"fahrenheit", rankine}, {"R", rankine}, {"rankine", rankine}});

        add(kMass, {{"mg", 1e-6}, {"oz", 0.028349523125}, {"slug", 14.5939029372064}, {"tonne", 1000.0}});
        add(kLength, {{"km", 1e3}, {"dm", 0.1}, {"um", 1e-6}, {"nm", 1e-9}, {"yd", 0.9144}, {"mi", 1609.344}});
        add(kTime, {{"ms", 1e-3}, {"us", 1e-6}, {"ns", 1e-9}, {"day", 86400.0}, {"week", 604800.0}});
        add(kForce, {{"N", 1.0}, {"kN", 1e3}, {"lbf", 4.4482216152605}, {"kgf", 9.80665}, {"dyne", 1e-5}});
        const double mmHg = 133.322387415;
        add(kPressure, {{"mPa", 1e-3}, {"hPa", 1e2}, {"GPa", 1e9}, {"mbar", 1e2}, {"torr", 101325.0 / 760.0},
                        {"mmHg", mmHg}, {"inHg", 25.4 * mmHg},
                        {"ksi", 1e3 * UnitConverter::toSI(1.0, UnitType::Pressure, "psi")}});
        add(kEnergy, {{"GJ", 1e9}, {"Wh", 3600.0}, {"kWh", 3.6e6}});
        add(kPower, {{"mW", 1e-3}, {"GW", 1e9}, {"hp", 745.699871582270}});
        add(kVolume, {{"L", 1e-3}, {"liter", 1e-3}, {"litre", 1e-3}, {"mL", 1e-6}, {"gal", 3.785411784e-3}, {"qt", 9.46352946e-4}});
        add(kArea, {{"acre", 4046.8564224}, {"ha", 1e4}});
        add(kVelocity, {{"mph", 0.44704}, {"knot", 1852.0 / 3600.0}});
        add(kAmount, {{"mol", 1.0}, {"mmol", 1e-3}, {"kmol", 1e3}, {"lbmol", 453.59237}});
        add(kCurrent, {{"A", 1.0}, {"mA", 1e-3}});
        add(kVoltage, {{"V", 1.0}, {"mV", 1e-3}, {"kV", 1e3}});
        add(kFrequency, {{"Hz", 1.0}, {"kHz", 1e3}, {"rpm", 2.0 * kPi / 60.0}});
        // angles carry no dimension
        add(kNone, {{"rad", 1.0}, {"deg", kPi / 180.0}, {"rev", 2.0 * kPi}});
        return t;
    }();
    return table;
}

std::string lowerCase(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

/// The names of simpleUnits() that differ only by case: lower-case spelling -> the exact spellings, sorted.
const std::unordered_map<std::string, std::vector<std::string>>& namesByLowerCase() {
    static const std::unordered_map<std::string, std::vector<std::string>> index = [] {
        std::unordered_map<std::string, std::vector<std::string>> byLower;
        for (const auto& entry : simpleUnits()) byLower[lowerCase(entry.first)].push_back(entry.first);
        for (auto& entry : byLower) std::sort(entry.second.begin(), entry.second.end());
        return byLower;
    }();
    return index;
}

/// Whether two strings differ by at most one inserted, deleted or substituted character.
bool withinOneEdit(const std::string& a, const std::string& b) {
    if (a.size() > b.size()) return withinOneEdit(b, a);
    if (b.size() - a.size() > 1) return false;
    size_t i = 0;
    while (i < a.size() && a[i] == b[i]) ++i;
    if (i == a.size()) return true;
    return a.size() == b.size() ? a.compare(i + 1, std::string::npos, b, i + 1, std::string::npos) == 0
                                : a.compare(i, std::string::npos, b, i + 1, std::string::npos) == 0;
}

/// The spelling of `name` in simpleUnits(): the exact one, else the only unit that differs from it by case.
/// nullptr when there is none; `candidates` receives the spellings when the match by case is ambiguous.
const UnitValue* findUnit(const std::string& name, std::vector<std::string>& candidates) {
    const auto& table = simpleUnits();
    const auto exact = table.find(name);
    if (exact != table.end()) return &exact->second;
    const auto folded = namesByLowerCase().find(lowerCase(name));
    if (folded == namesByLowerCase().end()) return nullptr;
    if (folded->second.size() > 1) {
        candidates = folded->second;
        return nullptr;
    }
    return &table.at(folded->second.front());
}

/// A unit string of EES: unit names with `-` or `*` between them as products, `/` as quotient (everything after a `/`
/// is in the denominator: `kJ/kg-K` is kJ/(kg*K)), `^n` for an integer power, parentheses, and the number 1
/// (`1/s`). Spaces are ignored around the names. Parsed by recursive descent.
class UnitStringParser {
public:
    explicit UnitStringParser(const std::string& text) : text_(text) {}

    /// @throws std::runtime_error with the reason when the text is not a unit.
    UnitValue parse() {
        const UnitValue unit = parseProduct();
        if (pos_ < text_.size()) fail("unexpected ')'");
        return unit;
    }

private:
    const std::string& text_;
    size_t pos_ = 0;

    [[noreturn]] void fail(const std::string& detail) const {
        throw std::runtime_error("invalid unit '" + text_ + "': " + detail);
    }

    void skipSpaces() {
        while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_;
    }

    static bool isNameStart(char c) {
        return std::isalpha(static_cast<unsigned char>(c)) || c == '_' || static_cast<unsigned char>(c) >= 0x80;
    }

    UnitValue parseProduct() {
        UnitValue result = parseFactor();
        bool denominator = false;
        for (;;) {
            skipSpaces();
            if (pos_ >= text_.size() || text_[pos_] == ')') return result;
            const char op = text_[pos_];
            if (op != '-' && op != '*' && op != '/')
                fail("'-', '*' or '/' expected before '" + text_.substr(pos_) + "'");
            ++pos_;
            if (op == '/') denominator = true;
            const UnitValue factor = parseFactor();
            if (denominator) {
                result.factor /= factor.factor;
                for (size_t i = 0; i < result.dim.size(); ++i) result.dim[i] -= factor.dim[i];
            } else {
                result.factor *= factor.factor;
                for (size_t i = 0; i < result.dim.size(); ++i) result.dim[i] += factor.dim[i];
            }
        }
    }

    /// A unit name, a group in parentheses or the number 1, with its power.
    UnitValue parseFactor() {
        skipSpaces();
        if (pos_ >= text_.size()) fail("a unit name is expected");
        UnitValue unit;
        const char c = text_[pos_];
        if (c == '(') {
            ++pos_;
            unit = parseProduct();
            skipSpaces();
            if (pos_ >= text_.size()) fail("')' is missing");
            ++pos_;   // the ')' that ended the product
        } else if (std::isdigit(static_cast<unsigned char>(c))) {
            const size_t start = pos_;
            while (pos_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_]))) ++pos_;
            if (text_.compare(start, pos_ - start, "1") != 0) fail("the only number allowed is 1, as in 1/s");
        } else if (isNameStart(c)) {
            const size_t start = pos_;
            while (pos_ < text_.size() && (isNameStart(text_[pos_]) || std::isdigit(static_cast<unsigned char>(text_[pos_])))) ++pos_;
            unit = simpleUnit(text_.substr(start, pos_ - start));
        } else {
            fail("unexpected '" + std::string(1, c) + "' where a unit name is expected");
        }
        skipSpaces();
        if (pos_ < text_.size() && text_[pos_] == '^') {
            ++pos_;
            return raised(unit, parseExponent());
        }
        return unit;
    }

    int parseExponent() {
        skipSpaces();
        int sign = 1;
        if (pos_ < text_.size() && (text_[pos_] == '-' || text_[pos_] == '+')) sign = text_[pos_++] == '-' ? -1 : 1;
        const size_t start = pos_;
        int value = 0;
        while (pos_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_])) && value < 1000)
            value = 10 * value + (text_[pos_++] - '0');
        const bool more = pos_ < text_.size() && (std::isalnum(static_cast<unsigned char>(text_[pos_])) || text_[pos_] == '.');
        if (pos_ == start || more) fail("the exponent after '^' must be an integer");
        return sign * value;
    }

    static UnitValue raised(UnitValue unit, int exponent) {
        unit.factor = std::pow(unit.factor, exponent);
        for (int& e : unit.dim) e *= exponent;
        return unit;
    }

    /// The unit written `name`, or an error that says what is wrong and what may have been meant.
    UnitValue simpleUnit(const std::string& name) const {
        std::vector<std::string> ambiguous;
        if (const UnitValue* unit = findUnit(name, ambiguous)) return *unit;
        if (!ambiguous.empty()) {
            std::string spellings;
            for (const auto& s : ambiguous) spellings += (spellings.empty() ? "" : ", ") + s;
            throw std::runtime_error("the unit '" + name + "' is ambiguous (" + spellings + "): write it with its exact capitalization");
        }
        std::string message = "unknown unit '" + name + "'";
        if (name != text_) message += " in '" + text_ + "'";
        // m3 for m^3
        size_t digits = name.size();
        while (digits > 0 && std::isdigit(static_cast<unsigned char>(name[digits - 1]))) --digits;
        std::vector<std::string> unused;
        if (digits > 0 && digits < name.size() && findUnit(name.substr(0, digits), unused)) {
            message += " (write '" + name.substr(0, digits) + "^" + name.substr(digits) + "' for a power)";
            throw std::runtime_error(message);
        }
        // a name that is one edit away from a unit
        std::vector<std::string> close;
        if (name.size() >= 3)
            for (const auto& entry : simpleUnits())
                if (withinOneEdit(lowerCase(name), lowerCase(entry.first))) close.push_back(entry.first);
        std::sort(close.begin(), close.end());
        if (close.size() > 3) close.resize(3);
        if (!close.empty()) {
            message += " (did you mean ";
            for (size_t i = 0; i < close.size(); ++i) message += (i ? ", '" : "'") + close[i] + "'";
            message += "?)";
        }
        throw std::runtime_error(message);
    }
};

/// What a dimension is called in a message: a name when it has one, else the SI base units (`m^3/kg`).
std::string describeDimension(const Dimension& dim) {
    static const std::vector<std::pair<Dimension, const char*>> names = {
        {kNone, "dimensionless"}, {kMass, "mass"}, {kLength, "length"}, {kTime, "time"},
        {kTemperature, "temperature"}, {kAmount, "amount of substance"}, {kCurrent, "electric current"},
        {kArea, "area"}, {kVolume, "volume"}, {kVelocity, "velocity"}, {{0, 1, -2, 0, 0, 0}, "acceleration"},
        {kFrequency, "frequency"}, {kForce, "force"}, {kPressure, "pressure"}, {kEnergy, "energy"},
        {kPower, "power"}, {kVoltage, "voltage"}, {kViscosity, "dynamic viscosity"},
        {{1, -3, 0, 0, 0, 0}, "density"}, {{1, 0, -1, 0, 0, 0}, "mass flow rate"},
        {{0, 3, -1, 0, 0, 0}, "volume flow rate"}, {{0, 2, -2, 0, 0, 0}, "specific energy"},
        {{0, 2, -2, -1, 0, 0}, "specific heat"}, {{1, 1, -3, -1, 0, 0}, "thermal conductivity"},
        {{1, 0, -3, -1, 0, 0}, "heat transfer coefficient"}, {{1, 0, -3, 0, 0, 0}, "heat flux"}};
    for (const auto& entry : names)
        if (entry.first == dim) return entry.second;
    static const char* const symbols[] = {"kg", "m", "s", "K", "mol", "A"};
    std::string numerator, denominator;
    int denominators = 0;
    for (size_t i = 0; i < dim.size(); ++i) {
        if (dim[i] == 0) continue;
        std::string part = symbols[i];
        if (std::abs(dim[i]) != 1) part += "^" + std::to_string(std::abs(dim[i]));
        if (dim[i] > 0) {
            numerator += (numerator.empty() ? "" : "*") + part;
        } else {
            denominator += (denominator.empty() ? "" : "*") + part;
            ++denominators;
        }
    }
    if (numerator.empty()) numerator = "1";
    if (denominators == 0) return numerator;
    return numerator + "/" + (denominators > 1 ? "(" + denominator + ")" : denominator);
}

}  // namespace

bool isUnitName(const std::string& name) {
    return simpleUnits().count(name) > 0;
}

double unitConversionFactor(const std::string& from, const std::string& to) {
    // a CONVERT with constant units is evaluated at every iteration of the solver: parse the pair once
    thread_local std::unordered_map<std::string, double> factors;
    const std::string key = from + '\x1f' + to;
    const auto known = factors.find(key);
    if (known != factors.end()) return known->second;
    try {
        const UnitValue source = UnitStringParser(from).parse();
        const UnitValue target = UnitStringParser(to).parse();
        if (source.dim != target.dim)
            throw std::runtime_error("incompatible units (" + describeDimension(source.dim) + " vs. " + describeDimension(target.dim) + ")");
        return factors.emplace(key, source.factor / target.factor).first->second;
    } catch (const std::runtime_error& e) {
        throw std::runtime_error("CONVERT('" + from + "','" + to + "'): " + e.what());
    }
}

namespace {

/// The size of one degree of the temperature scale `unit` in K (1 for C and K, 5/9 for F and R): the entry of
/// simpleUnits() (a temperature unit there is a difference of temperature), so the names and the case rule are
/// those of CONVERT.
/// @throws std::runtime_error when `unit` is not a temperature scale.
double temperatureDegree(const std::string& unit) {
    const size_t first = unit.find_first_not_of(" \t");
    const size_t last = unit.find_last_not_of(" \t");
    std::vector<std::string> ambiguous;
    const UnitValue* scale = first == std::string::npos ? nullptr : findUnit(unit.substr(first, last - first + 1), ambiguous);
    if (!scale || scale->dim != kTemperature)
        throw std::runtime_error("CONVERTTEMP: '" + unit + "' is not a temperature unit (expected C, K, F or R, or "
                                 "celsius, kelvin, fahrenheit or rankine)");
    return scale->factor;
}

}  // namespace

void checkTemperatureUnit(const std::string& unit) {
    temperatureDegree(unit);
}

double convertTemperature(double value, const std::string& from, const std::string& to, double* slope) {
    const double fromDegree = temperatureDegree(from);
    const double toDegree = temperatureDegree(to);
    if (slope) *slope = fromDegree / toDegree;
    return UnitConverter::fromSI(UnitConverter::toSI(value, UnitType::Temperature, from), UnitType::Temperature, to);
}

bool unitSettingActive(const std::string& setting, const UnitSystem& units) {
    const std::string s = normalize(setting);
    if (s == "si" || s == "mass" || s == "deg") return true;
    if (s == "eng" || s == "molar" || s == "rad") return false;
    for (const char* t : {"c", "k", "f", "r"})
        if (s == t) return normalize(units.temperature) == s;
    for (const char* p : {"pa", "kpa", "mpa", "bar", "psia", "atm"})
        if (s == p) return normalize(units.pressure) == s;
    for (const char* e : {"j", "kj", "btu", "kcal"})
        if (s == e) return normalize(units.energy) == s;
    throw std::runtime_error("UNITSYSTEM: unknown unit setting '" + setting +
        "' (expected SI, Eng, Mass, Molar, Deg, Rad, C, K, F, R, Pa, kPa, bar, MPa, psia, atm, J, kJ, Btu or kcal)");
}

std::string unitSettingReplacement(const std::string& setting, const UnitSystem& units) {
    if (unitSettingActive(setting, units)) return "";   // throws for an unknown setting
    const std::string s = normalize(setting);
    if (s == "eng") return "CoolSolve works in SI";
    if (s == "molar") return "CoolSolve works on a mass basis";
    if (s == "rad") return "angles are read in degrees";
    for (const char* t : {"c", "k", "f", "r"})
        if (s == t) return std::string("temperatures are read in ") + (normalize(units.temperature) == "k" ? "" : "°") + units.temperature;
    for (const char* p : {"pa", "kpa", "mpa", "bar", "psia", "atm"})
        if (s == p) return "pressures are read in " + units.pressure;
    return "energies are read in " + units.energy;   // J, kJ, Btu, kcal: the settings left
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
