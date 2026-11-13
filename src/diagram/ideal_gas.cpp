// The charts of CoolSolve's ideal gases (include/coolsolve/diagram/ideal_gas.h).
//
// CoolSolve's own: CoolPlot has no ideal-gas provider. The provider follows CoolPlot's provider interface
// (thermo/provider.py), the lines are CoolPlot's curve records (thermo/isolines.py: Curve).
#include "coolsolve/diagram/ideal_gas.h"

#include "CoolProp/DataStructures.h"

#include <cmath>
#include <limits>

namespace coolsolve {
namespace diagram {

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

/// Whether a quantity fixes the temperature of an ideal gas by itself.
bool fixesTemperature(const std::string& key) {
    return key == "T" || key == "h" || key == "u";
}

}  // namespace

// ---------------------------------------------------------------------------
// IdealGasModel
// ---------------------------------------------------------------------------

IdealGasModel::IdealGasModel(const ChartFluid& fluid) : name_(fluid.name()), T_(NaN), p_(NaN) {
    const std::shared_ptr<Fluid> registered = FluidRegistry::getFluid(fluid.idealGas);
    const auto gas = std::dynamic_pointer_cast<IdealGasFluid>(registered);
    if (!gas) throw DiagramError("'" + fluid.idealGas + "' is not an ideal gas of CoolSolve.");
    coolprop_ = gas->getCoolPropName();
    reference_ = gas->getReferenceState();
}

void IdealGasModel::update(const std::string& key1, double value1, const std::string& key2, double value2) {
    // To the CoolProp basis of the ideal-gas functions: the reference offsets off, a specific volume as a density
    auto input = [this](const std::string& given, double value) -> std::pair<CoolProp::parameters, double> {
        const std::string key = coolplot::quantity_key(given);
        if (key == "v") {
            if (!(value > 0.0)) throw std::runtime_error("A specific volume has to be positive.");
            return {CoolProp::iDmass, 1.0 / value};
        }
        if (key == "Q") throw std::runtime_error(name_ + " is an ideal gas: it has no vapour quality.");
        if (key == "h" || key == "u") value -= reference_.h_offset;
        if (key == "s") value -= reference_.s_offset;
        return {static_cast<CoolProp::parameters>(coolpropParameter(key)), value};
    };
    const auto [k1, v1] = input(key1, value1);
    const auto [k2, v2] = input(key2, value2);
    T_ = NaN;
    p_ = NaN;
    const double T = idealGasProperty(coolprop_, CoolProp::iT, k1, v1, k2, v2, false).value;
    const double p = idealGasProperty(coolprop_, CoolProp::iP, k1, v1, k2, v2, false).value;
    point_ = idealGasAt(coolprop_, T);
    T_ = T;
    p_ = p;
}

double IdealGasModel::get(const std::string& key) const {
    const std::string k = coolplot::quantity_key(key);
    const double R = point_.R;
    if (k == "T") return T_;
    if (k == "p") return p_;
    if (k == "h") return point_.h + reference_.h_offset;
    if (k == "u") return point_.h - R * T_ + reference_.h_offset;
    if (k == "s") return point_.s - R * std::log(p_ / kIdealGasStandardPressureSI) + reference_.s_offset;
    if (k == "rho") return p_ / (R * T_);
    if (k == "v") return R * T_ / p_;
    if (k == "Q") return NaN;
    throw DiagramError("'" + key + "' is not a quantity of a fluid state.");
}

double IdealGasModel::tTripleK() { return NaN; }
double IdealGasModel::pMinPa() { return NaN; }
double IdealGasModel::pMaxPa() { return NaN; }

// ---------------------------------------------------------------------------
// Charts and lines
// ---------------------------------------------------------------------------

void checkIdealGasDiagram(const coolplot::DiagramType& diagram) {
    if (fixesTemperature(diagram.x) && fixesTemperature(diagram.y))
        throw DiagramError("An ideal gas has no " + diagram.y + "-" + diagram.x + " chart: " + diagram.y + " and " + diagram.x +
                           " both depend on the temperature alone, so every state lies on one line. Its charts are p-h, "
                           "T-s, h-s, p-s, p-rho, T-rho and p-T.");
}

Curve idealGasLine(IdealGasModel& gas, const coolplot::DiagramType& diagram, const std::string& key, double value_SI,
                   const Domain& domain, int points) {
    Curve curve;
    curve.method = "ideal_gas";
    if (key == "Q") {
        curve.reason = gas.name() + " is an ideal gas: it has no two-phase region, so lines of constant quality cannot be drawn";
        return curve;
    }
    const coolplot::Family* family = diagram.family(key);
    if (!family) throw DiagramError("Lines of constant '" + key + "' cannot be drawn in a " + diagram.name() + " diagram.");
    const std::array<double, 4>& box = domain.tp_box;
    const bool byPressure = family->sweep == "p";
    const std::vector<double> sweep = byPressure ? logGrid(box[2], box[3], points) : linearGrid(box[0], box[1], points);
    const Range xRange = domain.rangeSI(diagram.x), yRange = domain.rangeSI(diagram.y);
    auto within = [](double v, const Range& r) {
        const double tolerance = 1e-9 * (r.second - r.first);
        return v >= r.first - tolerance && v <= r.second + tolerance;
    };
    std::vector<bool> inside(sweep.size(), false);
    curve.x.assign(sweep.size(), NaN);
    curve.y.assign(sweep.size(), NaN);
    for (size_t i = 0; i < sweep.size(); ++i) {
        try {
            gas.update(key, value_SI, family->sweep, sweep[i]);
            curve.x[i] = gas.get(diagram.x);
            curve.y[i] = gas.get(diagram.y);
        } catch (const std::exception&) {
            ++curve.n_failed;
            continue;
        }
        inside[i] = within(curve.x[i], xRange) && within(curve.y[i], yRange);
        if (!inside[i]) ++curve.n_outside;
    }
    curve.n_evaluated = static_cast<int>(sweep.size());
    // A point outside the axes is drawn only next to one inside: the line ends at the edge of the plot area
    for (size_t i = 0; i < sweep.size(); ++i) {
        const bool nextToInside = (i > 0 && inside[i - 1]) || (i + 1 < sweep.size() && inside[i + 1]);
        if (!inside[i] && !nextToInside) curve.x[i] = curve.y[i] = NaN;
    }
    return curve;
}

}  // namespace diagram
}  // namespace coolsolve
