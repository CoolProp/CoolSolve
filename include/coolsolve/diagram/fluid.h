// The fluid of a chart: the CoolProp fluid string of a model's fluid, and a CoolProp fluid model with its
// constants (critical point, triple point, limits of the equation of state, saturation range).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/native.py
// (fluid_string, port_fluid_string), spec/fluid.py (FluidSpec record) and thermo/fluid.py (Fluid).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#pragma once

#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace CoolProp {
class AbstractState;
}

namespace coolsolve {
namespace diagram {

using json = nlohmann::json;

/// A chart request the caller has to correct: an unknown fluid, a fluid not drawable yet, a spec value
/// refused (coolplot.spec raises ValueError for these). The server answers it with HTTP 400.
class DiagramError : public std::runtime_error {
public:
    explicit DiagramError(const std::string& message) : std::runtime_error(message) {}
};

/// A state in SI units by CoolPlot quantity key ("T", "p", "h", "s", "rho", "v", "u", "Q").
using State = std::map<std::string, double>;

/// The keys of a fluid state, in CoolPlot's order (coolplot.core.state.STATE_KEYS).
const std::vector<std::string>& stateKeys();

/// The CoolProp parameter (CoolProp::parameters) of a quantity key of a fluid state other than "v"
/// (thermo.fluid.parameter_index). @throws DiagramError for another key
int coolpropParameter(const std::string& key);

/**
 * @brief The fluid of a chart as CoolPlot describes it (coolplot.spec.FluidSpec, mole fractions only).
 */
struct ChartFluid {
    std::string backend = "HEOS";            ///< CoolProp backend: "HEOS", "BICUBIC&HEOS", "TTSE&HEOS", ...
    std::vector<std::string> components;     ///< CoolProp fluid names
    std::vector<double> fractions;           ///< mole fractions of a mixture; empty for a pure fluid

    /// The CoolProp fluid string, e.g. "HEOS::R134a" or "HEOS::R32[0.7]&R125[0.3]" (FluidSpec.__str__).
    std::string fluidString() const;
    /// The short name for messages and titles, e.g. "R134a" (FluidSpec.name).
    std::string name() const;
    /// The fluid string the port gets for a line of constant `key`: the lines of constant quality of a tabular
    /// backend (BICUBIC&HEOS, TTSE&HEOS) come from the equation of state behind the tables (native.port_fluid_string).
    std::string portFluidString(const std::string& key) const;
    /// The fluid record of a CoolPlot spec (FluidSpec.to_dict).
    json record() const;
    /// A chart fluid from a fluid record (FluidSpec.from_dict); throws DiagramError for a record it cannot use.
    static ChartFluid fromRecord(const json& record);
};

/**
 * @brief The chart fluid of a CoolSolve fluid name (FluidRegistry) with the model's CoolProp backend.
 *
 * Phase 1 draws the fluids CoolProp's port computes from a fluid string: `Real` fluids and CoolProp mixtures.
 * @throws DiagramError for an unknown name, and for ideal gases, humid air, incompressibles and NH3H2O, which are
 *         "not drawable yet" (docs/coolplot_integration.md §5, phase 3).
 */
ChartFluid chartFluid(const std::string& name, const std::string& backend = "HEOS");

/**
 * @brief A CoolProp fluid model with its own state and cached constants (coolplot.thermo.fluid.Fluid).
 *
 * Not thread safe: every calculation updates the state. The engine keeps one per fluid string and locks it.
 * Values that CoolProp does not know are NaN (Python: None).
 */
class FluidModel {
public:
    explicit FluidModel(const ChartFluid& fluid);
    ~FluidModel();
    FluidModel(const FluidModel&) = delete;
    FluidModel& operator=(const FluidModel&) = delete;

    /// The name for messages (FluidSpec.name).
    const std::string& name() const { return name_; }
    /// Whether the fluid has one component (a predefined blend counts by its components).
    bool isPure() const { return components_.size() <= 1; }

    /**
     * @brief Set the working state from two quantities (Fluid.update): a specific volume is given as a density;
     *        a mixture state that CoolProp cannot find is tried on a scratch state with the phase envelope
     *        built, and a mixture state reported single phase inside the dome is refused.
     * @throws std::exception (CoolProp's) when the state cannot be computed
     */
    void update(const std::string& key1, double value1, const std::string& key2, double value2);
    /// One quantity of the state of the last update(), in SI units ("v" is 1/rho).
    double get(const std::string& key) const;
    /// The full state from two quantities (Fluid.state_point): every key of stateKeys().
    State statePoint(const std::string& key1, double value1, const std::string& key2, double value2);

    double tTripleK();   ///< triple-point temperature; for a mixture the highest one of its components
    double tMinK();      ///< lowest temperature of the equation of state; for a mixture the highest one
    double tMaxK();      ///< highest temperature of the equation of state
    double pMinPa();     ///< lowest pressure of the equation of state
    double pMaxPa();     ///< highest pressure of the equation of state

    /// Whether CoolProp provides a critical point, i.e. a saturation curve and a dome.
    bool hasCriticalPoint();
    /// T, p, h, s, rho, u at the critical point (NaN where unknown). @throws std::runtime_error without one.
    const State& critical();
    /// For CoolProp's pseudo-pure model of a blend ("R407C"), the blend as a mixture ("R407C.mix"); else empty.
    std::string pseudoPureBlend();

    /**
     * @brief Lowest and highest saturation temperature ("T") or pressure ("p"), `margin` (relative to the
     *        critical values) inside the triple point and the critical point (Fluid.saturation_range).
     * @throws DiagramError for a fluid without a saturation curve or another key
     */
    std::pair<double, double> saturationRange(const std::string& key, double margin);

private:
    std::shared_ptr<CoolProp::AbstractState> makeState() const;
    std::shared_ptr<CoolProp::AbstractState> guideState();
    std::string inconsistency(CoolProp::AbstractState& state);
    double trivial(int key) const;
    double lowestTemperature(int key);

    std::string name_;
    std::string backend_;
    std::vector<std::string> components_;
    std::vector<double> fractions_;
    std::shared_ptr<CoolProp::AbstractState> state_;
    std::shared_ptr<CoolProp::AbstractState> current_;
    std::shared_ptr<CoolProp::AbstractState> guide_;
    bool guideTried_ = false;
    std::optional<State> critical_;
    std::optional<bool> hasCritical_;
    std::optional<std::string> pseudoPure_;
};

}  // namespace diagram
}  // namespace coolsolve
