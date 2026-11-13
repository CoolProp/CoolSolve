// The fluid of a chart: the CoolProp fluid string of a model's fluid, the interface the engine asks a fluid's
// states and limits through (a property provider), and the provider of CoolProp's fluids with its constants
// (critical point, triple point, limits of the equation of state, saturation range).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/native.py
// (fluid_string, port_fluid_string), spec/fluid.py (FluidSpec and HumidAirSpec records), thermo/provider.py (the
// provider interface) and thermo/fluid.py (Fluid).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
#pragma once

#include <array>
#include <cmath>
#include <limits>
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
class Fluid;

namespace diagram {

using json = nlohmann::json;

/// A chart request the caller has to correct: an unknown fluid, a fluid not drawable yet, a spec value
/// refused (coolplot.spec raises ValueError for these). The server answers it with HTTP 400.
class DiagramError : public std::runtime_error {
public:
    explicit DiagramError(const std::string& message) : std::runtime_error(message) {}
};

/// A state in SI units by CoolPlot quantity key: "T", "p", "h", "s", "rho", "v", "u", "Q" of a fluid; "T", "W",
/// "RH", "Twb", "Tdp", "h", "v", "Thx", "p" of humid air (per kg of dry air).
using State = std::map<std::string, double>;

/// The keys of a fluid state, in CoolPlot's order (coolplot.core.state.STATE_KEYS).
const std::vector<std::string>& stateKeys();

/// The CoolProp parameter (CoolProp::parameters) of a quantity key of a fluid state other than "v"
/// (thermo.fluid.parameter_index). @throws DiagramError for another key
int coolpropParameter(const std::string& key);

/**
 * @brief The fluid of a chart as CoolPlot describes it (coolplot.spec.FluidSpec): mole fractions of a mixture, mass
 *        fractions of an incompressible solution (CoolProp's INCOMP backend reads them so); or humid air at one total
 *        pressure (coolplot.spec.HumidAirSpec), whose charts are the psychrometric and Mollier charts.
 *
 * CoolSolve's ideal gases are drawn by CoolSolve's own provider (diagram/ideal_gas.h): `idealGas` names them, the
 * components hold the CoolProp fluid of their ideal-gas functions. Humid air (diagram/humid_air.h) has no
 * components: `p_Pa` is its total pressure.
 */
struct ChartFluid {
    std::string backend = "HEOS";            ///< CoolProp backend: "HEOS", "BICUBIC&HEOS", "TTSE&HEOS", "INCOMP", ...
    std::vector<std::string> components;     ///< CoolProp fluid names
    std::vector<double> fractions;           ///< one per component, or empty for a pure fluid
    std::string basis = "mole";              ///< what the fractions are: "mole", or "mass" for the INCOMP backend
    std::string idealGas;                    ///< CoolSolve's name of an ideal gas (FluidRegistry), empty for CoolProp's fluids
    double p_Pa = std::numeric_limits<double>::quiet_NaN();  ///< humid air: the total pressure; NaN for the other fluids

    /// The CoolProp fluid string, e.g. "HEOS::R134a", "HEOS::R32[0.7]&R125[0.3]", "INCOMP::MEG[0.3]"
    /// (FluidSpec.__str__); "IDEALGAS::N2" for an ideal gas and "HUMIDAIR::101325" for humid air, which CoolProp
    /// never sees (keys of the caches).
    std::string fluidString() const;
    /// The short name for messages and titles, e.g. "R134a", "MEG[0.3]" (FluidSpec.name); "N2" for an ideal gas,
    /// "humid air at 101325 Pa" (HumidAirSpec.name).
    std::string name() const;
    /// The fluid string the port gets for a line of constant `key`: the lines of constant quality of a tabular
    /// backend (BICUBIC&HEOS, TTSE&HEOS) come from the equation of state behind the tables (native.port_fluid_string).
    std::string portFluidString(const std::string& key) const;
    /// Whether CoolSolve's ideal-gas provider draws it.
    bool isIdealGas() const { return !idealGas.empty(); }
    /// Whether it is humid air (drawn by the humid-air provider in psychrometric charts).
    bool isHumidAir() const { return !std::isnan(p_Pa); }
    /// The kind of fluid, as GET /api/v1/coolprop/fluids names it: "Real", "Mixture", "Incompressible", "IdealGas",
    /// "HumidAir".
    std::string kind() const;
    /// The fluid record of a CoolPlot spec (FluidSpec.to_dict; HumidAirSpec.to_dict: {"model": "humid_air", "p_Pa"}).
    /// An ideal gas adds CoolSolve's namespaced {"coolsolve": {"model": "ideal_gas", "fluid": name}}: CoolPlot ignores
    /// it and would draw the real fluid.
    json record() const;
    /// A chart fluid from a fluid record (FluidSpec.from_dict, fluid_spec_from_dict); throws DiagramError for a record
    /// it cannot use.
    static ChartFluid fromRecord(const json& record);
};

/// The name of humid air on charts and in the states of a model: EES's substance name. FluidRegistry's humid-air fluid
/// calls itself "Air_ha" (HumidAirFluid::getName), which is the name of the real gas Air_ha.
constexpr const char* HUMID_AIR_FLUID = "AirH2O";

/// The name of a CoolSolve fluid in the Diagrams tab (the fluids listed, the states of a model): FluidRegistry's name,
/// HUMID_AIR_FLUID for humid air.
std::string diagramFluidName(const Fluid& fluid);

/// The mass fractions an incompressible solution of CoolProp can have (its fraction_min, fraction_max), e.g.
/// (0, 0.6) for "MEG". @throws DiagramError for a fluid CoolProp's INCOMP backend does not know
std::pair<double, double> solutionFractionRange(const std::string& coolpropKey);

/**
 * @brief The chart fluid of a CoolSolve fluid name (FluidRegistry) with the model's CoolProp backend.
 *
 * `Real` fluids and CoolProp mixtures take the backend; an incompressible solution is drawn with CoolProp's INCOMP
 * backend at the mass fraction `fractions` (one value, `basis` "mass" or empty); an ideal gas by CoolSolve's
 * ideal-gas provider (EES/JANAF reference states, diagram/ideal_gas.h); humid air (AirH2O) by the humid-air provider
 * (diagram/humid_air.h) at the total pressure `p_Pa`, 101325 Pa when it is not given.
 * @throws DiagramError for an unknown name, a solution without its fraction or with one outside its range, fractions
 *         or a total pressure given to another fluid, a total pressure outside the range of HAPropsSI, and for
 *         NH3H2O, which is "not drawable yet" (docs/coolplot_integration.md §5)
 */
ChartFluid chartFluid(const std::string& name, const std::string& backend = "HEOS", const std::vector<double>& fractions = {},
                      const std::string& basis = "", std::optional<double> p_Pa = std::nullopt);

/// The T-p box of the lines of an incompressible liquid where the spec gives no absolute pressure (CoolSolve's
/// default: CoolPlot asks for absolute limits): 1 bar to 100 bar. Its temperatures are those of the fluid.
constexpr double LIQUID_P_MIN_PA = 1e5;
constexpr double LIQUID_P_MAX_PA = 1e7;

/**
 * @brief What the engine asks of the fluid of a chart (CoolPlot's property provider, thermo/provider.py): full states
 *        from two quantities, the limits of the model, its critical point and its saturation range.
 *
 * Three providers exist: FluidModel (CoolProp, CoolPlot's thermo.fluid.Fluid), IdealGasModel (CoolSolve's ideal
 * gases, diagram/ideal_gas.h) and HumidAirModel (HAPropsSI, CoolPlot's thermo.humid_air.HumidAir,
 * diagram/humid_air.h). The defaults below are those of a fluid without a saturation curve. Not thread safe:
 * every calculation updates the state. Values that are not known are NaN (Python: None).
 */
class PropertyProvider {
public:
    virtual ~PropertyProvider() = default;

    /// The name for messages (FluidSpec.name).
    virtual const std::string& name() const = 0;
    /// Whether the fluid has one component (a predefined blend counts by its components).
    virtual bool isPure() const = 0;

    /// The keys of a state of this provider (provider.quantities): stateKeys() for a fluid.
    virtual const std::vector<std::string>& quantities() const { return stateKeys(); }
    /// The keys a state may be given by (provider.inputs): those of a state of a fluid; humid air not its pressure.
    virtual const std::vector<std::string>& inputs() const { return stateKeys(); }
    /// How the path between two states of a cycle is drawn (provider.path_kind): `kind` for a fluid ("auto" is
    /// resolved from the states later), "straight" for humid air.
    virtual std::string pathKind(const State& a, const State& b, const std::string& kind) const;
    /// Why a state computed from a solved model does not belong to this provider's chart (chart.check_state: humid
    /// air at another total pressure); empty when it does.
    virtual std::string stateProblem(const State& state) const;

    /**
     * @brief Set the working state from two quantities (Fluid.update); a specific volume is given as a density.
     * @throws std::exception when the state cannot be computed
     */
    virtual void update(const std::string& key1, double value1, const std::string& key2, double value2) = 0;
    /// One quantity of the state of the last update(), in SI units ("v" is 1/rho).
    virtual double get(const std::string& key) const = 0;
    /// The full state from two quantities (Fluid.state_point): every key of quantities().
    State statePoint(const std::string& key1, double value1, const std::string& key2, double value2);

    virtual double tTripleK() = 0;   ///< triple-point temperature; for a mixture the highest one of its components
    virtual double tMinK() = 0;      ///< lowest temperature of the model; for a mixture the highest one
    virtual double tMaxK() = 0;      ///< highest temperature of the model
    virtual double pMinPa() = 0;     ///< lowest pressure of the model
    virtual double pMaxPa() = 0;     ///< highest pressure of the model

    /// Whether the fluid has a critical point, i.e. a saturation curve and a dome.
    virtual bool hasCriticalPoint() { return false; }
    /// T, p, h, s, rho, u at the critical point (NaN where unknown). @throws std::runtime_error without one.
    virtual const State& critical();
    /// For CoolProp's pseudo-pure model of a blend ("R407C"), the blend as a mixture ("R407C.mix"); else empty.
    virtual std::string pseudoPureBlend() { return ""; }
    /**
     * @brief Lowest and highest saturation temperature ("T") or pressure ("p"), `margin` (relative to the
     *        critical values) inside the triple point and the critical point (Fluid.saturation_range).
     * @throws DiagramError for a fluid without a saturation curve or another key
     */
    virtual std::pair<double, double> saturationRange(const std::string& key, double margin);

    /// The pressure of a state given without one, NaN for none: an incompressible liquid is evaluated at
    /// IncompressibleFluid::kInternalPressureSI by the solver's property calls whose only state input is T.
    virtual double defaultPressure() const { return std::numeric_limits<double>::quiet_NaN(); }

    /**
     * @brief The T-p box (T_min_K, T_max_K, p_min_Pa, p_max_Pa) of the lines of a fluid without a saturation curve,
     *        where the spec gives no absolute bound (CoolSolve's: CoolPlot asks for absolute limits for these fluids,
     *        docs/coolplot_integration.md §5). Not used for a fluid with a saturation curve.
     */
    virtual std::array<double, 4> defaultBox() = 0;
};

/**
 * @brief A CoolProp fluid model with its own state and cached constants (coolplot.thermo.fluid.Fluid).
 *
 * Not thread safe: every calculation updates the state. The engine keeps one per fluid string and locks it.
 */
class FluidModel : public PropertyProvider {
public:
    explicit FluidModel(const ChartFluid& fluid);
    ~FluidModel() override;
    FluidModel(const FluidModel&) = delete;
    FluidModel& operator=(const FluidModel&) = delete;

    const std::string& name() const override { return name_; }
    bool isPure() const override { return components_.size() <= 1; }

    /**
     * @brief Set the working state from two quantities (Fluid.update): a specific volume is given as a density;
     *        a mixture state that CoolProp cannot find is tried on a scratch state with the phase envelope
     *        built, and a mixture state reported single phase inside the dome is refused.
     * @throws std::exception (CoolProp's) when the state cannot be computed
     */
    void update(const std::string& key1, double value1, const std::string& key2, double value2) override;
    double get(const std::string& key) const override;

    double tTripleK() override;
    /// Lowest temperature of the equation of state; for an incompressible solution not below its freezing point
    /// (CoolSolve's: CoolProp refuses the states below it).
    double tMinK() override;
    double tMaxK() override;
    double pMinPa() override;
    double pMaxPa() override;

    /// Whether CoolProp provides a critical point; incompressibles have none (Fluid.has_critical_point).
    bool hasCriticalPoint() override;
    const State& critical() override;
    std::string pseudoPureBlend() override;
    std::pair<double, double> saturationRange(const std::string& key, double margin) override;
    /// The temperatures of the fluid and the pressures LIQUID_P_MIN_PA to LIQUID_P_MAX_PA (an incompressible).
    std::array<double, 4> defaultBox() override;
    /// IncompressibleFluid::kInternalPressureSI for an incompressible, NaN otherwise.
    double defaultPressure() const override;

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
    std::string basis_;
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
