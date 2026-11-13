// The charts of CoolSolve's ideal gases (CS-FEAT-DIAGRAM-IDEAL): a property provider over CoolSolve's ideal-gas
// functions (ideal_gas.h, IdealGasFluid) with the EES/JANAF reference states of the solver, and the lines of a chart
// computed from it, so that the solved states of an ideal-gas model lie on their lines.
//
// CoolSolve's own: CoolPlot draws CoolProp's fluids only (its providers are thermo/fluid.py and thermo/humid_air.py),
// and CoolProp::Plot has no ideal-gas model. The provider follows CoolPlot's provider interface (thermo/provider.py,
// PropertyProvider in fluid.h), and the curves are CoolPlot's curve records (thermo/isolines.py: Curve).
// See docs/coolplot_integration.md §5 (phase 3).
#pragma once

#include "coolsolve/diagram/engine.h"  // Curve
#include "coolsolve/fluids.h"           // ReferenceState
#include "coolsolve/ideal_gas.h"        // IdealGasPoint

#include <array>
#include <string>

namespace coolsolve {
namespace diagram {

/// The T-p box of the lines of an ideal gas where the spec gives no absolute bound: 200 K to 2000 K, 0.1 bar to
/// 100 bar (T_min_K, T_max_K, p_min_Pa, p_max_Pa). Absolute limits in the spec replace any of them.
constexpr std::array<double, 4> IDEAL_GAS_BOX = {200.0, 2000.0, 1e4, 1e7};

/**
 * @brief The property provider of a CoolSolve ideal gas (IdealGasFluid): pv = RT, h and u functions of T only, on the
 *        EES/JANAF basis of the solver (IdealGasFluid::getReferenceState added to h, u and s).
 *
 * States come from any two independent quantities of T, p, h, s, rho (or v) and u, as the solver's ideal-gas
 * functions take them (two of T, h, u fix only the temperature). The model has no saturation curve, no vapour
 * quality (Q is NaN) and the temperature range of the ideal-gas solve (kIdealGasSolveTMinSI to kIdealGasSolveTMaxSI);
 * its pressure is unbounded. Not thread safe (one state); cheap to make.
 */
class IdealGasModel : public PropertyProvider {
public:
    /// The provider of `fluid`, a chart fluid with `idealGas` set (chartFluid). @throws DiagramError for another
    explicit IdealGasModel(const ChartFluid& fluid);

    const std::string& name() const override { return name_; }
    bool isPure() const override { return true; }
    /// @throws std::exception for a quantity that is not a state of an ideal gas (Q), two of T, h and u, or a state
    ///         the ideal-gas solve cannot reach
    void update(const std::string& key1, double value1, const std::string& key2, double value2) override;
    double get(const std::string& key) const override;

    double tTripleK() override;
    double tMinK() override { return kIdealGasSolveTMinSI; }
    double tMaxK() override { return kIdealGasSolveTMaxSI; }
    double pMinPa() override;
    double pMaxPa() override;
    std::array<double, 4> defaultBox() override { return IDEAL_GAS_BOX; }

private:
    std::string name_;
    std::string coolprop_;     ///< the CoolProp fluid of the ideal-gas functions
    ReferenceState reference_;  ///< EES - CoolProp offsets of h (and u) and s
    double T_;
    double p_;
    IdealGasPoint point_{};
};

/**
 * @brief Whether an ideal gas has a chart with these axes: not when both fix only the temperature (two of T, h and
 *        u: the T-h chart is a single line). @throws DiagramError naming the charts it has
 */
void checkIdealGasDiagram(const coolplot::DiagramType& diagram);

/**
 * @brief One line of constant `key` = `value_SI` of an ideal gas in a chart, as CoolPlot's curve record: `points`
 *        states spaced evenly in the quantity the family sweeps (coolplot::Family::sweep: the pressure on a
 *        logarithmic scale, the temperature linearly) over the T-p box of `domain`. The axes are the ranges of the
 *        domain: a point outside them is NaN unless it neighbours a point inside, so that the line reaches the edge
 *        of the plot area, where the scene clips it. A line of constant vapour quality has no point and its reason.
 */
Curve idealGasLine(IdealGasModel& gas, const coolplot::DiagramType& diagram, const std::string& key, double value_SI,
                   const Domain& domain, int points);

}  // namespace diagram
}  // namespace coolsolve
