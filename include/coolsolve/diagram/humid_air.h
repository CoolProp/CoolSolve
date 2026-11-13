// Humid air for psychrometric charts (CS-FEAT-PSYCHRO): the property provider of humid air at one total pressure,
// computed with CoolProp's HAPropsSI (the CoolProp the solver uses), the domain of a chart (its T-W or Thx-W region)
// and the lines of the psychrometric chart (humidity ratio W over dry-bulb temperature T, ASHRAE) and of the Mollier
// h,x chart (its temperature scale Thx over W, DIN 1946): lines of constant relative humidity RH, wet-bulb
// temperature Twb, enthalpy h and specific volume v (both per kg of dry air), and on the Mollier chart isotherms.
// The saturation line (RH = 1) bounds the chart: every line is cut where it meets saturation, the meeting point
// found exactly. All values in SI units: K, kg/kg, a fraction 0 to 1 for RH, J/kg and m3/kg per kg of dry air.
//
// The Mollier chart's vertical coordinate is Thx = 273.15 K + (h - r0 W) / cp_a, r0 = 2500 kJ/kg,
// cp_a = 1.006 kJ/(kg K), read as a temperature: at W = 0 it is the dry-bulb temperature, and mixing lines are
// straight in it.
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/humid_air.py (HumidAir,
// air_state, psychro_line and its helpers, chart_ranges, _dry_bulb_span, saturation_W), core/limits.py
// (HUMID_AIR_T_LIMITS_K, HUMID_AIR_P_LIMITS_PA, within) and core/state.py (AIR_STATE_KEYS, R0_J_KG,
// CP_DRY_AIR_J_KG_K).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
// See docs/coolplot_integration.md §5 (phase 3).
#pragma once

#include "coolsolve/diagram/engine.h"  // Curve

#include <string>
#include <vector>

namespace coolsolve {
namespace diagram {

/// What HAPropsSI covers: temperatures in K (core.limits.HUMID_AIR_T_LIMITS_K) ...
constexpr double HUMID_AIR_T_MIN_K = 173.15;
constexpr double HUMID_AIR_T_MAX_K = 623.15;
/// ... and total pressures in Pa (HUMID_AIR_P_LIMITS_PA).
constexpr double HUMID_AIR_P_MIN_PA = 1e4;
constexpr double HUMID_AIR_P_MAX_PA = 1e7;
/// The total pressure of a humid-air chart that is given none (HumidAirSpec, sea level).
constexpr double DEFAULT_HUMID_AIR_P_PA = 101325.0;

/// The enthalpy of evaporation of water at 0 °C and the specific heat of dry air of the Mollier chart's temperature
/// scale (core.state.R0_J_KG, CP_DRY_AIR_J_KG_K).
constexpr double R0_J_KG = 2500e3;
constexpr double CP_DRY_AIR_J_KG_K = 1006.0;

/// The keys of a state of humid air, in CoolPlot's order (core.state.AIR_STATE_KEYS): T, W, RH, Twb, Tdp, h, v, Thx, p.
const std::vector<std::string>& airStateKeys();

/// The keys a state of humid air may be given by: all but the total pressure (humid_air.HUMID_AIR).
const std::vector<std::string>& airInputKeys();

/// The total pressure `p_Pa` when HAPropsSI covers it. @throws DiagramError otherwise (HumidAir.__init__)
double checkedHumidAirPressure(double p_Pa);

/// "humid air at 101325 Pa" (HumidAirSpec.name), for messages and titles.
std::string humidAirName(double p_Pa);

/// The temperature scale of the Mollier chart, in K, of the state with enthalpy `h_J_kg` (per kg of dry air) and
/// humidity ratio `W` (humid_air.mollier_temperature).
double mollierTemperature(double h_J_kg, double W);

/// The enthalpy (J/kg of dry air) of the point (W, Thx) of a Mollier chart (humid_air.enthalpy_from_mollier).
double enthalpyFromMollier(double Thx_K, double W);

/**
 * @brief One property of humid air from two others at total pressure `p_Pa` (humid_air.humid_air): HAPropsSI with the
 *        names of CoolPlot's quantities ("T", "W", "RH", "Twb", "Tdp", "h", "v").
 * @throws std::runtime_error when HAPropsSI has no value (Python's HAPropsSI raises ValueError), DiagramError for a
 *         quantity humid air does not have
 */
double humidAir(const std::string& output, const std::string& key1, double value1, const std::string& key2, double value2,
                double p_Pa);

/// The humidity ratio of saturated air at `T_K`, NaN where HAPropsSI has none (above the boiling point) (saturation_W).
double saturationW(double T_K, double p_Pa);

/**
 * @brief The full state of humid air from two properties at total pressure `p_Pa` (humid_air.air_state), every key
 *        of airStateKeys(), per kg of dry air: e.g. airState("T", 298.15, "RH", 0.5, 101325). Thx is given together
 *        with W only; the pressure is the chart's, not one of the two.
 * @throws DiagramError for two properties that cannot be given together, std::runtime_error when HAPropsSI has no state
 */
State airState(const std::string& key1, double value1, const std::string& key2, double value2, double p_Pa);

/**
 * @brief The (min, max) of RH, Twb, h and v over the visible part of a chart, the region above saturation (fog) left
 *        out (humid_air.chart_ranges): from the corners of the region below saturation. RH is (0, 1); a quantity
 *        HAPropsSI gives at no corner is left out.
 */
std::vector<std::pair<std::string, Range>> chartRanges(double p_Pa, const Range& T_K, const Range& W);

/**
 * @brief A line of constant `key` (RH, Twb, h or v) in the T-W plane, in SI units (humid_air.psychro_line): sampled
 *        at `points` dry-bulb temperatures evenly over `T_K`. Points above saturation (fog) or below W = 0 are left
 *        out; the point where the line meets saturation is added exactly, so lines start on the saturation curve (a
 *        wet-bulb line at T = Twb), and so is the point where a falling line reaches W = 0. Above the boiling point at
 *        `p_Pa` there is no saturation and lines continue. `n_failed` counts the points where HAPropsSI found no state.
 * @throws DiagramError for another key
 */
Curve psychroLine(const std::string& key, double value_SI, double p_Pa, const Range& T_K, int points);

/**
 * @brief The property provider of humid air at one total pressure (coolplot.thermo.humid_air.HumidAir): states from
 *        two properties at the chart's pressure, the domain of a chart, its lines. Every call goes to HAPropsSI, whose
 *        working states are kept per thread by CoolProp: line() and domain() may run on several threads; update()
 *        keeps the last state, as the other providers do.
 *
 * It has no critical point, no vapour quality and no saturation range; its temperatures are those of HAPropsSI and
 * its pressure is the chart's. Its cycles are drawn with straight paths between their states.
 */
class HumidAirModel : public PropertyProvider {
public:
    /// The provider at total pressure `p_Pa`. @throws DiagramError outside the pressures of HAPropsSI
    explicit HumidAirModel(double p_Pa);

    const std::string& name() const override { return name_; }
    bool isPure() const override { return false; }
    const std::vector<std::string>& quantities() const override { return airStateKeys(); }
    const std::vector<std::string>& inputs() const override { return airInputKeys(); }
    /// Air processes are drawn as straight lines between their states (HumidAir.path_kind).
    std::string pathKind(const State& a, const State& b, const std::string& kind) const override;
    /// A state at another total pressure than the chart's (within 1e-6, chart.check_state).
    std::string stateProblem(const State& state) const override;

    /// The state of humid air from two properties (airState). @throws std::exception when it has none
    void update(const std::string& key1, double value1, const std::string& key2, double value2) override;
    double get(const std::string& key) const override;

    double tTripleK() override;
    double tMinK() override { return HUMID_AIR_T_MIN_K; }
    double tMaxK() override { return HUMID_AIR_T_MAX_K; }
    double pMinPa() override { return p_Pa_; }
    double pMaxPa() override { return p_Pa_; }
    std::array<double, 4> defaultBox() override { return {HUMID_AIR_T_MIN_K, HUMID_AIR_T_MAX_K, p_Pa_, p_Pa_}; }
    /// A state of a model given without its pressure is at the chart's total pressure.
    double defaultPressure() const override { return p_Pa_; }

    /// The total pressure, Pa.
    double pressure() const { return p_Pa_; }

    /**
     * @brief The calculation domain of a chart over `region`, bounds on T and W, or on Thx and W for the Mollier chart
     *        (HumidAir.domain): the T-p box (flat in pressure), the ranges of T, W, p, of Thx for the Mollier chart
     *        (whose dry-bulb temperatures are swept a little wider, _dry_bulb_span) and of RH, Twb, h and v over the
     *        part of the region below saturation (chartRanges), and a saturation curve.
     * @throws DiagramError for temperatures outside those of HAPropsSI, a negative humidity ratio, or a region entirely
     *         above saturation
     */
    Domain domain(const Region& region) const;

    /**
     * @brief A line of constant `key` over the temperatures of the domain's T-p box, in the axes of `diagram`
     *        (HumidAir.line): (T, W) on the psychrometric chart; (W, Thx) on the Mollier chart, where isotherms are
     *        lines too, from dry air to saturation and at most 1.2 times the top of the domain's W range.
     */
    Curve line(const std::string& key, double value_SI, const Domain& domain, int points, const coolplot::DiagramType& diagram) const;

private:
    Curve isotherm(double T_K, int points, const std::optional<Range>& W) const;
    Curve toMollier(Curve curve) const;

    double p_Pa_;
    std::string name_;
    State state_;
};

}  // namespace diagram
}  // namespace coolsolve
