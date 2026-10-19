#pragma once

#include "coolsolve/taylor_jet.h"

#include <string>

namespace coolsolve::nh3h2o {

/**
 * @file nh3h2o.h
 * @brief The thermodynamic core of the ammonia - water mixture: the IAPWS Guideline G4-01
 *        (2001) Helmholtz formulation of Tillner-Roth & Friend (1998).
 *
 * The molar Helmholtz free energy is f = R T Phi with the dimensionless Phi = Phi^o + Phi^r of
 * Eqs. (1)-(8) of the guideline: an ideal part in the common reduced variables tau0 = T0/T and
 * delta0 = rho/rho0 (Eq. 2, the pure-fluid ideal parts refit to the common reference T0 = 500 K,
 * rho0 = 15 mol/dm3, the reference state u = s = 0 of each pure saturated liquid at its triple
 * point already included - no other offset), and a residual part (1-x) Phi^r_1 + x Phi^r_2 +
 * Delta Phi^r (Eq. 3) in the mixture reduced variables tau = Tn(x)/T and delta = rho/rhon(x)
 * (Eqs. 4-7), with Phi^r_1 the residual part of IAPWS-95 for water, Phi^r_2 the residual part of
 * the Tillner-Roth & Baehr (1993) equation for ammonia (the guideline Appendix, Table A1) and
 * Delta Phi^r the departure function (Eq. 8, the composition factor x(1 - x^gamma)).
 *
 * The formulation is written once, as a template over the scalar type S (plain `double` or a
 * `taylor::Jet<N>` seeded in (T, rho, x)), so that every derivative is exact (analytical, by
 * truncated Taylor arithmetic); nothing here uses finite differences or a reference library.
 * The phase equilibrium (the saturation solver) and the 15 three-property flashes are built on
 * this core in `nh3h2o_flash.cpp` (declared at the end of this header); the EES fluid functions
 * and the procedure `CALL NH3H2O` (units, derivatives, limits of the trial points) are in
 * `evaluator_nh3h2o.cpp`.
 */

//@{
/// Reference constants of the guideline (Table 1 and Section 4): x is the mole fraction of
/// ammonia, component 1 is water and component 2 is ammonia, all quantities molar.
inline constexpr double kR = 8.314471;       ///< molar gas constant [J/(mol K)]
inline constexpr double kM1 = 18.015268e-3;  ///< molar mass of water [kg/mol]
inline constexpr double kM2 = 17.03026e-3;   ///< molar mass of ammonia [kg/mol]
inline constexpr double kTc1 = 647.096;      ///< critical temperature of water [K]
inline constexpr double kTc2 = 405.40;       ///< critical temperature of ammonia [K]
inline constexpr double kRhoc1 = 322.0 / kM1;  ///< critical molar density of water (322 kg/m3) [mol/m3]
inline constexpr double kRhoc2 = 225.0 / kM2;  ///< critical molar density of ammonia (225 kg/m3) [mol/m3]
inline constexpr double kT0 = 500.0;           ///< reducing temperature of the ideal part [K]
inline constexpr double kRho0 = 15000.0;       ///< reducing molar density of the ideal part (15 mol/dm3) [mol/m3]
//@}

/// The molar mass of the mixture M(x) = (1-x) M1 + x M2 [kg/mol].
inline double molarMass(double x) { return (1.0 - x) * kM1 + x * kM2; }

//@{
/// The pieces of the formulation, each written for S = double or S = Jet<1..3> (exact
/// derivatives through the jet arithmetic). Reduced variables: tau0 = T0/T, delta0 = rho/rho0,
/// tau = Tn(x)/T, delta = rho/rhon(x).

/// The ideal part Phi^o(tau0, delta0, x) of Eq. (2). At x = 0 and x = 1 exactly, the x ln x and
/// (1-x) ln(1-x) terms take their limit 0 as the zero of S (their x derivative is -infinity at
/// the endpoints; every quantity that Table 4 of the guideline derives from Phi needs no x
/// derivative there).
template <class S>
S idealHelmholtz(const S& tau0, const S& delta0, const S& x);

/// The residual part Phi^r_1(tau, delta) of IAPWS-95 for water (Wagner & Pruss 2002; the
/// published coefficients, as embedded in CoolProp's `dev/fluids/Water.json`). At delta = 1
/// exactly the non-analytic terms take the limit of all their existing derivatives (0, see
/// `taylor::pow`); at tau = delta = 1 exactly the term value is 0 and higher tau derivatives,
/// which diverge there, are set to the same limit convention.
template <class S>
S residualWater(const S& tau, const S& delta);

/// The residual part Phi^r_2(tau, delta) of pure ammonia, Eq. (A1) with Table A1 (Tillner-Roth
/// & Baehr 1993; the coefficients cross-checked digit for digit against CoolProp's
/// `dev/fluids/Ammonia.json`, EOS entry `TillnerRoth-DKV-1993`).
template <class S>
S residualAmmonia(const S& tau, const S& delta);

/// The reducing temperature Tn(x) of Eq. (4) [K].
template <class S>
S reducingTemperature(const S& x);

/// The reducing density rhon(x) of Eq. (5) [mol/m3].
template <class S>
S reducingDensity(const S& x);

/// The departure function Delta Phi^r(tau, delta, x) of Eq. (8): the bracket of coefficient
/// terms times the composition factor x(1 - x^gamma).
template <class S>
S departureFunction(const S& tau, const S& delta, const S& x);
//@}

/**
 * @brief The whole dimensionless Helmholtz free energy Phi = Phi^o + Phi^r of Eq. (1).
 *
 * @param T    temperature [K], positive
 * @param rho  molar density [mol/m3], positive
 * @param x    mole fraction of ammonia, in [0, 1]
 * @throws std::runtime_error for a non-finite argument, T <= 0, rho <= 0 or x outside [0, 1].
 */
template <class S>
S reducedHelmholtz(const S& T, const S& rho, const S& x);

/// The thermodynamic properties at (T, rho, x), molar SI units (J, Pa, K, mol, m3), from the
/// exact derivatives of f = R T Phi (the relations of Table 4 of the guideline, evaluated as
/// jets): p [Pa], f, u, h [J/mol], s, cv, cp [J/(mol K)] and the speed of sound w [m/s] with
/// the molar mass M(x).
struct State {
    double p;
    double f;
    double u;
    double h;
    double s;
    double cv;
    double cp;
    double w;
};

/// The properties of `State` at (T, rho, x).
/// @throws std::runtime_error as `reducedHelmholtz` does.
State state(double T, double rho, double x);

/// The quantities p, u, h, s as second-order jets in (T, rho, x) (value, gradient and Hessian
/// of each): they are the `derivative(...)` of the third-order Phi jet, so their own derivatives
/// are exact. Only the tests call it: the saturation solver and the flashes of `nh3h2o_flash.cpp`
/// build their own first-order jets of the phase quantities from the second-order jet of Phi.
struct StateJets {
    taylor::Jet<2> p;
    taylor::Jet<2> u;
    taylor::Jet<2> h;
    taylor::Jet<2> s;
};

/// The jets of `StateJets` at (T, rho, x).
/// @throws std::runtime_error as `reducedHelmholtz` does.
StateJets stateJets(double T, double rho, double x);

/// Mass fraction of ammonia for the mole fraction `x`, and its derivative.
double massFractionOf(double x);
double dMassFractionDx(double x);
/// Mole fraction of ammonia for the mass fraction `xi`, and its derivative.
double moleFractionOf(double xi);

//@{
/**
 * @brief The property designators of the three-property flash, numbered as the outputs of the EES procedure
 *        `NH3H2O(Code, In1, In2, In3 : T, P, x, h, s, u, v, q)`.
 *
 * Every state is fixed by three of them; the EES code is the three property numbers in ascending order
 * (123 = T, P, x; 138 = T, x, q; 234 = P, x, h; ...). `x` is the MASS fraction of ammonia and `q` the vapour
 * mass fraction (0..1 in the two-phase region, -0.001 sub-cooled liquid, 1.001 superheated vapour).
 */
enum class Prop : int { T = 1, P = 2, X = 3, H = 4, S = 5, U = 6, V = 7, Q = 8 };
//@}

/// The result of a flash: the eight state properties in SI units with mass basis (T [K], p [Pa], x [kg/kg], h and
/// u [J/kg], s [J/(kg K)], v [m3/kg], q) and, for each, its partial derivatives with respect to the three input
/// values in the order they were given, exact (implicit function theorem on the converged equations).
struct FlashResult {
    double value[8];           ///< indexed by `static_cast<int>(Prop) - 1`
    double derivative[8][3];   ///< d value[k] / d input[j]
    bool twoPhase = false;     ///< true inside the dome, q in [0, 1]
    double temperature() const { return value[0]; }
    double pressure() const { return value[1]; }
    double massFraction() const { return value[2]; }
    double enthalpy() const { return value[3]; }
    double entropy() const { return value[4]; }
    double internalEnergy() const { return value[5]; }
    double volume() const { return value[6]; }
    double quality() const { return value[7]; }
};

/// A vapour-liquid equilibrium: temperature [K], pressure [Pa], mole fractions of ammonia of the liquid and the vapour
/// and their molar densities [mol/m3].
struct Saturation {
    double T, p, xL, xV, rhoL, rhoV;
};

/// What is given to `saturation`: the bubble point (`a` = T or p, `b` = mole fraction of the liquid), the dew
/// point (`b` = mole fraction of the vapour) or the equilibrium pair at T and p (`a` = T, `b` = p).
enum class SatSpec { BubbleT, DewT, BubbleP, DewP, TP };

/// The vapour-liquid equilibrium of the given specification (SI units).
/// @throws std::runtime_error if there is none (outside the two-phase region, above the critical locus) or the
///         iteration fails.
Saturation saturation(SatSpec spec, double a, double b);

/**
 * @brief The state fixed by three properties (the 15 combinations of the EES procedure `NH3H2O`).
 *
 * `first`, `second`, `third` are three distinct properties with `first < second < third` in the order of
 * `Prop`, and `a`, `b`, `c` their values (SI units, mass basis). The accepted combinations are
 * T-P-x, T-P-q, T-x-v, T-x-q, T-h-q, T-s-q, T-u-q, T-v-q, P-x-h, P-x-s, P-x-q, P-h-q, P-s-q, P-u-q, P-v-q.
 *
 * @throws std::runtime_error naming the state when the combination is not supported, an input is outside its
 *         range (x in [0, 1], q in [0, 1], positive T, p, v), or no state exists / the iteration fails.
 */
FlashResult flash(Prop first, Prop second, Prop third, double a, double b, double c);

/// The specific heat at constant pressure [J/(kg K)] and the volume expansion coefficient (1/v)(dv/dT)_p [1/K] of
/// a single-phase state with their derivatives with respect to the three inputs of `flash` (exact; they need the
/// third derivatives of the Helmholtz energy).
struct SinglePhaseExtras {
    double cp, dcp[3];
    double expansion, dexpansion[3];
};

/// `SinglePhaseExtras` of the state fixed by the same three properties as `flash`.
/// @throws std::runtime_error if the state is two-phase or does not exist.
SinglePhaseExtras singlePhaseExtras(Prop first, Prop second, Prop third, double a, double b, double c);

}  // namespace coolsolve::nh3h2o
