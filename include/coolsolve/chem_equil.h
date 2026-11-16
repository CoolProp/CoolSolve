#pragma once

#include <array>

namespace coolsolve::chem_equil {

/**
 * @file chem_equil.h
 * @brief The chemical equilibrium of the 14 species of `CALL CHEM_EQUIL`: an ideal-gas mixture of H2, O2, H2O, CO, CO2,
 *        OH, H, O, N2, N, NO, NO2, CH4 and Ar at given pressure and temperature, in the gas phase only (no graphite).
 *
 * The mixture is described by the ratios of its atoms to its oxygen atoms (A\\O of argon, C\\O, H\\O, N\\O: the oxygen
 * has the abundance 1); the species follow from the minimum of the Gibbs energy at T and P under the conservation of the
 * five elements. The thermodynamic data are those of the NASA-Glenn polynomials of `nasa_glenn.h` for all 14 species
 * (the enthalpy anchored on the heat of formation, the entropy at 1 bar, R = 8.314462618 J/mol-K).
 *
 * Method: the element potentials of Gordon & McBride, "Computer Program for Calculation of Complex Chemical Equilibrium
 * Compositions and Applications", NASA RP-1311 (1994), sections 2 and 3, as in STANJAN (Reynolds, 1986). At the minimum
 * every mole fraction is `ln x_j = -mu0_j/RT - ln(P/P0) + sum_k a_kj pi_k` (`mu0_j/RT = h_j/RT - s_j/R`, a_kj the atoms
 * of element k in species j, pi_k the potential of the element, P0 = 1 bar), and the unknowns (the potentials of the
 * elements present and nu = ln N, N the moles of the mixture of the abundances above) are fixed by the conservation of
 * the elements and by `sum x_j = 1`. They are found by two nested problems that converge from any start:
 * - at a given nu, the potentials maximise the concave function `D(pi) = sum_k b_k pi_k - sum_j n_j` with
 *   `n_j = exp(nu - mu0_j/RT - ln(P/P0) + sum_k a_kj pi_k)` (the dual of the Gibbs minimum at constant volume): its gradient
 *   is the residual of the element balances, its Hessian `-M = -sum_j n_j a_j a_j^T`; the maximum exists for every nu,
 *   and Newton's method with a line search on D reaches it;
 * - nu solves `psi(nu) = ln(sum_j n_j) - nu = 0`, a decreasing function of slope `-b^T M^-1 b / sum_j n_j` between -1 and
 *   0: Newton's method inside a bisection bracket, the potentials following their tangent `-M^-1 b`.
 * The start is the complete combustion. Where the mixture is held by fewer independent species than there are elements
 * (the stoichiometric mixture at a low temperature: the species that tell an excess of oxygen from a lack of it are
 * 1E-40 of the mixture), M is singular to the rounding: a ridge of 1E-13 on its scaled diagonal bounds that direction,
 * which moves only those trace species. The derivatives with respect to the six inputs are those of the implicit
 * function theorem at the solution, on the bordered matrix `[[M, b], [b^T, 0]]`: exact, no finite difference.
 *
 * Pure C++: no dependency on CoolProp or on the evaluator, no allocation, deterministic (the same inputs give the same
 * result, there is no warm start). An element whose ratio is 0 is absent: the species that contain it are 0 with zero
 * derivatives, including the derivative with respect to that ratio (the function is not differentiable there: its
 * species grow like the ratio; use a small ratio such as 1E-10 when it is an unknown).
 */

/// The number of species of the equilibrium.
constexpr int kSpeciesCount = 14;

/// The number of inputs of the procedure: P, T, A\O, C\O, H\O, N\O in the order of the CALL.
constexpr int kInputCount = 6;

/// The inputs in the order of the CALL, SI: P [Pa], T [K], A\O, C\O, H\O, N\O [-].
using Inputs = std::array<double, kInputCount>;

/// The species in the order of the outputs of the CALL, with the names of the NASA-Glenn table.
inline constexpr const char* kSpeciesNames[kSpeciesCount] = {"H2", "O2", "H2O", "CO", "CO2", "OH", "H",
                                                             "O",  "N2", "N",   "NO", "NO2", "CH4", "Ar"};

/// The names of the inputs as the messages give them, in the order of `Inputs`.
inline constexpr const char* kInputNames[kInputCount] = {"P", "T", "A\\O", "C\\O", "H\\O", "N\\O"};

/// The equilibrium: the mole fractions and their derivatives with respect to the six inputs.
struct Result {
    double x[kSpeciesCount];                 ///< mole fractions [-], sum 1
    double dx[kSpeciesCount][kInputCount];   ///< d x_j / d input i in SI (P [1/Pa], T [1/K], ratios [-]); 0 when not asked
    int iterations;                          ///< the Newton steps taken
};

/// The temperatures at which the data of all the species that can be present have a value [K].
struct TemperatureRange {
    double tMin, tMax;
};

/**
 * @brief The range of T [K] that the NASA-Glenn data of the species present allow: 200 K to 6000 K when hydrogen or
 *        nitrogen is present (H2O and CH4, NO2 stop there), 20000 K otherwise (O, C, Ar). Only the ratios (inputs 2 to 5) are
 *        read.
 */
TemperatureRange temperatureRange(const Inputs& in);

/**
 * @brief The equilibrium mole fractions at `in`.
 *
 * @param in inputs in SI (see `Inputs`)
 * @param derivatives false skips the derivatives (`Result::dx` is then not set; they cost about 10 % of a call, which
 *        takes some 15 microseconds at 2400 K)
 * @throws std::domain_error naming the input when P is not positive, a ratio is negative or not finite, or the carbon
 *         cannot be placed in the 14 species (C\O >= 1 + H\O/4: no equilibrium without graphite)
 * @throws std::out_of_range naming T and its range when T is outside `temperatureRange`
 * @throws std::runtime_error when Newton's method does not converge (not met in the range of the sweeps of the tests)
 */
Result solve(const Inputs& in, bool derivatives = true);

}  // namespace coolsolve::chem_equil
