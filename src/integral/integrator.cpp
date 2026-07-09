/**
 * @file integrator.cpp
 * @brief Public integrator factory and helper functions.
 *
 * Mirrors the algebraic solver's `createSolver(SolverStrategy)` pattern
 * (`src/solver.cpp`).  Each concrete integrator is constructed by an internal
 * `make*()` function declared in `integrators_internal.h` and implemented in
 * its own translation unit.
 */
#include "coolsolve/integral/integrator.h"
#include "integrators_internal.h"
#include <algorithm>
#include <cctype>
#include <string>

namespace coolsolve {

std::unique_ptr<Integrator> createIntegrator(IntegratorOptions::Method method) {
    switch (method) {
        case IntegratorOptions::EulerExplicit: return makeEulerExplicit();
        case IntegratorOptions::EulerImplicit: return makeEulerImplicit();
        case IntegratorOptions::RK4:           return makeRK4();
        case IntegratorOptions::RK45:          return makeRK45();
    }
    return makeRK4();  // unreachable
}

std::string methodToString(IntegratorOptions::Method m) {
    switch (m) {
        case IntegratorOptions::EulerExplicit: return "EulerExplicit";
        case IntegratorOptions::EulerImplicit: return "EulerImplicit";
        case IntegratorOptions::RK4:           return "RK4";
        case IntegratorOptions::RK45:          return "RK45";
    }
    return "RK4";
}

bool parseIntegralMethod(const std::string& s, IntegratorOptions::Method& out) {
    std::string lower(s.size(), '\0');
    std::transform(s.begin(), s.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "euler" || lower == "eulerexplicit" || lower == "euler_explicit")
        out = IntegratorOptions::EulerExplicit;
    else if (lower == "eulerimplicit" || lower == "euler_implicit" || lower == "implicit")
        out = IntegratorOptions::EulerImplicit;
    else if (lower == "rk4")
        out = IntegratorOptions::RK4;
    else if (lower == "rk45" || lower == "dopri5" || lower == "dormandprince")
        out = IntegratorOptions::RK45;
    else
        return false;
    return true;
}

}  // namespace coolsolve
