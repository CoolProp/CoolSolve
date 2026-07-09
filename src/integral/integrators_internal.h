#pragma once
/**
 * @file integrators_internal.h
 * @brief Internal factory functions for each concrete integrator.
 *
 * The concrete integrator classes are defined entirely inside their
 * `src/integral/integrator_*.cpp` translation units (they are not part of the
 * public API).  This header lets the public `createIntegrator()` factory and
 * the Richardson wrapper refer to them without exposing class definitions.
 */
#include "coolsolve/integral/integrator.h"
#include <memory>

namespace coolsolve {

std::unique_ptr<Integrator> makeEulerExplicit();
std::unique_ptr<Integrator> makeEulerImplicit();
std::unique_ptr<Integrator> makeRK4();
std::unique_ptr<Integrator> makeRK45();

}  // namespace coolsolve
