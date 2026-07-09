#pragma once
/**
 * @file integrator_test_util.h
 * @brief Small march-loop helpers shared by the integrator unit tests.
 *
 * These replicate the outer time-marching loop that, in Phase 5, becomes the
 * responsibility of `IntegralSolver`.  Kept here so the pure-numerics tests
 * (Phase 0) can drive the `Integrator::step()` interface against analytical
 * solutions without any IR/evaluator coupling.
 */
#include "coolsolve/integral/integrator.h"
#include <catch2/catch_test_macros.hpp>
#include <cmath>

namespace coolsolve::test {

/// Fixed-step march: take `nSteps` steps of size `(tf-t0)/nSteps`.
inline Eigen::VectorXd marchFixed(Integrator& ig, double t0, double tf,
                                  const Eigen::VectorXd& y0,
                                  const RHSFunction& rhs, int nSteps,
                                  const IntegratorOptions& opt) {
    const double h = (tf - t0) / nSteps;
    double t = t0;
    Eigen::VectorXd y = y0;
    for (int i = 0; i < nSteps; ++i) {
        StepResult r = ig.step(t, y, rhs, h, opt);
        y = r.yNew;
        t += r.stepTaken;
    }
    return y;
}

/**
 * @brief Adaptive (RK45) march from `t0` to `tf`.
 *
 * Starts from a proposed `h0` step, accepts steps that satisfy the tolerance,
 * shrinks rejected steps, and grows accepted ones.  Honours `maxSteps`,
 * `minStep`, and `maxStep`.  Lands exactly on `tf` by clamping the final step.
 */
inline Eigen::VectorXd marchAdaptive(Integrator& ig, double t0, double tf,
                                     const Eigen::VectorXd& y0,
                                     const RHSFunction& rhs, double h0,
                                     const IntegratorOptions& opt,
                                     int* totalSteps = nullptr,
                                     int* rejectedSteps = nullptr) {
    const double span = tf - t0;
    const double hMin = opt.minStep > 0 ? opt.minStep : span * 1e-9;
    const double hMax = opt.maxStep > 0 ? opt.maxStep : span;
    double t = t0;
    double h = std::min(h0, hMax);
    Eigen::VectorXd y = y0;
    int taken = 0, rejected = 0;
    const int budget = opt.maxSteps > 0 ? opt.maxSteps * 50 : 100000;
    for (int i = 0; i < budget && t < tf - 1e-13; ++i) {
        if (t + h > tf) h = tf - t;  // land exactly on the endpoint
        StepResult r = ig.step(t, y, rhs, h, opt);
        if (r.accepted) {
            y = r.yNew;
            t += r.stepTaken;
            ++taken;
            // Adopt the integrator's recommended next step (clamped to bounds).
            h = std::min(hMax, std::max(hMin, r.nextStep > 0 ? r.nextStep : h));
        } else {
            ++rejected;
            h = std::max(hMin, r.nextStep > 0 ? r.nextStep : h * 0.5);
        }
    }
    if (totalSteps) *totalSteps = taken;
    if (rejectedSteps) *rejectedSteps = rejected;
    return y;
}

}  // namespace coolsolve::test
