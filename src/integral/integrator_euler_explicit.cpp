/**
 * @file integrator_euler_explicit.cpp
 * @brief Fixed-step explicit (forward) Euler integrator.
 *
 * Reference: Hairer, Nørsett & Wanner, *Solving ODE I* §II.1, and
 * Press et al., *Numerical Recipes* 3rd ed. §17.1.  First-order accurate,
 * cheapest possible step (one RHS evaluation).  Conditionally stable:
 * stable for |1 + h λ| ≤ 1 on the linear test equation, so it is a poor
 * choice for stiff systems (use EulerImplicit or RK45 there).
 */
#include "integrators_internal.h"

namespace coolsolve {

namespace {
class EulerExplicitIntegrator : public Integrator {
public:
    StepResult step(double t, const Eigen::VectorXd& y,
                    const RHSFunction& rhs, double h,
                    const IntegratorOptions& /*opt*/) override {
        StepResult r;
        Eigen::VectorXd f0 = rhs(t, y);
        r.yNew = y + h * f0;
        r.stepTaken = h;
        r.nextStep = h;
        r.errorEstimate = 0.0;
        r.accepted = true;
        r.rhsEvaluations = 1;
        return r;
    }
    const char* name() const override { return "EulerExplicit"; }
};
}  // namespace

std::unique_ptr<Integrator> makeEulerExplicit() {
    return std::make_unique<EulerExplicitIntegrator>();
}

}  // namespace coolsolve
