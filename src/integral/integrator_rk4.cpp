/**
 * @file integrator_rk4.cpp
 * @brief Fixed-step classic 4th-order Runge-Kutta integrator (RK4).
 *
 * Reference: the classical fourth-order Runge-Kutta method (Kutta 1901);
 * see Hairer, Nørsett & Wanner, *Solving ODE I* Table II.1.1, and
 * Press et al., *Numerical Recipes* 3rd ed. §17.1.  Four RHS evaluations
 * per step, fourth-order global accuracy.  This is the recommended default
 * integrator for non-stiff models (see `integralMethod = RK4`).
 */
#include "integrators_internal.h"

namespace coolsolve {

namespace {
class RK4Integrator : public Integrator {
public:
    StepResult step(double t, const Eigen::VectorXd& y,
                    const RHSFunction& rhs, double h,
                    const IntegratorOptions& /*opt*/) override {
        StepResult r;
        const double h2 = 0.5 * h;
        const double h6 = h / 6.0;

        Eigen::VectorXd k1 = rhs(t, y);
        Eigen::VectorXd k2 = rhs(t + h2, y + h2 * k1);
        Eigen::VectorXd k3 = rhs(t + h2, y + h2 * k2);
        Eigen::VectorXd k4 = rhs(t + h,   y + h   * k3);

        r.yNew = y + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        r.stepTaken = h;
        r.nextStep = h;
        r.errorEstimate = 0.0;
        r.accepted = true;
        r.rhsEvaluations = 4;
        return r;
    }
    const char* name() const override { return "RK4"; }
    int order() const override { return 4; }
};
}  // namespace

std::unique_ptr<Integrator> makeRK4() {
    return std::make_unique<RK4Integrator>();
}

}  // namespace coolsolve
