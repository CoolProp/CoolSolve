/**
 * @file integrator_euler_implicit.cpp
 * @brief Fixed-step implicit (backward) Euler integrator.
 *
 * Reference: Hairer & Wanner, *Solving ODE II* (Stiff and DAE), §IV.8, and
 * Press et al., *Numerical Recipes* 3rd ed. §17.5.  First-order accurate and
 * A-stable: unconditionally stable for any step size on the linear test
 * equation, which makes it the cheapest viable choice for stiff systems.
 *
 * Each step solves the nonlinear fixed-point equation
 *     G(y_{n+1}) = y_{n+1} − y_n − h·f(t_{n+1}, y_{n+1}) = 0
 * with a small Newton iteration using a finite-difference Jacobian.  The
 * explicit Euler value is used as the initial predictor.  In Phase 5 this
 * same coupling is delegated to the algebraic `Solver` when state and
 * algebraic variables are solved together; here the integrator stays
 * self-contained so the numerics can be tested analytically.
 */
#include "integrators_internal.h"
#include <cmath>

namespace coolsolve {

namespace {
class EulerImplicitIntegrator : public Integrator {
public:
    StepResult step(double t, const Eigen::VectorXd& y,
                    const RHSFunction& rhs, double h,
                    const IntegratorOptions& /*opt*/) override {
        StepResult r;
        const int n = static_cast<int>(y.size());
        const double t1 = t + h;

        // Predictor: explicit Euler.
        Eigen::VectorXd ynew = y + h * rhs(t, y);
        int evals = 1;

        Eigen::VectorXd G(n), f1(n), fcol(n);
        Eigen::MatrixXd Jf(n, n), JG(n, n);
        const int maxInner = 30;
        const double tol = 1e-12;
        const double fdEps = 1e-7;  // finite-difference perturbation scale

        for (int inner = 0; inner < maxInner; ++inner) {
            f1 = rhs(t1, ynew); ++evals;
            G = ynew - y - h * f1;
            double gnorm = G.lpNorm<Eigen::Infinity>();
            if (gnorm < tol) break;

            // Finite-difference Jacobian of f at (t1, ynew).
            for (int j = 0; j < n; ++j) {
                double step = fdEps * (std::abs(ynew(j)) + 1.0);
                Eigen::VectorXd yp = ynew; yp(j) += step;
                fcol = rhs(t1, yp); ++evals;
                Jf.col(j) = (fcol - f1) / step;
            }
            JG = Eigen::MatrixXd::Identity(n, n) - h * Jf;

            auto lu = JG.fullPivLu();
            Eigen::VectorXd dy = lu.solve(-G);
            if (!dy.allFinite()) break;
            ynew += dy;
            if (dy.lpNorm<Eigen::Infinity>() < tol) break;
        }

        r.yNew = ynew;
        r.stepTaken = h;
        r.nextStep = h;
        r.errorEstimate = 0.0;
        r.accepted = true;
        r.rhsEvaluations = evals;
        return r;
    }
    const char* name() const override { return "EulerImplicit"; }
};
}  // namespace

std::unique_ptr<Integrator> makeEulerImplicit() {
    return std::make_unique<EulerImplicitIntegrator>();
}

}  // namespace coolsolve
