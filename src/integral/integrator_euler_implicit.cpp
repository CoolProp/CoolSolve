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
 * explicit Euler value is used as the initial predictor.  A step whose
 * iteration does not converge is an error (the time and the residual G):
 * accepting it would put a point off the trajectory in the table.  In Phase 5 this
 * same coupling is delegated to the algebraic `Solver` when state and
 * algebraic variables are solved together; here the integrator stays
 * self-contained so the numerics can be tested analytically.
 */
#include "integrators_internal.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

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

        bool converged = false;
        double gnorm = 0.0;
        Eigen::VectorXd dy = Eigen::VectorXd::Constant(n, std::numeric_limits<double>::infinity());
        int inner = 0;
        for (; inner < maxInner; ++inner) {
            f1 = rhs(t1, ynew); ++evals;
            G = ynew - y - h * f1;
            gnorm = G.lpNorm<Eigen::Infinity>();
            if (gnorm < tol) { converged = true; break; }

            // Finite-difference Jacobian of f at (t1, ynew).
            for (int j = 0; j < n; ++j) {
                double step = fdEps * (std::abs(ynew(j)) + 1.0);
                Eigen::VectorXd yp = ynew; yp(j) += step;
                fcol = rhs(t1, yp); ++evals;
                Jf.col(j) = (fcol - f1) / step;
            }
            JG = Eigen::MatrixXd::Identity(n, n) - h * Jf;

            auto lu = JG.fullPivLu();
            dy = lu.solve(-G);
            if (!dy.allFinite()) break;
            ynew += dy;
            if (dy.lpNorm<Eigen::Infinity>() < tol) { converged = true; break; }
        }
        // The iterations are spent with a last correction of a few ulps of y: the iteration stalls at the rounding
        // of y (a state above ~1E3 cannot get a correction below 1E-12), the step is converged
        if (!converged && inner == maxInner) {
            converged = true;
            for (int i = 0; i < n; ++i)
                if (!(std::abs(dy(i)) <= 4.0 * std::numeric_limits<double>::epsilon() * std::abs(ynew(i)))) converged = false;
        }
        if (!converged) {
            std::ostringstream what;
            what << "EulerImplicit: the iteration of the step from t = " << t << " to " << t1 << " (h = " << h
                 << ") did not converge after " << std::min(inner + 1, maxInner)
                 << " iterations: the largest residual of y_new - y - h*f(t+h, y_new) is " << gnorm
                 << ". Reduce the step of INTEGRAL() (or integralFixedStep), or use RK4 or RK45";
            throw std::runtime_error(what.str());
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
