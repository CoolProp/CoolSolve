/**
 * @file integrator_euler_implicit.cpp
 * @brief Fixed-step implicit (backward) Euler integrator.
 *
 * Reference: Hairer & Wanner, *Solving ODE II* (Stiff and DAE), §IV.8, and
 * Press et al., *Numerical Recipes* 3rd ed. §17.5.  First-order accurate and
 * A-stable: unconditionally stable for any step size on the linear test
 * equation, which makes it the cheapest viable choice for stiff systems.
 *
 * Each step solves the nonlinear equation
 *     G(y_{n+1}) = y_{n+1} − y_n − h·f(t_{n+1}, y_{n+1}) = 0
 * from the explicit Euler value, with the simplified Newton iteration of the
 * BDF method (implicit_newton.h): a finite-difference Jacobian of f kept from
 * step to step and renewed when the iteration does not converge, a test
 * weighted with relTol / absTol. A fixed step cannot be shortened (BDF can):
 * when even a renewed Jacobian does not let the iteration converge (a kink of
 * f between the predictor and the solution, such as a control that starts to
 * act), the full Newton iteration is tried, with a Jacobian at every iterate.
 * A step that still does not converge is an error (the time and the residual
 * G): accepting it would put a point off the trajectory in the table.
 */
#include "implicit_newton.h"
#include "integrators_internal.h"

#include <sstream>
#include <stdexcept>

namespace coolsolve {

namespace {

/// Iterations of the simplified Newton iteration (those of BDF), then of the full Newton iteration
constexpr int kSimplifiedIterations = 4;
constexpr int kFullIterations = 30;

class EulerImplicitIntegrator : public Integrator {
public:
    StepResult step(double t, const Eigen::VectorXd& y,
                    const RHSFunction& rhs, double h,
                    const IntegratorOptions& opt) override {
        StepResult r;
        const double relTol = opt.relTol > 0.0 ? opt.relTol : 1e-6;
        const double absTol = opt.absTol > 0.0 ? opt.absTol : 1e-9;
        const double t1 = t + h;

        // Predictor: explicit Euler. In the form of implicit_newton.h, c = h and psi = yPredict − y.
        const Eigen::VectorXd f0 = rhs(t, y);
        r.rhsEvaluations = 1;
        const Eigen::VectorXd yPredict = y + h * f0;
        const Eigen::VectorXd psi = yPredict - y;
        const Eigen::VectorXd scale = (relTol * yPredict.cwiseAbs().array() + absTol).matrix();
        if (jacobian_.matrix().rows() != y.size()) {   // the first step: a Jacobian at the initial point
            r.rhsEvaluations += jacobian_.compute(rhs, t, y, f0, absTol);
            luValid_ = false;
        }

        // The Jacobian kept, then renewed at the predicted point, then at every iterate
        const auto factorise = [&] {
            lu_.compute(Eigen::MatrixXd::Identity(y.size(), y.size()) - h * jacobian_.matrix());
            luValid_ = true;
            luStep_ = h;
        };
        const RenewIterationMatrix atEveryIterate = [&](const Eigen::VectorXd& yk, const Eigen::VectorXd& fk) {
            const int evaluations = jacobian_.compute(rhs, t1, yk, fk, absTol);
            factorise();
            return evaluations;
        };
        const double tol = newtonTolerance(relTol);
        ImplicitStageResult stage;
        for (int attempt = 0; attempt < 3; ++attempt) {
            if (attempt == 1) {
                const Eigen::VectorXd f = rhs(t1, yPredict);
                r.rhsEvaluations += 1 + jacobian_.compute(rhs, t1, yPredict, f, absTol);
                luValid_ = false;
            }
            if (!luValid_ || luStep_ != h) factorise();
            stage = attempt < 2 ? solveImplicitStage(rhs, t1, yPredict, h, psi, lu_, scale, tol, kSimplifiedIterations)
                                : solveImplicitStage(rhs, t1, yPredict, h, psi, lu_, scale, tol, kFullIterations,
                                                     atEveryIterate);
            r.rhsEvaluations += stage.evaluations;
            if (stage.converged) break;
        }
        if (!stage.converged) {
            std::ostringstream what;
            what << "EulerImplicit: the iteration of the step from t = " << t << " to " << t1 << " (h = " << h
                 << ") did not converge (" << stage.failure << "): the largest residual of y_new - y - h*f(t+h, y_new)"
                 << " is " << stage.residual << ". Reduce the step of INTEGRAL() (or integralFixedStep), or use BDF "
                 << "(variable step for stiff models), RK4 or RK45";
            throw std::runtime_error(what.str());
        }

        r.yNew = stage.y;
        r.stepTaken = h;
        r.nextStep = h;
        r.errorEstimate = 0.0;
        r.accepted = true;
        return r;
    }
    const char* name() const override { return "EulerImplicit"; }

private:
    RhsJacobian jacobian_;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_;
    bool luValid_ = false;
    double luStep_ = 0.0;   ///< the step h of the factorisation of I − h·J
};

}  // namespace

std::unique_ptr<Integrator> makeEulerImplicit() {
    return std::make_unique<EulerImplicitIntegrator>();
}

}  // namespace coolsolve
