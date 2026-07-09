/**
 * @file integrator_rk45.cpp
 * @brief Variable-step Dormand-Prince RK45 (DOPRI5) embedded integrator.
 *
 * Reference: Dormand & Prince (1980), "A family of embedded Runge-Kutta
 * formulae", J. Comp. Appl. Math. 6(1):19–26; see also Hairer, Nørsett &
 * Wanner, *Solving ODE I* Table II.5.2 ("DOPRI5") and §II.4 on local error
 * estimation and step-size control.
 *
 * Seven stages produce both a 5th-order solution `y5` and a 4th-order solution
 * `y4`; their difference is a free local error estimate.  The step is accepted
 * when the scaled RMS error `err ≤ 1` and the next step size is updated with the
 * classic Hairer controller  h_new = h · fac · err^(−1/5)  (clamped).  DOPRI5 is
 * FSAL (First Same As Last): k7 of an accepted step equals k1 of the next one.
 * This implementation recomputes k1 each call for simplicity of the stateless
 * `step()` interface; the cost is one extra RHS evaluation per step, which is
 * irrelevant for the model-scale algebraic solves that dominate runtime.
 */
#include "integrators_internal.h"
#include <cmath>
#include <limits>

namespace coolsolve {

namespace {
class RK45Integrator : public Integrator {
public:
    StepResult step(double t, const Eigen::VectorXd& y,
                    const RHSFunction& rhs, double h,
                    const IntegratorOptions& opt) override {
        StepResult r;
        const int n = static_cast<int>(y.size());

        // ---- Dormand-Prince Butcher tableau ----
        const double c2 = 1.0 / 5.0,    c3 = 3.0 / 10.0,  c4 = 4.0 / 5.0,
                     c5 = 8.0 / 9.0;
        Eigen::VectorXd k1 = rhs(t, y);
        Eigen::VectorXd k2 = rhs(t + c2 * h, y + h * (1.0 / 5.0) * k1);
        Eigen::VectorXd k3 = rhs(t + c3 * h,
                                 y + h * (3.0 / 40.0 * k1 + 9.0 / 40.0 * k2));
        Eigen::VectorXd k4 = rhs(t + c4 * h,
                                 y + h * (44.0 / 45.0 * k1 - 56.0 / 15.0 * k2
                                          + 32.0 / 9.0 * k3));
        Eigen::VectorXd k5 = rhs(t + c5 * h,
                                 y + h * (19372.0 / 6561.0 * k1
                                          - 25360.0 / 2187.0 * k2
                                          + 64448.0 / 6561.0 * k3
                                          - 212.0 / 729.0 * k4));
        Eigen::VectorXd k6 = rhs(t + h,
                                 y + h * (9017.0 / 3168.0 * k1
                                          - 355.0 / 33.0 * k2
                                          + 46732.0 / 5247.0 * k3
                                          + 49.0 / 176.0 * k4
                                          - 5103.0 / 18656.0 * k5));

        // 5th-order solution (b = a7j row); also the FSAL stage-7 input.
        Eigen::VectorXd y5 = y + h * (35.0 / 384.0 * k1
                                      + 500.0 / 1113.0 * k3
                                      + 125.0 / 192.0 * k4
                                      - 2187.0 / 6784.0 * k5
                                      + 11.0 / 84.0 * k6);
        Eigen::VectorXd k7 = rhs(t + h, y5);

        // 4th-order embedded solution (bhat).
        Eigen::VectorXd y4 = y + h * (5179.0 / 57600.0 * k1
                                      + 7571.0 / 16695.0 * k3
                                      + 393.0 / 640.0 * k4
                                      - 92097.0 / 339200.0 * k5
                                      + 187.0 / 2100.0 * k6
                                      + 1.0 / 40.0 * k7);

        // Error estimate: y5 - y4 = sum(e_i k_i), with e = b - bhat.
        Eigen::VectorXd errVec = h * (71.0 / 57600.0 * k1
                                      - 71.0 / 16695.0 * k3
                                      + 71.0 / 1920.0 * k4
                                      - 17253.0 / 339200.0 * k5
                                      + 22.0 / 525.0 * k6
                                      - 1.0 / 40.0 * k7);

        // ---- Scaled RMS local error and step-size controller (Hairer II.4) ----
        double rel = opt.relTol > 0 ? opt.relTol : 1e-6;
        double abst = opt.absTol > 0 ? opt.absTol : 1e-9;
        double sumsq = 0.0;
        for (int i = 0; i < n; ++i) {
            double sc = abst + rel * std::max(std::abs(y(i)), std::abs(y5(i)));
            if (sc <= 0.0) sc = abst > 0 ? abst : 1e-12;
            double e = errVec(i) / sc;
            sumsq += e * e;
        }
        double err = (n > 0) ? std::sqrt(sumsq / static_cast<double>(n)) : 0.0;

        const double safety = 0.9;
        const double expo = -1.0 / 5.0;  // q = 4 (embedded order) => 1/(q+1)
        const double facMin = 0.2, facMax = 5.0;
        double errSafe = err > 0 ? err : std::numeric_limits<double>::min();
        double factor = safety * std::pow(errSafe, expo);
        if (!std::isfinite(factor)) factor = facMax;
        factor = std::min(facMax, std::max(facMin, factor));

        bool accepted = err <= 1.0;
        r.yNew = accepted ? y5 : y;  // on rejection, keep the previous state
        r.stepTaken = accepted ? h : 0.0;
        r.errorEstimate = err;
        r.accepted = accepted;
        r.nextStep = h * factor;
        r.rhsEvaluations = 7;
        return r;
    }
    const char* name() const override { return "RK45"; }
    int order() const override { return 5; }  // high-order member of the pair
};
}  // namespace

std::unique_ptr<Integrator> makeRK45() {
    return std::make_unique<RK45Integrator>();
}

}  // namespace coolsolve
