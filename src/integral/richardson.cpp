/**
 * @file richardson.cpp
 * @brief Richardson extrapolation wrapper for fixed-step integrators.
 *
 * Reference: Richardson extrapolation as applied to ODE integration, see e.g.
 * Press et al., *Numerical Recipes* §17.3 (step-doubling) and the EES manual
 * page `richardson_extrapolation.htm`.  For a fixed-step method of order `p`,
 * taking one step of size `h` (result `I_h`) and two half-steps of size `h/2`
 * (result `I_{h/2}`) lets us cancel the leading O(h^p) error term:
 *
 *     I ≈ (4 · I_{h/2} − I_h) / 3
 *
 * yielding an effective order `p+1`.  Cost: 3× the base method per step.
 * Only valid for fixed-step methods; the wrapper does not enforce this — it is
 * the caller's responsibility (see `IntegralSolver`).
 */
#include "coolsolve/integral/integrator.h"
#include <memory>

namespace coolsolve {

namespace {
class RichardsonIntegrator : public Integrator {
public:
    explicit RichardsonIntegrator(std::unique_ptr<Integrator> base)
        : base_(std::move(base)) {}

    StepResult step(double t, const Eigen::VectorXd& y,
                    const RHSFunction& rhs, double h,
                    const IntegratorOptions& opt) override {
        StepResult r;
        const double h2 = 0.5 * h;

        // One full step.
        StepResult sFull = base_->step(t, y, rhs, h, opt);
        // Two half-steps.
        StepResult s1 = base_->step(t, y, rhs, h2, opt);
        StepResult s2 = base_->step(t + h2, s1.yNew, rhs, h2, opt);

        // General Richardson combination cancelling the leading O(h^p) term:
        //   I_exact = I_h + C h^p + ...
        //   I_exact = I_{h/2} + C (h/2)^p + ...
        //   =>  I_rich = (2^p · I_{h/2} − I_h) / (2^p − 1)
        // The EES manual's (4·I_{h/2} − I_h)/3 is the p=2 special case.
        const int p = base_->order();
        const double pow2p = std::pow(2.0, p);
        const double wHalf = pow2p / (pow2p - 1.0);
        const double wFull = -1.0 / (pow2p - 1.0);

        r.yNew = wHalf * s2.yNew + wFull * sFull.yNew;
        r.stepTaken = h;
        r.nextStep = h;
        r.errorEstimate = 0.0;
        r.accepted = true;
        r.rhsEvaluations = sFull.rhsEvaluations + s1.rhsEvaluations + s2.rhsEvaluations;
        return r;
    }

    const char* name() const override { return base_ ? base_->name() : "Richardson"; }

private:
    std::unique_ptr<Integrator> base_;
};
}  // namespace

std::unique_ptr<Integrator> wrapRichardson(std::unique_ptr<Integrator> base,
                                           bool enable) {
    if (!enable) return base;
    return std::make_unique<RichardsonIntegrator>(std::move(base));
}

}  // namespace coolsolve
