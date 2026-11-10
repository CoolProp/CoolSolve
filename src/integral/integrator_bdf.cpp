/**
 * @file integrator_bdf.cpp
 * @brief Variable-order, variable-step BDF integrator for stiff systems (the NDF of Shampine and Reichelt).
 *
 * References: Shampine & Reichelt (1997), "The MATLAB ODE Suite", SIAM J. Sci. Comput. 18(1):1-22, §2 (the
 * numerical differentiation formulas NDF1-NDF5, their quasi-constant step size implementation with modified
 * divided differences, and the change of step by the matrices R and U); Hairer & Wanner, *Solving ODE II*, §III.5
 * and §IV.8. The algorithm is that of scipy's `BDF` (scipy/integrate/_ivp/bdf.py), on which the gains of this
 * method were measured (docs/solver_roadmap.md §6.1, step 6).
 *
 * The method keeps, from one step to the next, the backward differences D of the solution (scaled to the current
 * step), a Jacobian of the right-hand side and the factorisation of the Newton matrix. A step:
 * - predicts y from D, then solves the implicit formula with the simplified Newton iteration of implicit_newton.h,
 *   weighted with relTol / absTol; when it does not converge, the Jacobian is renewed at the predicted point and the
 *   iteration tried again; if it still does not converge the step is rejected and halved;
 * - estimates the local error from the Newton correction; above the tolerance the step is rejected and shortened;
 * - once accepted, after order + 1 steps of the same size, compares the errors of the orders k − 1, k and k + 1
 *   and takes the order and step that promise the largest step.
 * The first step is chosen from the right-hand side at the initial point (Hairer, Nørsett & Wanner, *Solving ODE
 * I*, §II.4), within the step proposed by the caller.
 */
#include "implicit_newton.h"
#include "integrators_internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace coolsolve {

namespace {

constexpr int kMaxOrder = 5;
constexpr int kNewtonMaxIterations = 4;
constexpr double kMinFactor = 0.2;    ///< smallest reduction of a step rejected by its error
constexpr double kMaxFactor = 10.0;   ///< largest growth of the step after an accepted one

/// The coefficients of the NDF of order k = 1..5 (index 0 unused): kappa of Shampine and Reichelt, gamma_k = sum 1/j.
struct NdfCoefficients {
    std::array<double, kMaxOrder + 1> gamma{}, alpha{}, errorConst{};
    NdfCoefficients() {
        const double kappa[kMaxOrder + 1] = {0.0, -0.1850, -1.0 / 9.0, -0.0823, -0.0415, 0.0};
        for (int k = 1; k <= kMaxOrder; ++k) gamma[k] = gamma[k - 1] + 1.0 / k;
        for (int k = 0; k <= kMaxOrder; ++k) {
            alpha[k] = (1.0 - kappa[k]) * gamma[k];
            errorConst[k] = kappa[k] * gamma[k] + 1.0 / (k + 1);
        }
    }
};
const NdfCoefficients kNdf;

/// The matrix R(order, factor) of Shampine and Reichelt: R(i, j) = prod_{m=1..i} (m − 1 − factor·j) / m.
Eigen::MatrixXd changeMatrix(int order, double factor) {
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(order + 1, order + 1);
    R.row(0).setOnes();
    for (int i = 1; i <= order; ++i)
        for (int j = 1; j <= order; ++j) R(i, j) = R(i - 1, j) * (i - 1 - factor * j) / i;
    return R;
}

class BDFIntegrator : public Integrator {
public:
    StepResult step(double t, const Eigen::VectorXd& y, const RHSFunction& rhs, double h,
                    const IntegratorOptions& opt) override {
        StepResult r;
        r.yNew = y;
        r.accepted = false;
        relTol_ = opt.relTol > 0.0 ? std::max(opt.relTol, 100.0 * std::numeric_limits<double>::epsilon()) : 1e-6;
        absTol_ = opt.absTol > 0.0 ? opt.absTol : 1e-9;

        // A point that is not the end of the last accepted step starts a new history
        const bool fresh = !started_ || t != t_ || y.size() != yLast_.size() || y != yLast_;
        if (fresh) start(t, y, rhs, h, r);
        const double hUse = fresh ? hAbs_ : h;   // the first step is chosen by the method, within h
        if (hUse != hAbs_) {
            changeDifferences(hUse / hAbs_);
            hAbs_ = hUse;
            nEqualSteps_ = 0;
            luValid_ = false;
        }

        // Predictor and the terms of the formula from the differences
        const int k = order_;
        const double tNew = t + hUse;
        const Eigen::VectorXd yPredict = D_.topRows(k + 1).colwise().sum().transpose();
        const Eigen::VectorXd predictScale = weights(yPredict);
        Eigen::VectorXd psi = Eigen::VectorXd::Zero(y.size());
        for (int j = 1; j <= k; ++j) psi += kNdf.gamma[j] * D_.row(j).transpose();
        psi /= kNdf.alpha[k];
        const double c = hUse / kNdf.alpha[k];

        // The implicit formula; a Jacobian that does not let the iteration converge is renewed once per step
        ImplicitStageResult stage;
        while (true) {
            if (!luValid_) {
                lu_.compute(Eigen::MatrixXd::Identity(y.size(), y.size()) - c * jacobian_.matrix());
                luValid_ = true;
            }
            stage = solveImplicitStage(rhs, tNew, yPredict, c, psi, lu_, predictScale, newtonTol_, kNewtonMaxIterations);
            r.rhsEvaluations += stage.evaluations;
            if (stage.converged || jacobianCurrent_) break;
            try {
                const Eigen::VectorXd f = rhs(tNew, yPredict);
                ++r.rhsEvaluations;
                r.rhsEvaluations += jacobian_.compute(rhs, tNew, yPredict, f, absTol_);
            } catch (const std::exception& e) {
                stage.failure = e.what();
                break;
            }
            jacobianCurrent_ = true;
            luValid_ = false;
        }
        if (!stage.converged) {   // halved at the next call, with a new factorisation
            r.nextStep = 0.5 * hUse;
            r.rejectReason = stage.failure;
            return r;
        }

        // Local error of the step, from the correction d = y − y_predict
        const double safety = 0.9 * (2 * kNewtonMaxIterations + 1) / (2 * kNewtonMaxIterations + stage.iterations);
        const Eigen::VectorXd scale = weights(stage.y);
        const double errorNorm = weightedRmsNorm(kNdf.errorConst[k] * stage.d, scale);
        r.errorEstimate = errorNorm;
        if (errorNorm > 1.0) {
            r.nextStep = hUse * std::max(kMinFactor, safety * std::pow(errorNorm, -1.0 / (k + 1)));
            return r;
        }

        // Accepted: the differences of the new point (D^{j+1} y_n = D^j y_n − D^j y_{n−1}, d = D^{k+1} y_n)
        ++nEqualSteps_;
        t_ = tNew;
        yLast_ = stage.y;
        D_.row(k + 2) = stage.d.transpose() - D_.row(k + 1);
        D_.row(k + 1) = stage.d.transpose();
        for (int i = k; i >= 0; --i) D_.row(i) += D_.row(i + 1);
        jacobianCurrent_ = false;
        hasStep_ = true;
        r.accepted = true;
        r.yNew = stage.y;
        r.stepTaken = hUse;
        r.nextStep = hUse;
        if (nEqualSteps_ < k + 1) return r;

        // The order (k − 1, k or k + 1) and the step that promise the largest step
        const double inf = std::numeric_limits<double>::infinity();
        const double norms[3] = {
            k > 1 ? weightedRmsNorm(kNdf.errorConst[k - 1] * D_.row(k).transpose(), scale) : inf,
            errorNorm,
            k < kMaxOrder ? weightedRmsNorm(kNdf.errorConst[k + 1] * D_.row(k + 2).transpose(), scale) : inf};
        int best = 0;
        double factors[3];
        for (int i = 0; i < 3; ++i) {
            factors[i] = std::pow(norms[i], -1.0 / (k + i));
            if (factors[i] > factors[best]) best = i;
        }
        order_ = k + best - 1;
        const double factor = std::min(kMaxFactor, safety * factors[best]);
        changeDifferences(factor);
        hAbs_ *= factor;
        nEqualSteps_ = 0;
        luValid_ = false;
        r.nextStep = hAbs_;
        return r;
    }

    const char* name() const override { return "BDF"; }

    bool denseOutput(double tq, Eigen::VectorXd& yq) const override {
        if (!hasStep_) return false;
        // The interpolating polynomial of the differences: D_0 + sum_j D_j prod_{m<j} (tq − t + m·h) / ((m + 1)·h)
        yq = D_.row(0).transpose();
        double p = 1.0;
        for (int j = 0; j < order_; ++j) {
            p *= (tq - (t_ - hAbs_ * j)) / (hAbs_ * (j + 1));
            yq += p * D_.row(j + 1).transpose();
        }
        return true;
    }

private:
    /// The weights of the error control and of the Newton test: absTol + relTol·|y|.
    Eigen::VectorXd weights(const Eigen::VectorXd& y) const {
        return (relTol_ * y.cwiseAbs().array() + absTol_).matrix();
    }

    /// A new history at (t, y): order 1, the first step (at most h) and a Jacobian at the initial point.
    void start(double t, const Eigen::VectorXd& y, const RHSFunction& rhs, double h, StepResult& r) {
        const Eigen::Index n = y.size();
        const Eigen::VectorXd f0 = rhs(t, y);
        ++r.rhsEvaluations;
        hAbs_ = std::min(h, initialStep(rhs, t, y, f0, h, r));
        r.rhsEvaluations += jacobian_.compute(rhs, t, y, f0, absTol_);
        D_ = Eigen::MatrixXd::Zero(kMaxOrder + 3, n);
        D_.row(0) = y.transpose();
        D_.row(1) = hAbs_ * f0.transpose();
        order_ = 1;
        nEqualSteps_ = 0;
        newtonTol_ = newtonTolerance(relTol_);
        luValid_ = false;
        jacobianCurrent_ = false;
        hasStep_ = false;
        started_ = true;
        t_ = t;
        yLast_ = y;
    }

    /// The first step: the size at which an explicit step of order 1 would make an error of about 1% of the
    /// tolerance, from f and its change over a trial step (one evaluation; the trial step when it fails).
    double initialStep(const RHSFunction& rhs, double t, const Eigen::VectorXd& y, const Eigen::VectorXd& f0,
                       double hMax, StepResult& r) const {
        const Eigen::VectorXd scale = weights(y);
        const double d0 = weightedRmsNorm(y, scale);
        const double d1 = weightedRmsNorm(f0, scale);
        const double h0 = std::min(hMax, (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1);
        Eigen::VectorXd f1;
        try {
            f1 = rhs(t + h0, y + h0 * f0);
            ++r.rhsEvaluations;
        } catch (const std::exception&) {
            return h0;
        }
        const double d2 = weightedRmsNorm(f1 - f0, scale) / h0;
        const double h1 = (d1 <= 1e-15 && d2 <= 1e-15) ? std::max(1e-6, h0 * 1e-3) : std::sqrt(0.01 / std::max(d1, d2));
        return std::min({100.0 * h0, h1, hMax});
    }

    /// Rescale the differences of order_ to a step factor times the current one (the same polynomial).
    void changeDifferences(double factor) {
        const Eigen::MatrixXd RU = changeMatrix(order_, factor) * changeMatrix(order_, 1.0);
        D_.topRows(order_ + 1) = (RU.transpose() * D_.topRows(order_ + 1)).eval();
    }

    bool started_ = false;
    bool hasStep_ = false;           ///< an accepted step since start: denseOutput is valid
    double relTol_ = 1e-6, absTol_ = 1e-9, newtonTol_ = 1e-3;
    double t_ = 0.0;                 ///< the end of the last accepted step (or the start of the history)
    Eigen::VectorXd yLast_;          ///< the state returned at t_
    Eigen::MatrixXd D_;              ///< the differences, row j = D^j, scaled to the step hAbs_
    double hAbs_ = 0.0;
    int order_ = 1;
    int nEqualSteps_ = 0;            ///< steps taken at the size and order of D_
    RhsJacobian jacobian_;
    bool jacobianCurrent_ = false;   ///< renewed during the current step: not renewed again before it is accepted
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_;
    bool luValid_ = false;           ///< lu_ factorises I − c·J for the current step and order
};

}  // namespace

std::unique_ptr<Integrator> makeBDF() {
    return std::make_unique<BDFIntegrator>();
}

}  // namespace coolsolve
