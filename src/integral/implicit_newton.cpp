/**
 * @file implicit_newton.cpp
 * @brief The Newton iteration shared by the implicit integrators (see implicit_newton.h).
 */
#include "implicit_newton.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>

namespace coolsolve {

namespace {
constexpr double kEps = std::numeric_limits<double>::epsilon();
// The adaptation of the difference steps of scipy's num_jac
const double kDiffSmall = std::pow(kEps, 0.75);   // a difference below this share of |f| is at its rounding level
const double kDiffBig = std::pow(kEps, 0.25);     // above this share, the step is larger than needed
constexpr double kFactorIncrease = 10.0;
constexpr double kFactorDecrease = 0.1;
constexpr double kMinFactor = 1e3 * kEps;
}  // namespace

double weightedRmsNorm(const Eigen::VectorXd& v, const Eigen::VectorXd& scale) {
    if (v.size() == 0) return 0.0;
    return std::sqrt(v.cwiseQuotient(scale).squaredNorm() / static_cast<double>(v.size()));
}

double newtonTolerance(double relTol) {
    return std::max(10.0 * kEps / relTol, std::min(0.03, std::sqrt(relTol)));
}

int RhsJacobian::compute(const RHSFunction& rhs, double t, const Eigen::VectorXd& y, const Eigen::VectorXd& f,
                         double threshold) {
    // Built aside: an evaluation that throws leaves the previous Jacobian and factors as they were
    const Eigen::Index n = y.size();
    Eigen::VectorXd factor = factor_.size() == n ? factor_ : Eigen::VectorXd::Constant(n, std::sqrt(kEps));
    Eigen::MatrixXd J(n, n);
    for (Eigen::Index j = 0; j < n; ++j) {
        // The step in the direction of f (towards where the solution goes), exactly representable
        const double yScale = (f(j) >= 0.0 ? 1.0 : -1.0) * std::max(threshold, std::abs(y(j)));
        double h = (y(j) + factor(j) * yScale) - y(j);
        while (h == 0.0) {
            factor(j) *= kFactorIncrease;
            h = (y(j) + factor(j) * yScale) - y(j);
        }
        Eigen::VectorXd yp = y;
        yp(j) += h;
        const Eigen::VectorXd fp = rhs(t, yp);
        const Eigen::VectorXd diff = fp - f;
        Eigen::Index iMax = 0;
        const double maxDiff = diff.cwiseAbs().maxCoeff(&iMax);
        const double fScale = std::max(std::abs(f(iMax)), std::abs(fp(iMax)));
        J.col(j) = diff / h;
        if (maxDiff < kDiffSmall * fScale) factor(j) *= kFactorIncrease;
        else if (maxDiff > kDiffBig * fScale) factor(j) *= kFactorDecrease;
        factor(j) = std::max(factor(j), kMinFactor);
    }
    J_ = std::move(J);
    factor_ = std::move(factor);
    return static_cast<int>(n);
}

ImplicitStageResult solveImplicitStage(const RHSFunction& rhs, double tNew, const Eigen::VectorXd& yPredict, double c,
                                       const Eigen::VectorXd& psi, const Eigen::PartialPivLU<Eigen::MatrixXd>& lu,
                                       const Eigen::VectorXd& scale, double tol, int maxIterations,
                                       const RenewIterationMatrix& renew) {
    ImplicitStageResult r;
    r.y = yPredict;
    r.d = Eigen::VectorXd::Zero(yPredict.size());
    std::optional<double> dyNormOld;
    for (int k = 0; k < maxIterations; ++k) {
        r.iterations = k + 1;
        Eigen::VectorXd f;
        try {
            f = rhs(tNew, r.y);
        } catch (const std::exception& e) {
            r.failure = e.what();
            return r;
        }
        ++r.evaluations;
        if (!f.allFinite()) {
            r.failure = "the right-hand side is not finite";
            return r;
        }
        const Eigen::VectorXd residual = c * f - psi - r.d;
        r.residual = residual.size() ? residual.lpNorm<Eigen::Infinity>() : 0.0;
        if (renew) {
            try {
                r.evaluations += renew(r.y, f);
            } catch (const std::exception& e) {
                r.failure = e.what();
                return r;
            }
        }
        const Eigen::VectorXd dy = lu.solve(residual);
        if (!dy.allFinite()) {
            r.failure = "the Newton correction is not finite (singular iteration matrix)";
            return r;
        }
        const double dyNorm = weightedRmsNorm(dy, scale);
        std::optional<double> rate;
        if (dyNormOld) rate = dyNorm / *dyNormOld;
        // The rate of the simplified iteration is that of all its iterations: one that will not converge is stopped.
        // The first corrections of the full iteration, made with the Jacobians of points far from the solution, do
        // not tell its rate near it (quadratic): it goes on.
        if (!renew && rate && (*rate >= 1.0 || std::pow(*rate, maxIterations - k) / (1.0 - *rate) * dyNorm > tol)) {
            std::ostringstream why;
            why << "the Newton iteration " << (*rate >= 1.0 ? "diverges" : "converges too slowly")
                << " (rate of convergence " << *rate << ")";
            r.failure = why.str();
            return r;
        }
        r.y += dy;
        r.d += dy;
        if (dyNorm == 0.0 || (rate && *rate < 1.0 && *rate / (1.0 - *rate) * dyNorm < tol)) {
            r.converged = true;
            return r;
        }
        dyNormOld = dyNorm;
    }
    r.failure = "the Newton iteration did not converge in " + std::to_string(maxIterations) + " iterations";
    return r;
}

}  // namespace coolsolve
