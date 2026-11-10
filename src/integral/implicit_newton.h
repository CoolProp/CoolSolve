#pragma once
/**
 * @file implicit_newton.h
 * @brief The Newton iteration shared by the implicit integrators (BDF, EulerImplicit).
 *
 * An implicit step solves y = y_predict + d with
 *     c·f(t_new, y) − psi − d = 0,
 * c a multiple of the step and psi a combination of the past values (for implicit Euler with the explicit Euler
 * predictor, c = h and psi = y_predict − y_n). The iteration is a simplified Newton iteration: the matrix I − c·J
 * is factorised with a Jacobian J of f kept from step to step, and renewed by the integrator when the iteration
 * does not converge. Its test is weighted like the error control: scale_i = absTol + relTol·|y_i|, and the
 * iteration stops when the error left, estimated from the rate of convergence, is below a fraction of it.
 *
 * Reference: Hairer & Wanner, *Solving ODE II*, §IV.8; Shampine & Reichelt (1997), "The MATLAB ODE Suite",
 * SIAM J. Sci. Comput. 18(1); the iteration of scipy's `BDF` (`solve_bdf_system`, `num_jac`).
 */
#include "coolsolve/integral/integrator.h"

#include <Eigen/Dense>
#include <functional>
#include <string>

namespace coolsolve {

/// sqrt(mean((v_i / scale_i)^2)): the norm of the error control and of the Newton test (0 for an empty vector).
double weightedRmsNorm(const Eigen::VectorXd& v, const Eigen::VectorXd& scale);

/// Tolerance of the Newton test on the weighted norm for the relative tolerance \p relTol:
/// max(10·eps / relTol, min(0.03, sqrt(relTol))), as scipy's BDF.
double newtonTolerance(double relTol);

/**
 * @brief Forward-difference Jacobian of the right-hand side, n evaluations of f.
 *
 * The perturbation of y_j is factor_j·max(threshold, |y_j|), in the direction of f_j. Each factor starts at
 * sqrt(eps) and is adapted from one Jacobian to the next on the size of the difference it produced (raised when
 * the difference is at the rounding level of f, lowered when it is large), as scipy's `num_jac` does, without its
 * second evaluation of the columns whose difference is too small.
 */
class RhsJacobian {
public:
    /// J at (t, y), with f = rhs(t, y) already evaluated; returns the number of evaluations of rhs (n).
    int compute(const RHSFunction& rhs, double t, const Eigen::VectorXd& y, const Eigen::VectorXd& f, double threshold);
    const Eigen::MatrixXd& matrix() const { return J_; }

private:
    Eigen::MatrixXd J_;
    Eigen::VectorXd factor_;
};

/// What an implicit stage gave.
struct ImplicitStageResult {
    bool converged = false;
    int iterations = 0;          ///< Newton iterations made
    int evaluations = 0;         ///< evaluations of the right-hand side
    Eigen::VectorXd y;           ///< the last iterate
    Eigen::VectorXd d;           ///< y − y_predict
    double residual = 0.0;       ///< largest |c·f − psi − d| at the last evaluation (the residual of the stage)
    std::string failure;         ///< why the iteration stopped without converging (empty when converged)
};

/// Renews the factorisation used by solveImplicitStage at the iterate y, where f = f(t_new, y); returns the
/// evaluations of the right-hand side it made. Throws when it cannot evaluate them.
using RenewIterationMatrix = std::function<int(const Eigen::VectorXd& y, const Eigen::VectorXd& f)>;

/**
 * @brief The simplified Newton iteration of an implicit stage (see the file comment), from y_predict.
 *
 * Converged when the weighted norm of the correction is 0, or when rate/(1 − rate)·‖dy‖ < \p tol (rate: the ratio
 * of two successive corrections). It stops without converging when the rate is 1 or more, when the error predicted
 * after the iterations left exceeds \p tol, after \p maxIterations, or when f is not finite or cannot be evaluated
 * (an exception of \p rhs, whose text is kept in `failure`).
 * @param lu  The factorisation of I − c·J
 * @param scale  The weights absTol + relTol·|y_predict|
 * @param renew  When given, renews \p lu at every iterate before its correction: the full Newton iteration, which
 *               stops on its rate only to converge (it runs its \p maxIterations otherwise)
 */
ImplicitStageResult solveImplicitStage(const RHSFunction& rhs, double tNew, const Eigen::VectorXd& yPredict, double c,
                                       const Eigen::VectorXd& psi, const Eigen::PartialPivLU<Eigen::MatrixXd>& lu,
                                       const Eigen::VectorXd& scale, double tol, int maxIterations,
                                       const RenewIterationMatrix& renew = nullptr);

}  // namespace coolsolve
