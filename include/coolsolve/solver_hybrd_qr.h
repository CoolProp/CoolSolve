#pragma once
/**
 * @file solver_hybrd_qr.h
 * @brief Incremental (rank-one) QR factorization update, following the
 *        MINPACK `hybrd` routines `r1updt` / `r1mpyq`.
 *
 * This is standalone numerical infrastructure: it has no dependency on any
 * `NonLinearSolver` and does not change the behaviour of any existing solver.
 * It exists to let a Broyden-style quasi-Newton update be applied directly to
 * a cached QR factorization of the Jacobian approximation, in O(n^2), instead
 * of forming the dense updated matrix and re-factorizing it from scratch in
 * O(n^3) every iteration (see docs/solver_roadmap.md, Tier 2.1).
 *
 * Reference: Golub, G. H., Van Loan, C. F. (2013), "Matrix Computations",
 * 4th ed., §12.5 ("Updating Matrix Factorizations"); Moré, Garbow, Hillstrom
 * (1980), MINPACK `r1updt.f` / `r1mpyq.f".
 */

#include <Eigen/Dense>

namespace coolsolve {

/**
 * @brief Compute the "thin" (square, n x n) QR factorization of a matrix.
 *
 * Convenience wrapper around Eigen::HouseholderQR so that all QR-related
 * code for the hybrd-style solvers lives in one place.
 *
 * @param A  Square input matrix (n x n).
 * @param Q  [out] Orthogonal factor (n x n), Q*Q^T = I.
 * @param R  [out] Upper-triangular factor (n x n), A = Q*R.
 */
void computeInitialQR(const Eigen::MatrixXd& A, Eigen::MatrixXd& Q, Eigen::MatrixXd& R);

/**
 * @brief Update an existing QR factorization for a rank-one change to the
 *        underlying matrix: A' = A + u*v^T, where A = Q*R.
 *
 * On return, Q and R are overwritten in place so that Q*R = A + u*v^T, with
 * Q orthogonal and R upper triangular. Cost is O(n^2), versus O(n^3) for
 * re-factorizing A' from scratch.
 *
 * Algorithm (Golub & Van Loan §12.5):
 *   1. w = Q^T * u.
 *   2. Apply n-1 Givens rotations (bottom-up) to zero w(n-1)..w(1), turning
 *      R into upper Hessenberg H and w into ||w||*e_0. The same rotations are
 *      applied to Q (accumulated on the right).
 *   3. Add ||w|| * v^T to row 0 of H (still upper Hessenberg).
 *   4. Apply n-1 Givens rotations (top-down) to eliminate the sub-diagonal,
 *      restoring upper-triangular R. Same rotations accumulated into Q.
 *
 * @param Q  [in,out] Orthogonal factor, n x n.
 * @param R  [in,out] Upper-triangular factor, n x n.
 * @param u  Left update vector (n).
 * @param v  Right update vector (n).
 */
void rank1QRUpdate(Eigen::MatrixXd& Q, Eigen::MatrixXd& R,
                    const Eigen::VectorXd& u, const Eigen::VectorXd& v);

}  // namespace coolsolve
