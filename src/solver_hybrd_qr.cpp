/**
 * @file solver_hybrd_qr.cpp
 * @brief Implementation of the incremental rank-one QR update
 *        (MINPACK `r1updt` / `r1mpyq` equivalent). See solver_hybrd_qr.h.
 */
#include "coolsolve/solver_hybrd_qr.h"
#include <cmath>

namespace coolsolve {

namespace {

/**
 * @brief Compute a Givens rotation (c, s) and resultant r such that
 *        c*a + s*b = r  and  -s*a + c*b = 0,  c^2 + s^2 = 1.
 *
 * Numerically robust formulation (Golub & Van Loan, Algorithm 5.1.3):
 * avoids over/underflow by dividing through by the larger magnitude term.
 */
void givensRotation(double a, double b, double& c, double& s, double& r) {
    if (b == 0.0) {
        c = (a >= 0.0) ? 1.0 : -1.0;
        s = 0.0;
        r = std::abs(a);
    } else if (a == 0.0) {
        c = 0.0;
        s = (b >= 0.0) ? 1.0 : -1.0;
        r = std::abs(b);
    } else if (std::abs(b) > std::abs(a)) {
        double t = a / b;
        double u = std::copysign(std::sqrt(1.0 + t * t), b);
        s = 1.0 / u;
        c = s * t;
        r = b * u;
    } else {
        double t = b / a;
        double u = std::copysign(std::sqrt(1.0 + t * t), a);
        c = 1.0 / u;
        s = c * t;
        r = a * u;
    }
}

/// Apply the rotation [[c, s], [-s, c]] on the left to rows (i, k) of M
/// (i.e. M(i,:) and M(k,:) are replaced by the rotated rows).
void applyRotationToRows(Eigen::MatrixXd& M, int i, int k, double c, double s) {
    for (Eigen::Index j = 0; j < M.cols(); ++j) {
        double mi = M(i, j), mk = M(k, j);
        M(i, j) = c * mi + s * mk;
        M(k, j) = -s * mi + c * mk;
    }
}

/// Apply the transposed rotation [[c, -s], [s, c]] on the right to columns
/// (i, k) of M (i.e. M(:,i) and M(:,k) are replaced by the rotated columns).
/// This is the accumulation step that keeps Q*R invariant when R is rotated
/// on the left by [[c, s], [-s, c]] on rows (i, k).
void applyRotationToColumns(Eigen::MatrixXd& M, int i, int k, double c, double s) {
    for (Eigen::Index r = 0; r < M.rows(); ++r) {
        double mi = M(r, i), mk = M(r, k);
        M(r, i) = c * mi + s * mk;
        M(r, k) = -s * mi + c * mk;
    }
}

}  // namespace

void computeInitialQR(const Eigen::MatrixXd& A, Eigen::MatrixXd& Q, Eigen::MatrixXd& R) {
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    Q = qr.householderQ();
    R = qr.matrixQR().triangularView<Eigen::Upper>();
}

void rank1QRUpdate(Eigen::MatrixXd& Q, Eigen::MatrixXd& R,
                    const Eigen::VectorXd& u, const Eigen::VectorXd& v) {
    const int n = static_cast<int>(R.rows());
    if (n <= 1) {
        // Trivial 0/1-dimensional case: fall back to direct re-factorization.
        Eigen::MatrixXd A = Q * R + u * v.transpose();
        computeInitialQR(A, Q, R);
        return;
    }

    Eigen::VectorXd w = Q.transpose() * u;

    // Phase 1: bottom-up elimination of w(n-1)..w(1), turning R into upper
    // Hessenberg and w into ||w||*e_0. Rotations accumulated into Q (right).
    for (int k = n - 1; k >= 1; --k) {
        if (w(k) == 0.0) continue;
        double c, s, r;
        givensRotation(w(k - 1), w(k), c, s, r);
        w(k - 1) = r;
        w(k) = 0.0;
        applyRotationToRows(R, k - 1, k, c, s);
        applyRotationToColumns(Q, k - 1, k, c, s);
    }

    // R is now upper Hessenberg. Add w(0) * v^T to row 0.
    R.row(0) += w(0) * v.transpose();

    // Phase 2: top-down elimination of the sub-diagonal, restoring upper
    // triangular R. Rotations accumulated into Q (right).
    for (int k = 0; k < n - 1; ++k) {
        if (R(k + 1, k) == 0.0) continue;
        double c, s, r;
        givensRotation(R(k, k), R(k + 1, k), c, s, r);
        applyRotationToRows(R, k, k + 1, c, s);
        applyRotationToColumns(Q, k, k + 1, c, s);
        R(k, k) = r;
        R(k + 1, k) = 0.0;
    }
}

}  // namespace coolsolve
