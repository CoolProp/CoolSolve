/**
 * @file test_hybrd_qr.cpp
 * @brief Correctness tests for the incremental rank-one QR update
 *        (`rank1QRUpdate`, MINPACK `r1updt`/`r1mpyq` equivalent).
 *
 * This is standalone numerical infrastructure (docs/solver_roadmap.md,
 * Tier 2.1): it is tested here in complete isolation from any
 * `NonLinearSolver`. Correctness is checked via two basis-independent
 * invariants rather than comparing R entries directly, since a QR
 * factorization is only unique up to the sign of each diagonal entry of R
 * (and the corresponding column of Q):
 *   1. Q remains orthogonal (Q^T Q = I).
 *   2. Q*R reconstructs the updated matrix A + u*v^T to near machine
 *      precision, tracked independently as a dense matrix.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "coolsolve/solver_hybrd_qr.h"
#include <random>

using namespace coolsolve;
using Catch::Matchers::WithinAbs;

namespace {

/// Deterministic PRNG per test run (not std::random_device) so failures are
/// reproducible.
std::mt19937& rng() {
    static std::mt19937 gen(20260701);
    return gen;
}

Eigen::MatrixXd randomMatrix(int n, double scale = 1.0) {
    std::uniform_real_distribution<double> dist(-scale, scale);
    Eigen::MatrixXd M(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M(i, j) = dist(rng());
    return M;
}

Eigen::VectorXd randomVector(int n, double scale = 1.0) {
    std::uniform_real_distribution<double> dist(-scale, scale);
    Eigen::VectorXd v(n);
    for (int i = 0; i < n; ++i) v(i) = dist(rng());
    return v;
}

}  // namespace

TEST_CASE("rank1QRUpdate - single update reconstructs A + u*v^T", "[hybrd][qr]") {
    for (int n : {2, 3, 5, 10, 30}) {
        Eigen::MatrixXd A = randomMatrix(n);
        // Ensure A is well-conditioned enough to factor (avoid singular draws).
        A += n * Eigen::MatrixXd::Identity(n, n);

        Eigen::MatrixXd Q, R;
        computeInitialQR(A, Q, R);

        // Sanity: Q*R must reconstruct A right after factorization.
        CHECK((Q * R - A).lpNorm<Eigen::Infinity>() < 1e-9);
        CHECK((Q.transpose() * Q - Eigen::MatrixXd::Identity(n, n)).lpNorm<Eigen::Infinity>() < 1e-9);

        Eigen::VectorXd u = randomVector(n);
        Eigen::VectorXd v = randomVector(n);
        Eigen::MatrixXd expected = A + u * v.transpose();

        rank1QRUpdate(Q, R, u, v);

        INFO("n=" << n);
        CHECK((Q * R - expected).lpNorm<Eigen::Infinity>() < 1e-8);
        CHECK((Q.transpose() * Q - Eigen::MatrixXd::Identity(n, n)).lpNorm<Eigen::Infinity>() < 1e-8);
        // R must be upper triangular (strictly lower part ~ 0).
        for (int i = 1; i < n; ++i)
            for (int j = 0; j < i; ++j)
                CHECK(std::abs(R(i, j)) < 1e-8);
    }
}

TEST_CASE("rank1QRUpdate - 10 sequential updates stay accurate", "[hybrd][qr]") {
    // Per docs/solver_roadmap.md §4.1: 50 random matrices x n in {2,5,10,30},
    // 10 sequential rank-1 updates each, comparing against a from-scratch QR
    // of the accumulated dense matrix (via reconstruction + orthogonality,
    // which is basis/sign-independent — see file header).
    for (int n : {2, 5, 10, 30}) {
        for (int trial = 0; trial < 50; ++trial) {
            Eigen::MatrixXd A = randomMatrix(n) + n * Eigen::MatrixXd::Identity(n, n);
            Eigen::MatrixXd Q, R;
            computeInitialQR(A, Q, R);

            Eigen::MatrixXd tracked = A;
            for (int step = 0; step < 10; ++step) {
                Eigen::VectorXd u = randomVector(n, 0.1);
                Eigen::VectorXd v = randomVector(n, 0.1);
                rank1QRUpdate(Q, R, u, v);
                tracked += u * v.transpose();
            }

            double scale = std::max(1.0, tracked.lpNorm<Eigen::Infinity>());
            INFO("n=" << n << " trial=" << trial);
            CHECK((Q * R - tracked).lpNorm<Eigen::Infinity>() < 1e-6 * scale);
            CHECK((Q.transpose() * Q - Eigen::MatrixXd::Identity(n, n)).lpNorm<Eigen::Infinity>() < 1e-8);
        }
    }
}

TEST_CASE("rank1QRUpdate - identity matrix, simple update", "[hybrd][qr]") {
    // Deterministic sanity check independent of randomness: A = I, update
    // with u = e_0, v = e_0 gives A' = I + e_0 e_0^T = diag(2, 1, 1).
    const int n = 3;
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd Q, R;
    computeInitialQR(A, Q, R);

    Eigen::VectorXd u = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
    u(0) = 1.0;
    v(0) = 1.0;

    rank1QRUpdate(Q, R, u, v);

    Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(n, n);
    expected(0, 0) = 2.0;
    CHECK((Q * R - expected).lpNorm<Eigen::Infinity>() < 1e-12);
}

TEST_CASE("rank1QRUpdate - n=1 edge case", "[hybrd][qr]") {
    Eigen::MatrixXd A(1, 1);
    A(0, 0) = 4.0;
    Eigen::MatrixXd Q, R;
    computeInitialQR(A, Q, R);

    Eigen::VectorXd u(1), v(1);
    u(0) = 2.0;
    v(0) = 3.0;
    rank1QRUpdate(Q, R, u, v);

    CHECK_THAT((Q * R)(0, 0), WithinAbs(10.0, 1e-9));  // 4 + 2*3 = 10
}
