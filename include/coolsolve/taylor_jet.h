#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

namespace coolsolve::taylor {

/**
 * @file taylor_jet.h
 * @brief A truncated multivariate Taylor polynomial in exactly 3 independent variables, the
 *        differentiation workhorse of the modules that need exact higher derivatives of a
 *        thermodynamic potential (e.g. `nh3h2o`).
 *
 * A `Jet<Order>` holds the Taylor coefficients of a function f around a point, in the increments
 * h = (h0, h1, h2) of the three variables: f = sum c(i,j,k) h0^i h1^j h2^k over the multi-indices
 * with i + j + k <= Order. In Taylor form (c(i,j,k) = the derivative divided by the factorial of
 * the multi-index) the product of two truncated polynomials is the plain polynomial product, done
 * with a precomputed index table (`kProductTable`), and everything else follows from the
 * one-variable chain rule: f(a) with f given through its derivatives up to third order at the
 * value of `a` (see `chain`). The derivative accessors convert back to ordinary derivatives.
 */

/// The number of Taylor coefficients of a jet of the given order (multi-indices i+j+k <= order).
constexpr int jetSize(int order) { return (order + 1) * (order + 2) * (order + 3) / 6; }

/// The linear storage index of the coefficient of h0^i h1^j h2^k: blocks of equal total degree,
/// lexicographic within a block. The first jetSize(Order) indices of any higher-order jet are
/// those of its truncation, which is what the truncating constructor of Jet copies.
constexpr int jetIndex(int i, int j, int k) {
    const int d = i + j + k;
    return d * (d + 1) * (d + 2) / 6 + i * (2 * d - i + 3) / 2 + j;
}

/// One entry of the product table: r.c[c] += a.c[a] * b.c[b].
struct ProductTerm {
    unsigned char a;
    unsigned char b;
    unsigned char c;
};

/// The number of index pairs whose product lands inside an order-`Order` jet.
constexpr int productTermCount(int order) {
    int count = 0;
    for (int i = 0; i <= order; ++i)
        for (int j = 0; i + j <= order; ++j)
            for (int k = 0; i + j + k <= order; ++k)
                for (int i2 = 0; i + i2 <= order; ++i2)
                    for (int j2 = 0; i + i2 + j + j2 <= order; ++j2)
                        for (int k2 = 0; i + i2 + j + j2 + k + k2 <= order; ++k2) ++count;
    return count;
}

/// The product table of a given order, built at compile time.
template <int Order>
constexpr std::array<ProductTerm, productTermCount(Order)> makeProductTable() {
    std::array<ProductTerm, productTermCount(Order)> table{};
    std::size_t n = 0;
    for (int i = 0; i <= Order; ++i)
        for (int j = 0; i + j <= Order; ++j)
            for (int k = 0; i + j + k <= Order; ++k)
                for (int i2 = 0; i + i2 <= Order; ++i2)
                    for (int j2 = 0; i + i2 + j + j2 <= Order; ++j2)
                        for (int k2 = 0; i + i2 + j + j2 + k + k2 <= Order; ++k2)
                            table[n++] = {static_cast<unsigned char>(jetIndex(i, j, k)),
                                          static_cast<unsigned char>(jetIndex(i2, j2, k2)),
                                          static_cast<unsigned char>(jetIndex(i + i2, j + j2, k + k2))};
    return table;
}

template <int Order>
inline constexpr std::array<ProductTerm, productTermCount(Order)> kProductTable = makeProductTable<Order>();

/// A truncated Taylor polynomial in 3 variables, `Order` in 0..3 (value, 3 first, 6 second and
/// 10 third Taylor coefficients). All arithmetic keeps the Taylor coefficients; the accessors
/// `gradient`, `hessian` and `thirdDerivative` return ordinary (factorial-scaled) derivatives.
template <int Order>
struct Jet {
    static_assert(Order >= 0, "the jet order is 0..3");

    /// Taylor coefficients c(i,j,k) in the storage order of `jetIndex`; default is the zero jet.
    std::array<double, jetSize(Order)> c{};

    Jet() = default;

    /// The constant `value`.
    static Jet constant(double value) {
        Jet jet;
        jet.c[0] = value;
        return jet;
    }

    /// The independent variable number `index` (0, 1 or 2) at the value `value`.
    static Jet variable(double value, int index) {
        Jet jet = constant(value);
        jet.c[jetIndex(index == 0 ? 1 : 0, index == 1 ? 1 : 0, index == 2 ? 1 : 0)] = 1.0;
        return jet;
    }

    /// The truncation of a higher-order jet to this order (drops the high-degree coefficients).
    template <int Higher, class = std::enable_if_t<(Higher >= Order)>>
    explicit Jet(const Jet<Higher>& jet) {
        for (int i = 0; i < jetSize(Order); ++i) c[i] = jet.c[i];
    }

    /// The value (the coefficient of h^0).
    double value() const { return c[0]; }

    /// The derivative with respect to the variable number `var` (0, 1 or 2).
    double gradient(int var) const { return c[jetIndex(var == 0 ? 1 : 0, var == 1 ? 1 : 0, var == 2 ? 1 : 0)]; }

    /// The second derivative with respect to the variables `i` and `j`.
    double hessian(int i, int j) const { return derivativeOf(i, j, -1); }

    /// The third derivative with respect to the variables `i`, `j` and `k`.
    double thirdDerivative(int i, int j, int k) const { return derivativeOf(i, j, k); }

private:
    /// The derivative with respect to up to three variables (`k` = -1 for two): the Taylor
    /// coefficient times the factorials of the powers of the multi-index.
    double derivativeOf(int i, int j, int k) const {
        const int e0 = (i == 0) + (j == 0) + (k == 0), e1 = (i == 1) + (j == 1) + (k == 1),
                  e2 = (i == 2) + (j == 2) + (k == 2);
        constexpr double factorial[4] = {1.0, 1.0, 2.0, 6.0};
        return factorial[e0] * factorial[e1] * factorial[e2] * c[jetIndex(e0, e1, e2)];
    }
};

template <int Order>
Jet<Order> operator+(const Jet<Order>& a, const Jet<Order>& b) {
    Jet<Order> r;
    for (int i = 0; i < jetSize(Order); ++i) r.c[i] = a.c[i] + b.c[i];
    return r;
}

template <int Order>
Jet<Order> operator-(const Jet<Order>& a, const Jet<Order>& b) {
    Jet<Order> r;
    for (int i = 0; i < jetSize(Order); ++i) r.c[i] = a.c[i] - b.c[i];
    return r;
}

template <int Order>
Jet<Order> operator-(const Jet<Order>& a) {
    Jet<Order> r;
    for (int i = 0; i < jetSize(Order); ++i) r.c[i] = -a.c[i];
    return r;
}

template <int Order>
Jet<Order> operator*(const Jet<Order>& a, const Jet<Order>& b) {
    Jet<Order> r;
    for (const ProductTerm& term : kProductTable<Order>) r.c[term.c] += a.c[term.a] * b.c[term.b];
    return r;
}

/// The reciprocal of `b`, by the chain rule for f(x) = 1/x (the value of `b` must not be 0).
template <int Order>
Jet<Order> inverse(const Jet<Order>& b) {
    const double b0 = b.value();
    return chain(b, 1.0 / b0, -1.0 / (b0 * b0), 2.0 / (b0 * b0 * b0), -6.0 / (b0 * b0 * b0 * b0));
}

template <int Order>
Jet<Order> operator/(const Jet<Order>& a, const Jet<Order>& b) {
    return a * inverse(b);
}

template <int Order>
Jet<Order> operator+(const Jet<Order>& a, double b) {
    Jet<Order> r = a;
    r.c[0] += b;
    return r;
}

template <int Order>
Jet<Order> operator+(double a, const Jet<Order>& b) {
    return b + a;
}

template <int Order>
Jet<Order> operator-(const Jet<Order>& a, double b) {
    Jet<Order> r = a;
    r.c[0] -= b;
    return r;
}

template <int Order>
Jet<Order> operator-(double a, const Jet<Order>& b) {
    return -b + a;
}

template <int Order>
Jet<Order> operator*(const Jet<Order>& a, double b) {
    Jet<Order> r;
    for (int i = 0; i < jetSize(Order); ++i) r.c[i] = a.c[i] * b;
    return r;
}

template <int Order>
Jet<Order> operator*(double a, const Jet<Order>& b) {
    return b * a;
}

template <int Order>
Jet<Order> operator/(const Jet<Order>& a, double b) {
    return a * (1.0 / b);
}

template <int Order>
Jet<Order> operator/(double a, const Jet<Order>& b) {
    return a * inverse(b);
}

/// Comparison of the value part only (exact; for tests and domain guards).
template <int Order>
bool operator==(const Jet<Order>& a, const Jet<Order>& b) {
    return a.value() == b.value();
}

template <int Order>
bool operator==(const Jet<Order>& a, double b) {
    return a.value() == b;
}

template <int Order>
bool operator==(double a, const Jet<Order>& b) {
    return a == b.value();
}

template <int Order>
bool operator!=(const Jet<Order>& a, const Jet<Order>& b) {
    return !(a == b);
}

template <int Order>
bool operator!=(const Jet<Order>& a, double b) {
    return !(a == b);
}

template <int Order>
bool operator!=(double a, const Jet<Order>& b) {
    return !(a == b);
}

/// The partial derivative with respect to the variable number `var`, as a jet of one order less.
template <int Order>
Jet<Order - 1> derivative(const Jet<Order>& a, int var) {
    static_assert(Order >= 1, "derivative needs a jet of order >= 1");
    Jet<Order - 1> r;
    for (int i = 0; i + 1 <= Order; ++i)
        for (int j = 0; i + j + 1 <= Order; ++j)
            for (int k = 0; i + j + k + 1 <= Order; ++k) {
                const int factor = var == 0 ? i + 1 : var == 1 ? j + 1 : k + 1;
                r.c[jetIndex(i, j, k)] =
                    factor * a.c[jetIndex(var == 0 ? i + 1 : i, var == 1 ? j + 1 : j, var == 2 ? k + 1 : k)];
            }
    return r;
}

/**
 * @brief The one-variable chain rule of order 3: f(a) = f0 + f1 h + (f2/2) h^2 + (f3/6) h^3 with
 *        h = a - value(a) and f0..f3 the derivatives of f at value(a).
 *
 * The powers of h are jet products, so the truncation is the polynomial one; higher-order f
 * derivatives would be needed only for orders beyond 3.
 */
template <int Order>
Jet<Order> chain(const Jet<Order>& a, double f0, double f1, double f2, double f3) {
    Jet<Order> h = a;
    h.c[0] = 0.0;
    Jet<Order> r = Jet<Order>::constant(f0);
    if (Order >= 1) r = r + h * f1;
    if (Order >= 2) {
        const Jet<Order> h2 = h * h;
        r = r + h2 * (f2 / 2.0);
        if (Order >= 3) r = r + (h2 * h) * (f3 / 6.0);
    }
    return r;
}

/// exp of a jet, by the chain rule.
template <int Order>
Jet<Order> exp(const Jet<Order>& a) {
    const double e = std::exp(a.value());
    return chain(a, e, e, e, e);
}

/// log of a jet; the value of `a` must be positive.
template <int Order>
Jet<Order> log(const Jet<Order>& a) {
    const double a0 = a.value();
    return chain(a, std::log(a0), 1.0 / a0, -1.0 / (a0 * a0), 2.0 / (a0 * a0 * a0));
}

/**
 * @brief pow of a jet with a real exponent, by the chain rule; the value of `a` must be positive.
 *
 * At the value 0 with an exponent greater than 1 the zero jet is returned: through a composition
 * that is nilpotent at 0 (a squared factor such as (x - x0)^2, the form this code base feeds it)
 * that is the limit of every derivative of the outer function that exists, which covers the
 * x^e and ((delta-1)^2)^e factors of the Helmholtz formulations at their endpoints. Exponents in
 * (0, 1] at the value 0 have an infinite first derivative and produce NaN deliberately.
 */
template <int Order>
Jet<Order> pow(const Jet<Order>& a, double e) {
    const double a0 = a.value();
    if (a0 == 0.0 && e > 1.0) return Jet<Order>();
    const double p0 = std::pow(a0, e);
    return chain(a, p0, e * p0 / a0, e * (e - 1.0) * p0 / (a0 * a0), e * (e - 1.0) * (e - 2.0) * p0 / (a0 * a0 * a0));
}

/// pow of a jet with an integer exponent (small exponents by repeated multiplication).
template <int Order>
Jet<Order> pow(const Jet<Order>& a, int n) {
    if (n == 0) return Jet<Order>::constant(1.0);
    if (n == 1) return a;
    if (n == 2) return a * a;
    if (n == 3) return a * a * a;
    if (n < 0) return inverse(pow(a, -n));
    return pow(a, static_cast<double>(n));
}

/// sqrt of a jet; the value of `a` must be positive.
template <int Order>
Jet<Order> sqrt(const Jet<Order>& a) {
    return pow(a, 0.5);
}

}  // namespace coolsolve::taylor
