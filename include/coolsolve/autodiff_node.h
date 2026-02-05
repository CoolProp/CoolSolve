#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <stdexcept>

namespace coolsolve {

// ============================================================================
// Forward-Mode Automatic Differentiation Value
// ============================================================================

/**
 * @brief A dual number for forward-mode automatic differentiation.
 * 
 * Each ADValue stores:
 * - value: The computed value
 * - gradient: The partial derivatives w.r.t. all independent variables
 * 
 * The gradient vector size equals the number of variables in the current block.
 */
struct ADValue {
    double value = 0.0;
    std::vector<double> gradient;
    
    ADValue() = default;
    explicit ADValue(double v) : value(v) {}
    ADValue(double v, size_t numVars) : value(v), gradient(numVars, 0.0) {}
    ADValue(double v, std::vector<double> grad) : value(v), gradient(std::move(grad)) {}
    
    // Initialize as an independent variable (x_j has gradient = unit vector e_j)
    static ADValue independent(double value, size_t varIndex, size_t numVars) {
        ADValue ad(value, numVars);
        ad.gradient[varIndex] = 1.0;
        return ad;
    }
    
    // Initialize as a constant (gradient = zero vector)
    static ADValue constant(double value, size_t numVars) {
        return ADValue(value, numVars);
    }
};

// ============================================================================
// ADValue Arithmetic Operations (Forward-Mode AD Rules)
// ============================================================================

// Addition: z = x + y
// z.value = x.value + y.value
// z.gradient = x.gradient + y.gradient
inline ADValue operator+(const ADValue& x, const ADValue& y) {
    ADValue z(x.value + y.value, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = x.gradient[i] + y.gradient[i];
    }
    return z;
}

// Subtraction: z = x - y
inline ADValue operator-(const ADValue& x, const ADValue& y) {
    ADValue z(x.value - y.value, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = x.gradient[i] - y.gradient[i];
    }
    return z;
}

// Multiplication: z = x * y
// z.gradient = x.gradient * y.value + y.gradient * x.value
inline ADValue operator*(const ADValue& x, const ADValue& y) {
    ADValue z(x.value * y.value, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = x.gradient[i] * y.value + y.gradient[i] * x.value;
    }
    return z;
}

// Division: z = x / y
// z.gradient = (x.gradient * y.value - y.gradient * x.value) / (y.value^2)
inline ADValue operator/(const ADValue& x, const ADValue& y) {
    double yv2 = y.value * y.value;
    ADValue z(x.value / y.value, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = (x.gradient[i] * y.value - y.gradient[i] * x.value) / yv2;
    }
    return z;
}

// Unary negation: z = -x
inline ADValue operator-(const ADValue& x) {
    ADValue z(-x.value, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = -x.gradient[i];
    }
    return z;
}

// Power: z = x^y (where y is a constant double)
// z.gradient = y * x^(y-1) * x.gradient
inline ADValue pow(const ADValue& x, double y) {
    double val = std::pow(x.value, y);
    double dval = y * std::pow(x.value, y - 1);
    ADValue z(val, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// Power: z = x^y (where both are ADValues)
// z = exp(y * ln(x))
// z.gradient = x^y * (y.gradient * ln(x) + y * x.gradient / x)
inline ADValue pow(const ADValue& x, const ADValue& y) {
    double val = std::pow(x.value, y.value);
    double lnx = std::log(x.value);
    ADValue z(val, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = val * (y.gradient[i] * lnx + y.value * x.gradient[i] / x.value);
    }
    return z;
}

// ============================================================================
// Transcendental Functions
// ============================================================================

// sin(x)
inline ADValue sin(const ADValue& x) {
    ADValue z(std::sin(x.value), x.gradient.size());
    double dval = std::cos(x.value);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// cos(x)
inline ADValue cos(const ADValue& x) {
    ADValue z(std::cos(x.value), x.gradient.size());
    double dval = -std::sin(x.value);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// tan(x)
inline ADValue tan(const ADValue& x) {
    double c = std::cos(x.value);
    ADValue z(std::tan(x.value), x.gradient.size());
    double dval = 1.0 / (c * c);  // sec^2(x)
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// exp(x)
inline ADValue exp(const ADValue& x) {
    double val = std::exp(x.value);
    ADValue z(val, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = val * x.gradient[i];
    }
    return z;
}

// log(x) - natural logarithm
inline ADValue log(const ADValue& x) {
    ADValue z(std::log(x.value), x.gradient.size());
    double dval = 1.0 / x.value;
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// log10(x)
inline ADValue log10(const ADValue& x) {
    static const double ln10 = std::log(10.0);
    ADValue z(std::log10(x.value), x.gradient.size());
    double dval = 1.0 / (x.value * ln10);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// sqrt(x)
inline ADValue sqrt(const ADValue& x) {
    double val = std::sqrt(x.value);
    ADValue z(val, x.gradient.size());
    double dval = 0.5 / val;
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// abs(x)
inline ADValue abs(const ADValue& x) {
    ADValue z(std::abs(x.value), x.gradient.size());
    double sign = (x.value >= 0) ? 1.0 : -1.0;
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = sign * x.gradient[i];
    }
    return z;
}

// asin(x)
inline ADValue asin(const ADValue& x) {
    ADValue z(std::asin(x.value), x.gradient.size());
    double dval = 1.0 / std::sqrt(1.0 - x.value * x.value);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// acos(x)
inline ADValue acos(const ADValue& x) {
    ADValue z(std::acos(x.value), x.gradient.size());
    double dval = -1.0 / std::sqrt(1.0 - x.value * x.value);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// atan(x)
inline ADValue atan(const ADValue& x) {
    ADValue z(std::atan(x.value), x.gradient.size());
    double dval = 1.0 / (1.0 + x.value * x.value);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// sinh(x)
inline ADValue sinh(const ADValue& x) {
    ADValue z(std::sinh(x.value), x.gradient.size());
    double dval = std::cosh(x.value);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// cosh(x)
inline ADValue cosh(const ADValue& x) {
    ADValue z(std::cosh(x.value), x.gradient.size());
    double dval = std::sinh(x.value);
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// tanh(x)
inline ADValue tanh(const ADValue& x) {
    double th = std::tanh(x.value);
    ADValue z(th, x.gradient.size());
    double dval = 1.0 - th * th;  // sech^2(x)
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = dval * x.gradient[i];
    }
    return z;
}

// min(x, y) - note: not differentiable at x == y
inline ADValue min(const ADValue& x, const ADValue& y) {
    if (x.value <= y.value) {
        return x;
    }
    return y;
}

// max(x, y) - note: not differentiable at x == y
inline ADValue max(const ADValue& x, const ADValue& y) {
    if (x.value >= y.value) {
        return x;
    }
    return y;
}

// ============================================================================
// Mixed operations with scalars
// ============================================================================

inline ADValue operator+(const ADValue& x, double c) {
    ADValue z(x.value + c, x.gradient);
    return z;
}

inline ADValue operator+(double c, const ADValue& x) {
    return x + c;
}

inline ADValue operator-(const ADValue& x, double c) {
    ADValue z(x.value - c, x.gradient);
    return z;
}

inline ADValue operator-(double c, const ADValue& x) {
    ADValue z(c - x.value, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = -x.gradient[i];
    }
    return z;
}

inline ADValue operator*(const ADValue& x, double c) {
    ADValue z(x.value * c, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = x.gradient[i] * c;
    }
    return z;
}

inline ADValue operator*(double c, const ADValue& x) {
    return x * c;
}

inline ADValue operator/(const ADValue& x, double c) {
    ADValue z(x.value / c, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = x.gradient[i] / c;
    }
    return z;
}

inline ADValue operator/(double c, const ADValue& x) {
    double xv2 = x.value * x.value;
    ADValue z(c / x.value, x.gradient.size());
    for (size_t i = 0; i < z.gradient.size(); ++i) {
        z.gradient[i] = -c * x.gradient[i] / xv2;
    }
    return z;
}

// ============================================================================
// Standard Math Function Type
// ============================================================================

// A type alias for math functions that work with ADValues
using ADFunction = std::function<ADValue(const std::vector<ADValue>&)>;

// Registry of standard math functions (sin, cos, etc.)
ADValue evaluateStandardFunction(const std::string& funcName, 
                                  const std::vector<ADValue>& args);

}  // namespace coolsolve
