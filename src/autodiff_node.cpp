#include "coolsolve/autodiff_node.h"
#include <stdexcept>
#include <algorithm>

namespace coolsolve {

// ============================================================================
// Standard Math Function Evaluation
// ============================================================================

ADValue evaluateStandardFunction(const std::string& funcName, 
                                  const std::vector<ADValue>& args) {
    // Convert function name to lowercase for case-insensitive matching
    std::string name = funcName;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    // Constants (0 arguments)
    if (name == "pi") {
        size_t n = args.empty() ? 0 : args[0].gradient.size();
        // If we don't have n, we'll have to rely on the caller or a default
        // But in this context, we can return a constant with 0 gradient if n is unknown
        return ADValue::constant(3.14159265358979323846, n);
    }

    if (args.empty()) {
        throw std::runtime_error("Function " + funcName + " requires at least one argument");
    }
    
    // Variadic functions (min, max)
    if (name == "min") {
        ADValue result = args[0];
        for (size_t i = 1; i < args.size(); ++i) {
            result = min(result, args[i]);
        }
        return result;
    }
    if (name == "max") {
        ADValue result = args[0];
        for (size_t i = 1; i < args.size(); ++i) {
            result = max(result, args[i]);
        }
        return result;
    }
    
    // Aggregation functions (variadic)
    if (name == "sum" || name == "sum2d") {
        ADValue result = args[0];
        for (size_t i = 1; i < args.size(); ++i) {
            result = result + args[i];
        }
        return result;
    }
    if (name == "average") {
        ADValue result = args[0];
        for (size_t i = 1; i < args.size(); ++i) {
            result = result + args[i];
        }
        return result / ADValue::constant(static_cast<double>(args.size()), result.gradient.size());
    }
    if (name == "product") {
        ADValue result = args[0];
        for (size_t i = 1; i < args.size(); ++i) {
            result = result * args[i];
        }
        return result;
    }
    if (name == "stddev") {
        // Standard deviation: sqrt(sum((x_i - mean)^2) / N)
        double n = static_cast<double>(args.size());
        size_t nVars = args[0].gradient.size();
        // Compute mean
        ADValue mean = args[0];
        for (size_t i = 1; i < args.size(); ++i) mean = mean + args[i];
        mean = mean / ADValue::constant(n, nVars);
        // Compute sum of squared deviations
        ADValue sumSq = ADValue::constant(0.0, nVars);
        for (size_t i = 0; i < args.size(); ++i) {
            ADValue diff = args[i] - mean;
            sumSq = sumSq + diff * diff;
        }
        return sqrt(sumSq / ADValue::constant(n, nVars));
    }

    // Single argument functions
    if (args.size() == 1) {
        const ADValue& x = args[0];
        
        // EES convention: trig functions use degrees
        constexpr double deg2rad = 3.14159265358979323846 / 180.0;
        constexpr double rad2deg = 180.0 / 3.14159265358979323846;
        
        if (name == "sin") { ADValue xrad = x * ADValue::constant(deg2rad, x.gradient.size()); return sin(xrad); }
        if (name == "cos") { ADValue xrad = x * ADValue::constant(deg2rad, x.gradient.size()); return cos(xrad); }
        if (name == "tan") { ADValue xrad = x * ADValue::constant(deg2rad, x.gradient.size()); return tan(xrad); }
        if (name == "exp") return exp(x);
        if (name == "ln" || name == "log") return log(x);
        if (name == "log10") return log10(x);
        if (name == "sqrt") return sqrt(x);
        if (name == "abs") return abs(x);
        if (name == "asin" || name == "arcsin") return asin(x) * ADValue::constant(rad2deg, x.gradient.size());
        if (name == "acos" || name == "arccos") return acos(x) * ADValue::constant(rad2deg, x.gradient.size());
        if (name == "atan" || name == "arctan") return atan(x) * ADValue::constant(rad2deg, x.gradient.size());
        if (name == "sinh") return sinh(x);
        if (name == "cosh") return cosh(x);
        if (name == "tanh") return tanh(x);
        if (name == "arcsinh" || name == "asinh") return asinh(x);
        if (name == "arccosh" || name == "acosh") return acosh(x);
        if (name == "arctanh" || name == "atanh") return atanh(x);
        if (name == "erf") return erf(x);
        if (name == "erfc") return erfc(x);
        
        // Rounding functions (gradient = 0, not differentiable at integers)
        if (name == "ceil") return ADValue::constant(std::ceil(x.value), x.gradient.size());
        if (name == "floor") return ADValue::constant(std::floor(x.value), x.gradient.size());
        if (name == "round") return ADValue::constant(std::round(x.value), x.gradient.size());
        if (name == "trunc") return ADValue::constant(std::trunc(x.value), x.gradient.size());
        if (name == "sign") {
            double s = (x.value > 0.0) ? 1.0 : (x.value < 0.0) ? -1.0 : 0.0;
            return ADValue::constant(s, x.gradient.size());
        }
    }
    
    // Two argument functions
    if (args.size() == 2) {
        const ADValue& x = args[0];
        const ADValue& y = args[1];
        
        if (name == "pow") return pow(x, y);
        if (name == "mod") {
            // mod(x, y) = x - floor(x/y) * y; gradient passes through (not differentiable at multiples)
            double val = std::fmod(x.value, y.value);
            return ADValue::constant(val, x.gradient.size());
        }
        // LMTD_f(DT1, DT2) = (DT1 - DT2) / ln(DT1/DT2)
        // EES built-in log-mean temperature difference.  Handles the
        // singular case DT1 ≈ DT2 with a smooth limit.
        if (name == "lmtd_f" || name == "lmtd") {
            double d1 = x.value;
            double d2 = y.value;
            // Both DTs must have the same sign and be non-zero
            if (d1 * d2 <= 0.0) {
                // Not physical — return the arithmetic mean as a penalty
                return (x + y) / ADValue::constant(2.0, x.gradient.size());
            }
            // Smooth limit at d1 = d2: LMTD -> d1
            // Relative closeness check
            double ratio = d1 / d2;
            if (std::abs(ratio - 1.0) < 1e-6) {
                // Series: LMTD ≈ (d1+d2)/2 * [1 - ((d1-d2)^2)/(12*((d1+d2)/2)^2) + ...]
                return (x + y) / ADValue::constant(2.0, x.gradient.size());
            }
            // General case — do it entirely in AD
            return (x - y) / log(x / y);
        }
        if (name == "atan2") {
            // atan2(y, x) = atan(y/x) with proper quadrant handling
            // Returns result in degrees (EES convention)
            double val = std::atan2(x.value, y.value) * 180.0 / 3.14159265358979323846;
            double denom = x.value * x.value + y.value * y.value;
            double rad2deg = 180.0 / 3.14159265358979323846;
            ADValue z(val, x.gradient.size());
            for (size_t i = 0; i < z.gradient.size(); ++i) {
                z.gradient[i] = rad2deg * (y.value * x.gradient[i] - x.value * y.gradient[i]) / denom;
            }
            return z;
        }
    }
    
    // Three argument functions
    if (args.size() == 3) {
        // IF(condition, true_value, false_value)
        if (name == "if") {
            const ADValue& cond = args[0];
            return (cond.value > 0.0) ? args[1] : args[2];
        }
    }
    
    // If we get here, the function is not implemented
    throw std::runtime_error("Unknown or unsupported function: " + funcName + 
                            " with " + std::to_string(args.size()) + " arguments");
}

}  // namespace coolsolve
