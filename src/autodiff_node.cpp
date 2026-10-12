#include "coolsolve/autodiff_node.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <sstream>

namespace coolsolve {

namespace {

// ----------------------------------------------------------------------------
// Time functions of the EES library SETP.LIB (Duffie & Beckman, Solar Engineering of Thermal
// Processes): nDay_, NumDay_, AveDay_, EqnTime_. They count in a non-leap year (365 days).
// ----------------------------------------------------------------------------

/// Days of each month of a non-leap year (Table 1.6.1).
constexpr int monthLength[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
/// Day number n of the "average day" of each month (Table 1.6.1): the day whose extraterrestrial
/// radiation is about the monthly mean (not valid above the Arctic or below the Antarctic circle).
constexpr int averageDay[12] = {17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344};

std::string showValue(double v) {
    std::ostringstream os;
    os << v;
    return os.str();
}

/// `fname` takes exactly `expected` arguments (`parameters`: their names, for the message).
void checkDateArgCount(const char* fname, const char* parameters, size_t expected, size_t got) {
    if (got != expected)
        throw std::runtime_error(std::string(fname) + " takes " + std::to_string(expected) +
                                 (expected == 1 ? " argument (" : " arguments (") + parameters + "), got " +
                                 std::to_string(got));
}

/// 0-based index of `month` rounded to the nearest integer (the library's `monthx = round(month)`,
/// so 2.5 is month 3); an error unless the rounded month is in 1..12. Shared by the month functions.
int roundedMonthIndex(const char* fname, double month) {
    const double rounded = std::round(month);
    // the negated test also rejects NaN
    if (!(rounded >= 1.0 && rounded <= 12.0))
        throw std::runtime_error(std::string(fname) + ": month must be between 1 and 12 after rounding, got " +
                                 showValue(month));
    return static_cast<int>(rounded) - 1;
}

/// nDay_(month, day): day of the year, the days of the months before `month` plus `day`. `day` is
/// at least 1 and at most the length of the month plus 0.001 (so 28.0005 is accepted in February),
/// used as given; anything else is an error, not a value. The month is a piecewise-constant selector
/// with no derivative: the gradient is that of `day`.
ADValue nDayOfYear(const std::vector<ADValue>& args) {
    checkDateArgCount("nDay_", "month, day", 2, args.size());
    const int month = roundedMonthIndex("nDay_", args[0].value);
    const double day = args[1].value;
    // the negated tests also reject NaN
    if (!(day >= 1.0)) throw std::runtime_error("nDay_: day must be at least 1, got " + showValue(day));
    if (!(day <= monthLength[month] + 0.001))
        throw std::runtime_error("nDay_: day " + showValue(day) + " exceeds the " + std::to_string(monthLength[month]) +
                                 " days of month " + std::to_string(month + 1));
    int daysBefore = 0;
    for (int m = 0; m < month; ++m) daysBefore += monthLength[m];
    return args[1] + ADValue::constant(daysBefore, args[1].gradient.size());
}

/// NumDay_(month) and AveDay_(month): the entry of `table` for the month. Piecewise constant, so the
/// derivative is 0.
ADValue monthTableValue(const char* fname, const int (&table)[12], const std::vector<ADValue>& args) {
    checkDateArgCount(fname, "month", 1, args.size());
    return ADValue::constant(table[roundedMonthIndex(fname, args[0].value)], args[0].gradient.size());
}

/// EqnTime_(n): equation of time in minutes, Duffie & Beckman eq. 1.5.3, as a function of the day
/// number n, with B = (n-1)*360/365 degrees. Smooth, with no range check (as in the library).
ADValue equationOfTime(const std::vector<ADValue>& args) {
    checkDateArgCount("EqnTime_", "n", 1, args.size());
    constexpr double pi = 3.14159265358979323846;
    const ADValue b = (args[0] - 1.0) * (2.0 * pi / 365.0);   // radians
    const ADValue b2 = b * 2.0;
    return (0.000075 + 0.001868 * cos(b) - 0.032077 * sin(b) - 0.014615 * cos(b2) - 0.04089 * sin(b2)) * 229.2;
}

}  // namespace

// ============================================================================
// Colebrook-White friction factor (built-in procedure Colebrook of the EES library)
// ============================================================================

ADValue colebrookFrictionFactor(const ADValue& re, const ADValue& rr) {
    const double reynolds = re.value, roughness = rr.value;
    // the negated tests also reject NaN
    if (!(reynolds > 0.0) || !std::isfinite(reynolds))
        throw std::runtime_error("Colebrook: the Reynolds number Re must be positive (and finite), got " +
                                 showValue(reynolds));
    // RR >= 3.7 leaves 1/sqrt(f) = -2*log10(RR/3.7 + ...) no positive solution
    if (!(roughness >= 0.0) || !(roughness < 3.7))
        throw std::runtime_error("Colebrook: the relative roughness RR = e/D must be at least 0 and below 3.7, got " +
                                 showValue(roughness));

    // g(x) = x + 2*log10(a + b*x) = 0 with x = 1/sqrt(f): increasing and concave on x > 0, where its root
    // lies (g(0+) = 2*log10(a) < 0). Newton from the Swamee-Jain approximation; after the first step the
    // iterates approach the root from below, a step that leaves x > 0 is replaced by a halving.
    const double a = roughness / 3.7, b = 2.51 / reynolds;
    constexpr double twoOverLn10 = 2.0 / 2.30258509299404568402;
    double x = -2.0 * std::log10(a + 5.74 / std::pow(reynolds, 0.9));
    if (!(x > 0.0) || !std::isfinite(x)) x = 1.0;
    bool converged = false;
    for (int iteration = 0; iteration < 100 && !converged; ++iteration) {
        const double s = a + b * x;
        const double g = x + 2.0 * std::log10(s);
        const double dg = 1.0 + twoOverLn10 * b / s;
        double next = x - g / dg;
        if (!(next > 0.0)) next = 0.5 * x;
        converged = std::abs(next - x) <= 1e-14 * next;
        x = next;
    }
    if (!converged || !std::isfinite(x))
        throw std::runtime_error("Colebrook: the Newton iterations did not converge for Re = " + showValue(reynolds) +
                                 ", RR = " + showValue(roughness));

    // Derivatives by the implicit function theorem at the solution: dx = -(g_Re dRe + g_RR dRR)/g_x, f = x^-2
    const double s = a + b * x;
    const double gX = 1.0 + twoOverLn10 * b / s;
    const double gRe = -twoOverLn10 * b * x / (s * reynolds);   // d(b*x)/dRe = -b*x/Re
    const double gRR = twoOverLn10 / (3.7 * s);
    const double dfdx = -2.0 / (x * x * x);
    const double dfdRe = -dfdx * gRe / gX, dfdRR = -dfdx * gRR / gX;

    ADValue f(1.0 / (x * x), std::max(re.gradient.size(), rr.gradient.size()));
    for (size_t i = 0; i < re.gradient.size(); ++i) f.gradient[i] += dfdRe * re.gradient[i];
    for (size_t i = 0; i < rr.gradient.size(); ++i) f.gradient[i] += dfdRR * rr.gradient[i];
    return f;
}

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

    // Inline conditional (two EES forms, selected by the number of arguments):
    //   IF(cond, a, b)      a if cond > 0, else b
    //   IF(A, B, X, Y, Z)   X if A < B, Y if A = B (exact floating-point equality, as in EES), Z if A > B
    // Value and gradient are those of the selected branch; the tested arguments (cond, or A and B)
    // are piecewise-constant selectors and contribute no derivative.
    if (name == "if") {
        if (args.size() == 3) return (args[0].value > 0.0) ? args[1] : args[2];
        if (args.size() == 5) {
            const double a = args[0].value, b = args[1].value;
            return (a < b) ? args[2] : (a == b) ? args[3] : args[4];
        }
        throw std::runtime_error("IF takes 3 or 5 arguments: IF(cond, a, b) (a if cond > 0, else b) or "
                                 "IF(A, B, X, Y, Z) (X if A < B, Y if A = B, Z if A > B); got " +
                                 std::to_string(args.size()));
    }

    // EES time functions (library SETP.LIB), checked for their own number of arguments
    if (name == "nday_") return nDayOfYear(args);
    if (name == "numday_") return monthTableValue("NumDay_", monthLength, args);
    if (name == "aveday_") return monthTableValue("AveDay_", averageDay, args);
    if (name == "eqntime_") return equationOfTime(args);

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
    
    // If we get here, the function is not implemented
    throw std::runtime_error("Unknown or unsupported function: " + funcName + 
                            " with " + std::to_string(args.size()) + " arguments");
}

}  // namespace coolsolve
