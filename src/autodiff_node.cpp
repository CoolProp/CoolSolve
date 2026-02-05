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

    // Single argument functions
    if (args.size() == 1) {
        const ADValue& x = args[0];
        
        if (name == "sin") return sin(x);
        if (name == "cos") return cos(x);
        if (name == "tan") return tan(x);
        if (name == "exp") return exp(x);
        if (name == "ln" || name == "log") return log(x);
        if (name == "log10") return log10(x);
        if (name == "sqrt") return sqrt(x);
        if (name == "abs") return abs(x);
        if (name == "asin" || name == "arcsin") return asin(x);
        if (name == "acos" || name == "arccos") return acos(x);
        if (name == "atan" || name == "arctan") return atan(x);
        if (name == "sinh") return sinh(x);
        if (name == "cosh") return cosh(x);
        if (name == "tanh") return tanh(x);
    }
    
    // Two argument functions
    if (args.size() == 2) {
        const ADValue& x = args[0];
        const ADValue& y = args[1];
        
        if (name == "pow") return pow(x, y);
        if (name == "atan2") {
            // atan2(y, x) = atan(y/x) with proper quadrant handling
            // For simplicity, we compute it numerically and use chain rule
            // d(atan2(y,x))/dy = x / (x^2 + y^2)
            // d(atan2(y,x))/dx = -y / (x^2 + y^2)
            double val = std::atan2(x.value, y.value);  // Note: atan2(y, x) in standard math
            double denom = x.value * x.value + y.value * y.value;
            ADValue z(val, x.gradient.size());
            for (size_t i = 0; i < z.gradient.size(); ++i) {
                z.gradient[i] = (y.value * x.gradient[i] - x.value * y.gradient[i]) / denom;
            }
            return z;
        }
    }
    
    // If we get here, the function is not implemented
    throw std::runtime_error("Unknown or unsupported function: " + funcName + 
                            " with " + std::to_string(args.size()) + " arguments");
}

}  // namespace coolsolve
