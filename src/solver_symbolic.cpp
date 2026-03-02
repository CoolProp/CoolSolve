/**
 * @file solver_symbolic.cpp
 * @brief Symbolic Block Reduction — implementation.
 *
 * Reduces the size of algebraic blocks before iterative solving by:
 *  1. Extracting explicit equations (all RHS vars are external).
 *  2. Inverting CoolProp calls so that an unknown input becomes the output.
 *  3. Substituting single-use variables into other equations.
 *
 * See symbolic_reduction.h for the full description.
 */

#include "coolsolve/symbolic_reduction.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace coolsolve {

// ============================================================================
// Internal helpers
// ============================================================================

/// Case-insensitive string comparison.
static bool ciEqual(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    }
    return true;
}

/// Case-insensitive set membership test.
static bool ciContains(const std::set<std::string, CaseInsensitiveLess>& s,
                       const std::string& key) {
    return s.find(key) != s.end();
}

static bool ciContainsVec(const std::vector<std::string>& v,
                          const std::string& key) {
    for (const auto& s : v) {
        if (ciEqual(s, key)) return true;
    }
    return false;
}

/// Lower-case a string.
static std::string toLower(const std::string& s) {
    std::string r = s;
    std::transform(r.begin(), r.end(), r.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return r;
}

// ============================================================================
// CoolProp inversion maps
// ============================================================================

/// Map lower-case CoolProp function names to the canonical input parameter name
/// that represents the same property.  For example, enthalpy → "H".
static const std::unordered_map<std::string, std::string>& funcToInputMap() {
    static const std::unordered_map<std::string, std::string> m = {
        {"enthalpy",    "H"},
        {"h",           "H"},
        {"entropy",     "S"},
        {"s",           "S"},
        {"temperature", "T"},
        {"t",           "T"},
        {"pressure",    "P"},
        {"p",           "P"},
        {"density",     "D"},
        {"rho",         "D"},
        {"d",           "D"},
        {"quality",     "Q"},
        {"x",           "Q"},
        {"q",           "Q"},
        {"internalenergy", "U"},
        {"u",           "U"},
        {"intenergy",   "U"},
    };
    return m;
}

/// Map lower-case input parameter names to the corresponding CoolProp output
/// function name.  For example, "h" → "enthalpy".
static const std::unordered_map<std::string, std::string>& inputToFuncMap() {
    static const std::unordered_map<std::string, std::string> m = {
        {"h",           "enthalpy"},
        {"enthalpy",    "enthalpy"},
        {"s",           "entropy"},
        {"entropy",     "entropy"},
        {"t",           "temperature"},
        {"temperature", "temperature"},
        {"p",           "pressure"},
        {"pressure",    "pressure"},
        {"d",           "density"},
        {"rho",         "density"},
        {"density",     "density"},
        {"q",           "quality"},
        {"x",           "quality"},
        {"quality",     "quality"},
        {"u",           "internalenergy"},
        {"internalenergy", "internalenergy"},
        {"intenergy",   "internalenergy"},
    };
    return m;
}

/// Valid CoolProp input pairs (as sets of two canonical parameter letters).
/// These are the pairs that CoolProp's AbstractState::update() supports.
/// Each pair is stored as {lower(p1), lower(p2)} with p1 < p2 lexicographically.
static bool isValidPairCanonical(const std::string& a, const std::string& b) {
    // Canonical single-letter upper-case names: D H P Q S T U
    std::string lo = (a < b) ? a : b;
    std::string hi = (a < b) ? b : a;

    // Supported pairs (from CoolProp documentation and generate_update_pair):
    //   PT, HP, PS, HS, DP, DT, QT, QP, UP, UT, DS, HS, HU…
    // We use a whitelist of commonly supported pairs.
    static const std::set<std::pair<std::string,std::string>> valid = {
        {"P", "T"},   // PT_INPUTS
        {"H", "P"},   // HmassP_INPUTS
        {"P", "S"},   // PSmass_INPUTS
        {"H", "S"},   // HmassSmass_INPUTS
        {"D", "P"},   // DmassP_INPUTS
        {"D", "T"},   // DmassT_INPUTS
        {"Q", "T"},   // QT_INPUTS
        {"P", "Q"},   // PQ_INPUTS
        {"D", "H"},   // DmassHmass_INPUTS
        {"D", "S"},   // DmassSmass_INPUTS
        {"H", "T"},   // HmassT — NOTE: NOT supported by CoolProp for most backends
        {"P", "U"},   // PUmass_INPUTS
        {"D", "U"},   // DmassUmass_INPUTS
    };

    // Exclude H,T — CoolProp generally does not support this pair
    if ((lo == "H" && hi == "T") || (lo == "T" && hi == "H"))
        return false;

    return valid.count({lo, hi}) > 0 || valid.count({hi, lo}) > 0;
}

// ============================================================================
// Public helpers
// ============================================================================

std::string coolpropFuncToInputName(const std::string& funcName) {
    auto it = funcToInputMap().find(toLower(funcName));
    if (it != funcToInputMap().end()) return it->second;
    return "";
}

std::string coolpropInputToFuncName(const std::string& inputName) {
    auto it = inputToFuncMap().find(toLower(inputName));
    if (it != inputToFuncMap().end()) return it->second;
    return "";
}

bool isValidCoolPropInputPair(const std::string& param1Name,
                              const std::string& param2Name) {
    // Canonicalise to single upper-case letters
    std::string c1 = coolpropFuncToInputName(coolpropInputToFuncName(param1Name));
    if (c1.empty()) {
        // Try direct mapping for input names like "T", "P"
        auto it = funcToInputMap().find(toLower(param1Name));
        if (it != funcToInputMap().end()) c1 = it->second;
    }
    // For simple single-letter names
    if (c1.empty() && param1Name.size() <= 2) {
        c1 = param1Name;
        std::transform(c1.begin(), c1.end(), c1.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    }

    std::string c2 = coolpropFuncToInputName(coolpropInputToFuncName(param2Name));
    if (c2.empty()) {
        auto it = funcToInputMap().find(toLower(param2Name));
        if (it != funcToInputMap().end()) c2 = it->second;
    }
    if (c2.empty() && param2Name.size() <= 2) {
        c2 = param2Name;
        std::transform(c2.begin(), c2.end(), c2.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    }

    if (c1.empty() || c2.empty()) return false;
    return isValidPairCanonical(c1, c2);
}

// ============================================================================
// Collect variables used in an expression
// ============================================================================

static void collectExprVars(const ExprPtr& expr,
                            std::set<std::string, CaseInsensitiveLess>& out) {
    if (!expr) return;
    if (expr->is<Variable>()) {
        out.insert(expr->as<Variable>().flattenedName());
    } else if (expr->is<UnaryOp>()) {
        collectExprVars(expr->as<UnaryOp>().operand, out);
    } else if (expr->is<BinaryOp>()) {
        collectExprVars(expr->as<BinaryOp>().left, out);
        collectExprVars(expr->as<BinaryOp>().right, out);
    } else if (expr->is<FunctionCall>()) {
        const auto& fc = expr->as<FunctionCall>();
        for (const auto& a : fc.args) collectExprVars(a, out);
        for (const auto& [n, e] : fc.namedArgs) collectExprVars(e, out);
    }
    // NumberLiteral, StringLiteral: no variables
}

// ============================================================================
// Detect CoolProp function call in an expression
// ============================================================================

/**
 * @brief Check if an expression is a CoolProp thermodynamic function call.
 *
 * Returns true if the expression is a FunctionCall with:
 *  - A recognized CoolProp function name (enthalpy, temperature, etc.)
 *  - Exactly one positional argument (the fluid name)
 *  - Exactly two named arguments (the input properties)
 */
static bool isCoolPropCall(const ExprPtr& expr) {
    if (!expr || !expr->is<FunctionCall>()) return false;
    const auto& fc = expr->as<FunctionCall>();
    if (fc.namedArgs.size() != 2) return false;
    if (fc.args.empty()) return false;
    // Check if function name is a recognized thermo function
    return !coolpropFuncToInputName(fc.name).empty();
}

// ============================================================================
// Analyse one block
// ============================================================================

BlockReductionResult analyseBlockReduction(
    const Block& block,
    const IR& ir,
    const StructuralAnalysisResult& analysis) {

    BlockReductionResult result;

    // Nothing to reduce for size <= 1
    if (block.size() <= 1) return result;

    const auto& equations = ir.getEquations();

    // --- Build working copies of the block's equation/variable sets ---
    std::set<std::string, CaseInsensitiveLess> blockVars(
        block.variables.begin(), block.variables.end());
    std::set<int> blockEqs(block.equationIds.begin(), block.equationIds.end());

    // Set of variables that are "known" — external to the block or already
    // extracted by a previous reduction step.
    // We start with ALL model variables minus the block variables.
    std::set<std::string, CaseInsensitiveLess> knownVars;
    for (const auto& [name, info] : ir.getVariables()) {
        if (!ciContains(blockVars, name)) {
            knownVars.insert(name);
        }
    }

    // Map from equation ID → its matched output variable
    std::map<int, std::string> eqToOutput;
    for (int eqId : block.equationIds) {
        if (eqId >= 0 && eqId < static_cast<int>(analysis.matching.size())) {
            eqToOutput[eqId] = analysis.matching[eqId];
        }
    }

    // Track which equations and variables have been extracted
    std::set<int> extractedEqs;
    std::set<std::string, CaseInsensitiveLess> extractedVars;

    // Collect pre-solve and post-solve steps
    std::vector<ReductionStep> preSolve;
    std::vector<ReductionStep> postSolve;

    // =======================================================================
    // Iterative reduction loop — repeat until no more reductions possible
    // =======================================================================
    bool changed = true;
    int maxPasses = static_cast<int>(block.size());  // Safety bound
    int pass = 0;

    while (changed && pass < maxPasses) {
        changed = false;
        ++pass;

        // -------------------------------------------------------------------
        // Phase 1: CoolProp call inversion
        // -------------------------------------------------------------------
        // For each non-extracted equation in the block, check if it's of the
        // form:  output_var = thermo_func(fluid, A=a, B=b)
        // where output_var is the matched output, one of {a, b} is a block
        // variable, and the other is known.  If the output property can serve
        // as an input and the unknown input can become the new output, invert.
        for (int eqId : block.equationIds) {
            if (extractedEqs.count(eqId)) continue;
            if (eqId < 0 || eqId >= static_cast<int>(equations.size())) continue;

            const auto& eq = equations[eqId];

            // The equation must be: LHS = RHS  with LHS being a simple variable
            if (!eq.lhs || !eq.lhs->is<Variable>()) continue;
            const std::string lhsVar = eq.lhs->as<Variable>().flattenedName();

            // LHS must be the matched output variable
            auto outIt = eqToOutput.find(eqId);
            if (outIt == eqToOutput.end()) continue;
            if (!ciEqual(outIt->second, lhsVar)) continue;

            // RHS must be a CoolProp call
            if (!isCoolPropCall(eq.rhs)) continue;

            const auto& fc = eq.rhs->as<FunctionCall>();
            const std::string outputPropName = coolpropFuncToInputName(fc.name);
            if (outputPropName.empty()) continue;

            // LHS variable must still be in the block (not yet extracted)
            bool lhsIsBlockVar = ciContains(blockVars, lhsVar) &&
                                 !ciContains(extractedVars, lhsVar);

            // Analyse the two named arguments
            for (int argIdx = 0; argIdx < 2; ++argIdx) {
                const auto& [inputName, inputExpr] = fc.namedArgs[argIdx];
                const int otherIdx = 1 - argIdx;
                const auto& [otherName, otherExpr] = fc.namedArgs[otherIdx];

                // Collect variables in each named arg
                std::set<std::string, CaseInsensitiveLess> argVars, otherArgVars;
                collectExprVars(inputExpr, argVars);
                collectExprVars(otherExpr, otherArgVars);

                // The "target" input must contain exactly one block variable
                // and no other block variables.
                std::string targetBlockVar;
                bool targetIsSimpleBlockVar = false;

                if (argVars.size() == 1) {
                    const std::string& v = *argVars.begin();
                    if (ciContains(blockVars, v) && !ciContains(extractedVars, v)) {
                        // Check the expression is just a simple variable (not a complex expr)
                        if (inputExpr->is<Variable>()) {
                            targetBlockVar = v;
                            targetIsSimpleBlockVar = true;
                        }
                    }
                }
                if (!targetIsSimpleBlockVar) continue;

                // The "other" input must be fully known (all its variables are
                // known or already extracted from the block)
                bool otherKnown = true;
                for (const auto& v : otherArgVars) {
                    if (ciContains(blockVars, v) && !ciContains(extractedVars, v) &&
                        !ciContains(knownVars, v)) {
                        otherKnown = false;
                        break;
                    }
                }
                if (!otherKnown) continue;

                // The LHS (output) must be known OR be the other block unknown
                // that we're keeping in the block.  For inversion, we need the
                // output value to become an INPUT to the inverted call, so it
                // must be known.
                bool lhsKnown = ciContains(knownVars, lhsVar) ||
                                ciContains(extractedVars, lhsVar);
                // If LHS is not known, but it's a block variable that could be
                // determined by another equation, we can still invert IF we make
                // this a post-solve step (after the reduced block determines lhsVar).
                bool isPostSolve = false;
                if (!lhsKnown && lhsIsBlockVar) {
                    // LHS is a block var that will be solved by the reduced block.
                    // We can post-solve the target variable after.
                    isPostSolve = true;
                    lhsKnown = true;  // treat as known for inversion feasibility
                }
                if (!lhsKnown) continue;

                // Check if the inverted call forms a valid CoolProp input pair:
                // new inputs = (outputPropName, otherInputName)
                // new output func = the function for targetInput
                std::string canonicalOther = otherName;  // e.g. "P", "T"
                if (!isValidCoolPropInputPair(outputPropName, canonicalOther))
                    continue;

                // Build the inverted function name
                std::string invertedFunc = coolpropInputToFuncName(inputName);
                if (invertedFunc.empty()) continue;

                // --- Inversion is feasible ---
                ReductionStep step;
                step.variable = targetBlockVar;
                step.equationId = eqId;
                step.inverted = true;
                step.invertedFuncName = invertedFunc;
                step.fluidExpr = fc.args.empty() ? nullptr : fc.args[0];

                if (step.fluidExpr && step.fluidExpr->is<StringLiteral>()) {
                    step.fluidName = step.fluidExpr->as<StringLiteral>().value;
                }

                // New named args: (outputPropName, lhsVarExpr), (otherName, otherExpr)
                // Build a Variable expression for the LHS var
                auto lhsExpr = std::make_shared<Expression>();
                lhsExpr->node = Variable{lhsVar, {}};

                step.newNamedArgs.push_back({outputPropName, lhsExpr});
                step.newNamedArgs.push_back({otherName, otherExpr});

                if (isPostSolve) {
                    postSolve.push_back(std::move(step));
                } else {
                    preSolve.push_back(std::move(step));
                }

                // Mark equation and target variable as extracted
                extractedEqs.insert(eqId);
                extractedVars.insert(targetBlockVar);
                knownVars.insert(targetBlockVar);

                // If LHS was also a block variable and it's now the only
                // user in this equation, it stays in the block for other eqs.
                // But if this equation was its ONLY equation, we might need to
                // handle that.  For now, just mark the equation extracted.

                result.inversionsApplied++;
                changed = true;
                break;  // Move to next equation
            }
        }

        // -------------------------------------------------------------------
        // Phase 2: Explicit extraction
        // -------------------------------------------------------------------
        // Find equations where the matched output variable's RHS only references
        // known variables (external or already-extracted).
        for (int eqId : block.equationIds) {
            if (extractedEqs.count(eqId)) continue;
            if (eqId < 0 || eqId >= static_cast<int>(equations.size())) continue;

            const auto& eq = equations[eqId];

            // Must be: LHS_var = RHS where LHS is the matched output
            if (!eq.lhs || !eq.lhs->is<Variable>()) continue;
            const std::string lhsVar = eq.lhs->as<Variable>().flattenedName();
            auto outIt = eqToOutput.find(eqId);
            if (outIt == eqToOutput.end()) continue;
            if (!ciEqual(outIt->second, lhsVar)) continue;

            // LHS must still be a block variable
            if (!ciContains(blockVars, lhsVar) || ciContains(extractedVars, lhsVar))
                continue;

            // Check if ALL variables in the RHS are known
            std::set<std::string, CaseInsensitiveLess> rhsVars;
            collectExprVars(eq.rhs, rhsVars);

            bool allKnown = true;
            for (const auto& v : rhsVars) {
                if (ciContains(blockVars, v) && !ciContains(extractedVars, v) &&
                    !ciContains(knownVars, v)) {
                    allKnown = false;
                    break;
                }
            }
            if (!allKnown) continue;

            // --- Extract this equation ---
            ReductionStep step;
            step.variable = lhsVar;
            step.equationId = eqId;
            step.inverted = false;

            preSolve.push_back(std::move(step));

            extractedEqs.insert(eqId);
            extractedVars.insert(lhsVar);
            knownVars.insert(lhsVar);

            result.extractionsApplied++;
            changed = true;
        }

        // -------------------------------------------------------------------
        // Phase 3: Equation substitution
        // -------------------------------------------------------------------
        // If variable x is the output of equation eq_A (x = f(...)) and x
        // appears in exactly one other block equation eq_B, we can conceptually
        // substitute f(...) for x in eq_B and remove eq_A + x.
        //
        // We don't literally rewrite the AST here; instead, if after
        // substitution eq_A becomes "extractable" (all its RHS vars are known),
        // we extract it (which effectively turns x into a pre/post-solve step).
        //
        // For now we use a simpler heuristic: count how many remaining block
        // equations reference each remaining block variable.  If a variable
        // appears in only its own output equation, it can be post-solved.
        for (int eqId : block.equationIds) {
            if (extractedEqs.count(eqId)) continue;
            if (eqId < 0 || eqId >= static_cast<int>(equations.size())) continue;

            const auto& eq = equations[eqId];
            auto outIt = eqToOutput.find(eqId);
            if (outIt == eqToOutput.end()) continue;

            const std::string& outVar = outIt->second;
            if (ciContains(extractedVars, outVar)) continue;

            // Must be: LHS = RHS with LHS being the output variable
            if (!eq.lhs || !eq.lhs->is<Variable>()) continue;
            if (!ciEqual(eq.lhs->as<Variable>().flattenedName(), outVar)) continue;

            // Count how many OTHER remaining block equations reference outVar
            int refCount = 0;
            for (int otherEqId : block.equationIds) {
                if (otherEqId == eqId || extractedEqs.count(otherEqId)) continue;
                if (otherEqId < 0 || otherEqId >= static_cast<int>(equations.size())) continue;
                const auto& otherEq = equations[otherEqId];
                if (ciContains(otherEq.variables, outVar)) {
                    ++refCount;
                }
            }

            // If outVar appears in zero other equations, it only appears in its
            // own defining equation.  Check if ALL RHS variables of this equation
            // will be known after the reduced block is solved.
            if (refCount == 0) {
                // This variable doesn't participate in any coupling — it can be
                // post-solved after the block (its RHS vars will be known then).
                std::set<std::string, CaseInsensitiveLess> rhsVars;
                collectExprVars(eq.rhs, rhsVars);

                // All RHS vars must be either known now, extracted, or remaining
                // in the reduced block (which will be solved).  Since refCount==0,
                // this equation doesn't help other equations, so we extract it.
                bool allResolvable = true;
                bool allKnownNow = true;
                for (const auto& v : rhsVars) {
                    if (ciContains(blockVars, v) && !ciContains(extractedVars, v) &&
                        !ciContains(knownVars, v)) {
                        allKnownNow = false;
                        // But will it be solved by the reduced block? Yes, if it's
                        // a remaining block variable.
                        // It's fine as post-solve.
                    }
                }

                ReductionStep step;
                step.variable = outVar;
                step.equationId = eqId;
                step.inverted = false;

                if (allKnownNow) {
                    preSolve.push_back(std::move(step));
                } else {
                    postSolve.push_back(std::move(step));
                }

                extractedEqs.insert(eqId);
                extractedVars.insert(outVar);

                if (allKnownNow) {
                    knownVars.insert(outVar);
                }

                result.substitutionsApplied++;
                changed = true;
            }
        }
    }

    // =======================================================================
    // Build result
    // =======================================================================
    if (extractedEqs.empty()) {
        // No reductions found
        result.reduced = false;
        result.reducedEquationIds = block.equationIds;
        result.reducedVariables = block.variables;
        return result;
    }

    result.reduced = true;
    result.preSolveSteps = std::move(preSolve);
    result.postSolveSteps = std::move(postSolve);

    // Build remaining equations/variables
    for (int eqId : block.equationIds) {
        if (!extractedEqs.count(eqId)) {
            result.reducedEquationIds.push_back(eqId);
        }
    }
    for (const auto& v : block.variables) {
        if (!ciContains(extractedVars, v)) {
            result.reducedVariables.push_back(v);
        }
    }

    return result;
}

}  // namespace coolsolve
