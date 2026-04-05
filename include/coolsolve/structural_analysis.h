#pragma once

#include "ir.h"
#include "diagnostic.h"
#include <vector>
#include <optional>

namespace coolsolve {

// ============================================================================
// Block Information (Strongly Connected Component)
// ============================================================================

struct Block {
    int id;                                // Block number (0-based)
    std::vector<int> equationIds;          // Equations in this block
    std::vector<std::string> variables;    // Variables determined by this block
    
    // If block size is 1, equation can be solved explicitly
    bool isExplicit() const { return equationIds.size() == 1; }
    
    // Size of the algebraic loop (0 for explicit equations)
    size_t size() const { return equationIds.size(); }
};

// ============================================================================
// Structural Analysis Result
// ============================================================================

struct StructuralAnalysisResult {
    bool success;
    std::string errorMessage;
    
    // Matching: equationId -> variableId (the output variable for that equation)
    std::vector<std::string> matching;  // matching[eqIdx] = variable name
    
    // Blocks in topological order (solve Block 0 first, then Block 1, etc.)
    std::vector<Block> blocks;
    
    // Statistics
    int totalEquations = 0;
    int totalVariables = 0;
    int totalBlocks = 0;
    int largestBlockSize = 0;
    int explicitEquationCount = 0;  // Number of Block-0 style equations
    
    // Comparison with EES reference
    struct ComparisonResult {
        bool matches = false;
        int eesBlockCount = 0;
        int ourBlockCount = 0;
        std::vector<std::string> differences;
    };
    std::optional<ComparisonResult> eesComparison;

    // Diagnostics (unmatched variables, block info, etc.)
    DiagnosticCollector diagnostics;
};

// ============================================================================
// Structural Analyzer
// ============================================================================

class StructuralAnalyzer {
public:
    // Perform structural analysis on the IR
    static StructuralAnalysisResult analyze(IR& ir);
    
    // Compare our result with EES residuals file
    static StructuralAnalysisResult::ComparisonResult compareWithEES(
        const StructuralAnalysisResult& result,
        const std::string& residualsFilePath);

    /// Re-decompose a reduced block into sub-SCCs after symbolic reduction.
    /// Builds a local dependency graph on the reduced equations/variables,
    /// runs Tarjan's SCC algorithm, and returns the sub-blocks in
    /// topological solve order.  If the reduced block is still a single
    /// SCC, returns a vector of size 1.
    static std::vector<Block> redecomposeBlock(
        const std::vector<int>& reducedEquationIds,
        const std::vector<std::string>& reducedVariables,
        const IR& ir,
        const StructuralAnalysisResult& analysis);
    
private:
    // Hopcroft-Karp algorithm for maximum bipartite matching
    static std::vector<int> hopcroftKarp(
        int numEquations,
        int numVariables,
        const std::vector<std::set<std::string, CaseInsensitiveLess>>& incidenceMatrix,
        const std::map<std::string, int, CaseInsensitiveLess>& varNameToIndex);
    
    // Tarjan's algorithm for finding Strongly Connected Components
    static std::vector<std::vector<int>> tarjanSCC(
        int numNodes,
        const std::vector<std::vector<int>>& adjacencyList);
    
    // Build dependency graph from matching
    static std::vector<std::vector<int>> buildDependencyGraph(
        const IR& ir,
        const std::vector<int>& matching);
};

// ============================================================================
// Tearing (feedback vertex set for algebraic loops)
// ============================================================================

/**
 * Result of tear set computation for a single block.
 * Used when enableTearing is true: tear variables are iterated with Newton;
 * the acyclic part is solved in topoOrderNonTearEqIds order.
 */
struct BlockTearSetResult {
    std::vector<std::string> tearVarNames;       // Variables to tear (output of tear equations)
    std::vector<int> topoOrderNonTearEqIds;      // Global equation IDs in dependency order (acyclic)
    std::vector<int> tearEquationIds;           // Global equation IDs that define tear variables
};

/**
 * Compute a tear set for a block using a greedy feedback-vertex-set heuristic.
 * Breaks the block's equation dependency graph so the remaining equations form a DAG.
 * @param block Block (SCC) to tear
 * @param ir IR for incidence (which variables appear in which equation)
 * @return Tear set result, or empty tearVarNames if block size < 2
 */
BlockTearSetResult computeBlockTearSet(const Block& block, const IR& ir);

// ============================================================================
// Residuals File Parser
// ============================================================================

struct EESResidualsInfo {
    int totalEquations = 0;
    
    struct EquationEntry {
        int blockId;
        std::string equationText;
        std::string outputVariable;  // The "bold" variable
    };
    std::vector<EquationEntry> entries;
    
    // Block statistics from EES
    std::map<int, std::vector<int>> blockToEquations;  // blockId -> list of equation indices
};

// Parse an EES .residuals file
std::optional<EESResidualsInfo> parseResidualsFile(const std::string& filepath);

// ============================================================================
// Output Generation
// ============================================================================

// Generate a report similar to EES residuals format
std::string generateResidualsReport(const IR& ir, const StructuralAnalysisResult& result);

// Generate JSON output for automated testing
std::string generateAnalysisJSON(const IR& ir, const StructuralAnalysisResult& result);

}  // namespace coolsolve
