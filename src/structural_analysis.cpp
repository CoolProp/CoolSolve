#include "coolsolve/structural_analysis.h"
#include "coolsolve/constants.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <stack>
#include <regex>
#include <limits>

namespace coolsolve {

// ============================================================================
// Hopcroft-Karp Algorithm for Maximum Bipartite Matching
// ============================================================================

std::vector<int> StructuralAnalyzer::hopcroftKarp(
    int numEquations,
    int numVariables,
    const std::vector<std::set<std::string, CaseInsensitiveLess>>& incidenceMatrix,
    const std::map<std::string, int, CaseInsensitiveLess>& varNameToIndex)
{
    // pairU[eq] = variable matched to equation eq, -1 if unmatched
    // pairV[var] = equation matched to variable var, -1 if unmatched
    std::vector<int> pairU(numEquations, -1);
    std::vector<int> pairV(numVariables, -1);
    std::vector<int> dist(numEquations + 1);
    
    const int NIL = numEquations;  // Dummy node
    const int INF = std::numeric_limits<int>::max();
    
    // Build adjacency list
    std::vector<std::vector<int>> adj(numEquations);
    for (int eq = 0; eq < numEquations; ++eq) {
        for (const auto& varName : incidenceMatrix[eq]) {
            auto it = varNameToIndex.find(varName);
            if (it != varNameToIndex.end()) {
                adj[eq].push_back(it->second);
            }
        }
    }
    
    // BFS to build layers
    auto bfs = [&]() -> bool {
        std::queue<int> Q;
        for (int eq = 0; eq < numEquations; ++eq) {
            if (pairU[eq] == -1) {
                dist[eq] = 0;
                Q.push(eq);
            } else {
                dist[eq] = INF;
            }
        }
        dist[NIL] = INF;
        
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();
            
            if (dist[u] < dist[NIL]) {
                for (int v : adj[u]) {
                    int next = (pairV[v] == -1) ? NIL : pairV[v];
                    if (dist[next] == INF) {
                        dist[next] = dist[u] + 1;
                        if (next != NIL) {
                            Q.push(next);
                        }
                    }
                }
            }
        }
        return dist[NIL] != INF;
    };
    
    // DFS to find augmenting paths
    std::function<bool(int)> dfs = [&](int u) -> bool {
        if (u == NIL) return true;
        
        for (int v : adj[u]) {
            int next = (pairV[v] == -1) ? NIL : pairV[v];
            if (dist[next] == dist[u] + 1 && dfs(next)) {
                pairV[v] = u;
                pairU[u] = v;
                return true;
            }
        }
        dist[u] = INF;
        return false;
    };
    
    // Main Hopcroft-Karp loop
    while (bfs()) {
        for (int eq = 0; eq < numEquations; ++eq) {
            if (pairU[eq] == -1) {
                dfs(eq);
            }
        }
    }
    
    return pairU;  // Returns matching: eq -> var index
}

// ============================================================================
// Tarjan's Algorithm for Strongly Connected Components
// ============================================================================

std::vector<std::vector<int>> StructuralAnalyzer::tarjanSCC(
    int numNodes,
    const std::vector<std::vector<int>>& adjacencyList)
{
    std::vector<int> index(numNodes, -1);
    std::vector<int> lowlink(numNodes);
    std::vector<bool> onStack(numNodes, false);
    std::stack<int> S;
    std::vector<std::vector<int>> SCCs;
    int currentIndex = 0;
    
    std::function<void(int)> strongconnect = [&](int v) {
        index[v] = currentIndex;
        lowlink[v] = currentIndex;
        currentIndex++;
        S.push(v);
        onStack[v] = true;
        
        for (int w : adjacencyList[v]) {
            if (index[w] == -1) {
                strongconnect(w);
                lowlink[v] = std::min(lowlink[v], lowlink[w]);
            } else if (onStack[w]) {
                lowlink[v] = std::min(lowlink[v], index[w]);
            }
        }
        
        if (lowlink[v] == index[v]) {
            std::vector<int> scc;
            int w;
            do {
                w = S.top();
                S.pop();
                onStack[w] = false;
                scc.push_back(w);
            } while (w != v);
            SCCs.push_back(scc);
        }
    };
    
    for (int v = 0; v < numNodes; ++v) {
        if (index[v] == -1) {
            strongconnect(v);
        }
    }
    
    // Note: With our dependency graph definition (adj[eq] = equations that eq depends on),
    // Tarjan outputs SCCs in correct topological order (dependencies first).
    // No reversal needed.
    return SCCs;
}

// ============================================================================
// Build Dependency Graph
// ============================================================================

std::vector<std::vector<int>> StructuralAnalyzer::buildDependencyGraph(
    const IR& ir,
    const std::vector<int>& matching)
{
    const auto& equations = ir.getEquations();
    int numEq = static_cast<int>(equations.size());
    
    // Create reverse mapping: variable name -> equation that determines it
    std::map<std::string, int, CaseInsensitiveLess> varToEq;
    const auto& vars = ir.getVariables();
    
    // IMPORTANT: Build indexToVar excluding constants, matching how Hopcroft-Karp
    // builds its allVars list. The matching indices are into this filtered list.
    std::vector<std::string> indexToVar;
    indexToVar.reserve(vars.size());
    for (const auto& [name, info] : vars) {
        if (!Constants::isConstant(name)) {
            indexToVar.push_back(name);
        }
    }
    
    for (int eq = 0; eq < numEq; ++eq) {
        if (matching[eq] >= 0 && matching[eq] < static_cast<int>(indexToVar.size())) {
            varToEq[indexToVar[matching[eq]]] = eq;
        }
    }
    
    // Build adjacency list: eq -> list of equations it depends on
    std::vector<std::vector<int>> adj(numEq);
    
    for (int eq = 0; eq < numEq; ++eq) {
        // Get variables used by this equation (except its output variable)
        for (const auto& varName : equations[eq].variables) {
            auto it = varToEq.find(varName);
            if (it != varToEq.end() && it->second != eq) {
                // This equation depends on the equation that determines varName
                adj[eq].push_back(it->second);
            }
        }
    }
    
    // ---- CALL output consistency edges ----
    // When a procedure CALL has an output variable that is matched to a
    // different equation (e.g. CALL output DT_min is matched to the
    // constraint equation "DELTAT_pp = DT_min" rather than to the CALL
    // equation itself), the two equations must be in the same SCC so the
    // solver can enforce that the CALL-computed output equals the
    // constrained value.  Add an artificial back-edge from the non-CALL
    // equation to the primary CALL equation to create this cycle.
    for (int eq = 0; eq < numEq; ++eq) {
        if (!equations[eq].procedureCall) continue;
        
        const auto& call = *equations[eq].procedureCall;
        // Only process the primary CALL equation (the first one for this
        // source line — all sibling output equations share the same
        // procedureCall, so deduplicate by checking if the previous
        // equation is the same CALL).
        if (eq > 0 && equations[eq - 1].procedureCall &&
            equations[eq - 1].sourceLine == equations[eq].sourceLine) {
            continue;  // Skip sibling CALL equations; only process the first
        }
        
        // For each output variable of this CALL
        for (const auto& outVar : call.outputVars) {
            std::string outName = outVar.flattenedName();
            auto it = varToEq.find(outName);
            if (it == varToEq.end()) continue;
            
            int otherEq = it->second;
            // If this output is matched to a different equation (not one of
            // the sibling CALL equations for this same source line), add a
            // back-edge from that equation to this CALL equation.
            if (otherEq != eq) {
                bool isSiblingCall = false;
                // Check if otherEq is a sibling CALL equation (same source line)
                if (equations[otherEq].procedureCall &&
                    equations[otherEq].sourceLine == equations[eq].sourceLine) {
                    isSiblingCall = true;
                }
                if (!isSiblingCall) {
                    // Add back-edge: otherEq depends on this CALL
                    adj[otherEq].push_back(eq);
                }
            }
        }
    }
    
    return adj;
}

// ============================================================================
// Tear Set Computation (Greedy Feedback Vertex Set)
// ============================================================================

namespace {

// Build equation dependency graph for a single block (local indices 0..n-1).
// eq i depends on eq j if the output variable of eq j appears in eq i.
std::vector<std::vector<int>> buildBlockDependencyGraph(
    const Block& block,
    const IR& ir)
{
    const auto& equations = ir.getEquations();
    const size_t n = block.equationIds.size();
    std::vector<std::vector<int>> adj(n);

    for (size_t i = 0; i < n; ++i) {
        const int globalEqId = block.equationIds[i];
        if (globalEqId < 0 || globalEqId >= static_cast<int>(equations.size())) continue;
        const auto& varsInEq = equations[static_cast<size_t>(globalEqId)].variables;
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            if (varsInEq.count(block.variables[j])) {
                adj[i].push_back(static_cast<int>(j));
            }
        }
    }
    return adj;
}

// Find any cycle in the graph (only among nodes where removed[i] is false).
std::vector<int> findAnyCycle(const std::vector<std::vector<int>>& adj,
                              const std::vector<bool>& removed) {
    const int n = static_cast<int>(adj.size());
    std::vector<int> color(n, 0);
    std::vector<int> parent(n, -1);
    std::vector<int> cycleNodes;

    std::function<bool(int)> dfs = [&](int v) -> bool {
        if (removed[v]) return false;
        color[v] = 1;
        for (int w : adj[v]) {
            if (removed[w]) continue;
            if (color[w] == 0) {
                parent[w] = v;
                if (dfs(w)) return true;
            } else if (color[w] == 1) {
                cycleNodes.push_back(w);
                for (int u = v; u != w; u = parent[u]) {
                    cycleNodes.push_back(u);
                }
                return true;
            }
        }
        color[v] = 2;
        return false;
    };

    for (int i = 0; i < n; ++i) {
        if (removed[i] || color[i] != 0) continue;
        if (dfs(i)) {
            return cycleNodes;
        }
    }
    return cycleNodes;
}

// Topological sort of non-removed nodes. adj[i] = equations i depends on (j in adj[i] => j before i).
// inDegree[i] = count of (non-removed) deps of i. When we add v, decrement inDegree for i such that v in adj[i].
void topologicalSort(
    const std::vector<std::vector<int>>& adj,
    const std::vector<bool>& removed,
    std::vector<int>& order)
{
    const int n = static_cast<int>(adj.size());
    std::vector<int> inDegree(n, 0);
    std::vector<std::vector<int>> revAdj(n);
    for (int i = 0; i < n; ++i) {
        if (removed[i]) continue;
        for (int j : adj[i]) {
            if (!removed[j]) {
                inDegree[i]++;
                revAdj[j].push_back(i);
            }
        }
    }
    std::queue<int> q;
    for (int i = 0; i < n; ++i) {
        if (!removed[i] && inDegree[i] == 0) q.push(i);
    }
    order.clear();
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        order.push_back(v);
        for (int i : revAdj[v]) {
            if (removed[i]) continue;
            inDegree[i]--;
            if (inDegree[i] == 0) q.push(i);
        }
    }
}

}  // namespace

BlockTearSetResult computeBlockTearSet(const Block& block, const IR& ir) {
    BlockTearSetResult result;
    if (block.size() < 2) {
        return result;
    }

    std::vector<std::vector<int>> adj = buildBlockDependencyGraph(block, ir);
    const int n = static_cast<int>(adj.size());
    std::vector<bool> removed(n, false);

    // Greedy feedback vertex set: repeatedly find a cycle and remove the node with max out-degree in that cycle
    while (true) {
        std::vector<int> cycle = findAnyCycle(adj, removed);
        if (cycle.empty()) break;
        // Build active adjacency (excluding removed)
        std::vector<std::vector<int>> activeAdj(n);
        for (int i = 0; i < n; ++i) {
            if (removed[i]) continue;
            for (int j : adj[i]) {
                if (!removed[j]) activeAdj[i].push_back(j);
            }
        }
        std::vector<int> activeCycle;
        for (int v : cycle) {
            if (!removed[v]) activeCycle.push_back(v);
        }
        if (activeCycle.empty()) break;
        int best = activeCycle[0];
        int bestOut = static_cast<int>(activeAdj[best].size());
        for (size_t i = 1; i < activeCycle.size(); ++i) {
            int v = activeCycle[i];
            int out = static_cast<int>(activeAdj[v].size());
            if (out > bestOut) {
                bestOut = out;
                best = v;
            }
        }
        removed[best] = true;
    }

    // Tear equations = removed nodes. Tear vars = their output variables.
    for (int i = 0; i < n; ++i) {
        if (removed[i]) {
            result.tearEquationIds.push_back(block.equationIds[i]);
            result.tearVarNames.push_back(block.variables[i]);
        }
    }

    // Topological order of non-tear equations (global equation IDs)
    std::vector<int> order;
    topologicalSort(adj, removed, order);
    for (int localIdx : order) {
        result.topoOrderNonTearEqIds.push_back(block.equationIds[localIdx]);
    }

    return result;
}

// ============================================================================
// Re-Decompose a Reduced Block into Sub-SCCs
// ============================================================================

std::vector<Block> StructuralAnalyzer::redecomposeBlock(
    const std::vector<int>& reducedEquationIds,
    const std::vector<std::string>& reducedVariables,
    const IR& ir,
    const StructuralAnalysisResult& analysis)
{
    const size_t k = reducedEquationIds.size();
    if (k <= 1) {
        // Trivial: size 0 or 1 — return a single block
        Block b;
        b.id = 0;
        b.equationIds = reducedEquationIds;
        b.variables = reducedVariables;
        return {b};
    }

    const auto& equations = ir.getEquations();

    // Build local output mapping: localIndex → output variable name
    // using the global matching from the original structural analysis.
    std::vector<std::string> localOutput(k);
    for (size_t i = 0; i < k; ++i) {
        int gid = reducedEquationIds[i];
        if (gid >= 0 && gid < static_cast<int>(analysis.matching.size())) {
            localOutput[i] = analysis.matching[gid];
        }
    }

    // Build reverse map: variable name → local equation index
    std::map<std::string, int, CaseInsensitiveLess> varToLocalIdx;
    for (size_t i = 0; i < k; ++i) {
        if (!localOutput[i].empty()) {
            varToLocalIdx[localOutput[i]] = static_cast<int>(i);
        }
    }

    // Build set of reduced variables for fast lookup
    std::set<std::string, CaseInsensitiveLess> reducedVarSet(
        reducedVariables.begin(), reducedVariables.end());

    // Build local adjacency list: localEq i → localEq j
    // if equation i references the output variable of equation j
    std::vector<std::vector<int>> adj(k);
    for (size_t i = 0; i < k; ++i) {
        int gid = reducedEquationIds[i];
        if (gid < 0 || gid >= static_cast<int>(equations.size())) continue;
        const auto& vars = equations[gid].variables;
        for (const auto& v : vars) {
            // Skip if variable is not in the reduced block
            if (!reducedVarSet.count(v)) continue;
            // Skip if v is the output of equation i itself
            if (CaseInsensitiveLess()(v, localOutput[i]) == false &&
                CaseInsensitiveLess()(localOutput[i], v) == false) continue;
            // Find which local equation produces v
            auto it = varToLocalIdx.find(v);
            if (it != varToLocalIdx.end() && it->second != static_cast<int>(i)) {
                adj[i].push_back(it->second);
            }
        }
    }

    // Run Tarjan's SCC algorithm on the local graph
    auto sccs = tarjanSCC(static_cast<int>(k), adj);

    // Build sub-blocks in topological order
    // tarjanSCC already returns SCCs in topological order (dependencies first)
    std::vector<Block> subBlocks;
    subBlocks.reserve(sccs.size());
    for (size_t si = 0; si < sccs.size(); ++si) {
        Block b;
        b.id = static_cast<int>(si);
        for (int localIdx : sccs[si]) {
            b.equationIds.push_back(reducedEquationIds[localIdx]);
            b.variables.push_back(localOutput[localIdx]);
        }
        subBlocks.push_back(std::move(b));
    }

    return subBlocks;
}

// ============================================================================
// Main Analysis Function
// ============================================================================

StructuralAnalysisResult StructuralAnalyzer::analyze(IR& ir) {
    StructuralAnalysisResult result;
    result.success = true;
    
    const auto& equations = ir.getEquations();
    const auto& variables = ir.getVariables();
    const auto& incidenceMatrix = ir.getIncidenceMatrix();
    
    int numEq = static_cast<int>(equations.size());
    
    // We match against ALL variables (including strings)
    std::vector<std::string> allVars;
    std::map<std::string, int, CaseInsensitiveLess> varNameToIndex;
    for (const auto& [name, info] : variables) {
        // Exclude constants from matching
        if (!Constants::isConstant(name)) {
            varNameToIndex[name] = static_cast<int>(allVars.size());
            allVars.push_back(name);
        }
    }
    
    int numVar = static_cast<int>(allVars.size());
    
    result.totalEquations = numEq;
    result.totalVariables = numVar;
    
    // Check if system is square
    if (numEq != numVar) {
        result.success = false;
        result.errorMessage = "There are " + std::to_string(numEq) + 
                             " equations and " + std::to_string(numVar) + 
                             " unknowns. The system is not square and cannot be solved";
        // Return early - cannot proceed with non-square system
        return result;
    }
    
    // Run Hopcroft-Karp matching
    std::vector<int> matchingIndices = hopcroftKarp(numEq, numVar, incidenceMatrix, varNameToIndex);
    
    // Convert matching to variable names and check for unmatched equations
    result.matching.resize(numEq);
    std::vector<int> unmatchedEquations;
    for (int eq = 0; eq < numEq; ++eq) {
        if (matchingIndices[eq] >= 0 && matchingIndices[eq] < numVar) {
            result.matching[eq] = allVars[matchingIndices[eq]];
        } else {
            // This equation could not be matched to any variable
            unmatchedEquations.push_back(eq);
        }
    }
    
    // Check for structural issues: unmatched equations indicate overdetermination
    if (!unmatchedEquations.empty()) {
        std::ostringstream errMsg;
        errMsg << "Structural singularity detected: " << unmatchedEquations.size() 
               << " equation(s) could not be matched to any variable. ";
        errMsg << "This typically means the system is overdetermined - ";
        errMsg << "some variables are defined by multiple conflicting equations.\n";
        errMsg << "Unmatched equations:\n";
        
        for (int eqId : unmatchedEquations) {
            if (eqId < static_cast<int>(equations.size())) {
                errMsg << "  - Equation " << eqId << ": " << equations[eqId].originalText << "\n";
            }
        }
        
        // Add warning but don't fail - let the solver catch the inconsistency
        // This allows us to still attempt solving and report which equations fail
        result.errorMessage = errMsg.str();
        // Note: we don't set result.success = false here because the system might still
        // be solvable if the redundant equations happen to be consistent
    }
    
    // Update IR with matching information
    auto& eqsMutable = const_cast<std::vector<EquationInfo>&>(ir.getEquations());
    for (int eq = 0; eq < numEq; ++eq) {
        eqsMutable[eq].outputVariable = result.matching[eq];
    }
    
    // Build dependency graph
    auto dependencyGraph = buildDependencyGraph(ir, matchingIndices);
    
    // Run Tarjan's SCC algorithm
    auto sccs = tarjanSCC(numEq, dependencyGraph);
    
    // Create blocks from SCCs
    int blockId = 0;
    for (const auto& scc : sccs) {
        Block block;
        block.id = blockId++;
        block.equationIds = scc;
        
        for (int eqId : scc) {
            if (eqId < numEq && !result.matching[eqId].empty()) {
                block.variables.push_back(result.matching[eqId]);
            }
            
            // Update equation's block ID in IR
            if (eqId < numEq) {
                eqsMutable[eqId].blockId = block.id;
            }
        }
        
        result.blocks.push_back(block);
        
        if (block.size() > static_cast<size_t>(result.largestBlockSize)) {
            result.largestBlockSize = static_cast<int>(block.size());
        }
        
        if (block.isExplicit()) {
            result.explicitEquationCount++;
        }
    }
    
    result.totalBlocks = blockId;
    
    return result;
}

// ============================================================================
// Comparison with EES Residuals
// ============================================================================

StructuralAnalysisResult::ComparisonResult StructuralAnalyzer::compareWithEES(
    const StructuralAnalysisResult& result,
    const std::string& residualsFilePath)
{
    StructuralAnalysisResult::ComparisonResult comparison;
    
    auto eesInfo = parseResidualsFile(residualsFilePath);
    if (!eesInfo) {
        comparison.differences.push_back("Could not parse EES residuals file");
        return comparison;
    }
    
    comparison.eesBlockCount = static_cast<int>(eesInfo->blockToEquations.size());
    comparison.ourBlockCount = result.totalBlocks;
    
    // Compare total equations
    if (eesInfo->totalEquations != result.totalEquations) {
        comparison.differences.push_back(
            "Equation count mismatch: EES=" + std::to_string(eesInfo->totalEquations) +
            ", ours=" + std::to_string(result.totalEquations));
    }
    
    // Check if we have at least as many blocks as EES (more is better - more decomposition)
    if (comparison.ourBlockCount < comparison.eesBlockCount) {
        comparison.differences.push_back(
            "Fewer blocks than EES: EES=" + std::to_string(comparison.eesBlockCount) +
            ", ours=" + std::to_string(comparison.ourBlockCount));
    }
    
    comparison.matches = comparison.differences.empty();
    return comparison;
}

// ============================================================================
// Residuals File Parser
// ============================================================================

std::optional<EESResidualsInfo> parseResidualsFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return std::nullopt;
    }
    
    EESResidualsInfo info;
    std::string line;
    
    // First line: "There are a total of N equations..."
    if (!std::getline(file, line)) {
        return std::nullopt;
    }
    
    std::regex totalRegex(R"(There are a total of (\d+) equations)");
    std::smatch match;
    if (std::regex_search(line, match, totalRegex)) {
        info.totalEquations = std::stoi(match[1].str());
    }
    
    // Skip header line
    std::getline(file, line);
    
    // Parse equation entries
    // Format: " Block\tRel. Res.\tAbs. Res.\tUnits\tCalls\tTime(ms)\tEquations"
    std::regex entryRegex(R"(^\s*(\d+)\s+[\d.E+-]+\s+[\d.E+-]+\s+\S+\s+\d+\s+\d+\s+(.+)$)");
    
    while (std::getline(file, line)) {
        if (line.find("Variables shown in bold") != std::string::npos) {
            break;
        }
        
        if (std::regex_search(line, match, entryRegex)) {
            EESResidualsInfo::EquationEntry entry;
            entry.blockId = std::stoi(match[1].str());
            entry.equationText = match[2].str();
            
            // Try to extract the "bold" output variable from the equation
            // In EES, the output variable is shown in bold, but in text it's typically
            // the variable that appears alone on one side of the equation
            size_t eqPos = entry.equationText.find('=');
            if (eqPos != std::string::npos) {
                std::string lhs = entry.equationText.substr(0, eqPos);
                // Remove spaces
                lhs.erase(std::remove_if(lhs.begin(), lhs.end(), ::isspace), lhs.end());
                entry.outputVariable = lhs;
            }
            
            info.entries.push_back(entry);
            info.blockToEquations[entry.blockId].push_back(
                static_cast<int>(info.entries.size()) - 1);
        }
    }
    
    return info;
}

// ============================================================================
// Output Generation
// ============================================================================

std::string generateResidualsReport(const IR& ir, const StructuralAnalysisResult& result) {
    std::ostringstream oss;
    
    oss << "There are a total of " << result.totalEquations << " equations in the Main program.\n";
    oss << "Block\tRel. Res.\tAbs. Res.\tUnits\tCalls\tTime(ms)\tEquations\n";
    
    // Group equations by block
    std::map<int, std::vector<const EquationInfo*>> blockEquations;
    for (const auto& eq : ir.getEquations()) {
        blockEquations[eq.blockId].push_back(&eq);
    }
    
    for (const auto& [blockId, eqs] : blockEquations) {
        for (const auto* eq : eqs) {
            oss << " " << blockId;
            oss << "\t 0.000E+00";
            oss << "\t 0.000E+00";
            oss << "\t OK";
            oss << "\t1\t0\t";
            oss << eq->originalText << "\n";
        }
    }
    
    oss << "\nVariables shown in bold font are determined by the equation(s) in each block.\n";
    
    return oss.str();
}

std::string generateAnalysisJSON(const IR& ir, const StructuralAnalysisResult& result) {
    nlohmann::json j;
    
    // Statistics
    j["stats"]["totalEquations"] = result.totalEquations;
    j["stats"]["totalVariables"] = result.totalVariables;
    j["stats"]["totalBlocks"] = result.totalBlocks;
    j["stats"]["largestBlockSize"] = result.largestBlockSize;
    j["stats"]["explicitEquationCount"] = result.explicitEquationCount;
    j["stats"]["success"] = result.success;
    
    if (!result.errorMessage.empty()) {
        j["error"] = result.errorMessage;
    }
    
    // Blocks
    nlohmann::json blocksJson = nlohmann::json::array();
    for (const auto& block : result.blocks) {
        nlohmann::json blockJson;
        blockJson["id"] = block.id;
        blockJson["size"] = block.size();
        blockJson["isExplicit"] = block.isExplicit();
        blockJson["equationIds"] = block.equationIds;
        blockJson["variables"] = block.variables;
        blocksJson.push_back(blockJson);
    }
    j["blocks"] = blocksJson;
    
    // Matching
    nlohmann::json matchingJson = nlohmann::json::array();
    const auto& equations = ir.getEquations();
    for (size_t i = 0; i < result.matching.size(); ++i) {
        nlohmann::json m;
        m["equationId"] = static_cast<int>(i);
        m["equationText"] = (i < equations.size()) ? equations[i].originalText : "";
        m["outputVariable"] = result.matching[i];
        m["blockId"] = (i < equations.size()) ? equations[i].blockId : -1;
        matchingJson.push_back(m);
    }
    j["matching"] = matchingJson;
    
    // EES Comparison (if available)
    if (result.eesComparison) {
        nlohmann::json compJson;
        compJson["matches"] = result.eesComparison->matches;
        compJson["eesBlockCount"] = result.eesComparison->eesBlockCount;
        compJson["ourBlockCount"] = result.eesComparison->ourBlockCount;
        compJson["differences"] = result.eesComparison->differences;
        j["eesComparison"] = compJson;
    }
    
    return j.dump(2);
}

}  // namespace coolsolve
