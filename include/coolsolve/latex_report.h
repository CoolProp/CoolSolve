#pragma once

/**
 * @file latex_report.h
 * @brief Comprehensive LaTeX report generator for CoolSolve.
 *
 * Generates a self-contained LaTeX document containing:
 *   - Model equations (grouped by solution block)
 *   - Variable table with solved values and units
 *   - Block structure summary and solver statistics
 *   - Timing profile
 *
 * The generated .tex file is compiled to PDF by the GUI frontend using
 * the configured LaTeX compiler (default: pdflatex).  The backend only
 * produces the .tex source — no LaTeX installation is required on the
 * solver side.
 */

#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/solver.h"
#include "coolsolve/runner.h"
#include <string>

namespace coolsolve {

/**
 * @brief Generate a comprehensive LaTeX report as a string.
 *
 * The returned string is a complete .tex document that can be compiled
 * with pdflatex (or another LaTeX engine).  It does NOT require any
 * external files — all content is inline.
 *
 * @param ir           Intermediate representation (equations, variables)
 * @param analysis     Structural analysis result (blocks, matching)
 * @param solveResult  Solver output (variable values, per-block stats)
 * @param timing       Pipeline timing data
 * @param modelName    User-facing model name (used in the title)
 * @param inputFile    Path of the source .eescode file
 * @return Complete LaTeX document as a string
 */
std::string generateLatexReport(
    const IR& ir,
    const StructuralAnalysisResult& analysis,
    const SolveResult& solveResult,
    const CoolSolveRunner::PipelineTiming& timing,
    const std::string& modelName,
    const std::string& inputFile
);

/**
 * @brief Write the LaTeX report to a file.
 *
 * Convenience wrapper that calls generateLatexReport() and writes
 * the result to @p outputPath.
 *
 * @return true if the file was written successfully
 */
bool writeLatexReport(
    const std::string& outputPath,
    const IR& ir,
    const StructuralAnalysisResult& analysis,
    const SolveResult& solveResult,
    const CoolSolveRunner::PipelineTiming& timing,
    const std::string& modelName,
    const std::string& inputFile
);

}  // namespace coolsolve
