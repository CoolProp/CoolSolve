// The badge lines of a chart: where its properties come from (the references of the equation of state of its
// fluid, from CoolProp's BibTeX keys and its bibliography CoolPropBibTeXLibrary.bib, embedded at build time) and
// what drew it (the versions of CoolProp and CoolPlot).
//
// Mirrors CoolPlot (github.com/CoolProp/CoolPlot) at commit 7338f48 (v1.5.0): thermo/references.py
// (latex_to_text, parse_bibtex, cite, fluid_references, blend_name, composition_text, versions_line) and
// engine/engine.py (Engine.source_lines).
// CoolPlot is MIT licensed: Copyright 2014-2019 Jorrit Wronski and the other developers.
// See docs/coolplot_integration.md §5 (phase 4).
#pragma once

#include "coolsolve/diagram/fluid.h"

#include <string>
#include <vector>

namespace coolsolve {
namespace diagram {

/// The BibTeX key of the humid air model behind HAPropsSI (references.HUMID_AIR_KEY).
extern const char* const HUMID_AIR_KEY;

/// How far (in mass fraction) a mixture may be from a predefined blend of CoolProp to carry its name
/// (references.BLEND_TOLERANCE).
constexpr double BLEND_TOLERANCE = 0.005;

/**
 * @brief Plain text from the LaTeX of a BibTeX field (references.latex_to_text): accents resolved to the letters of
 *        Unicode (UTF-8), braces and commands removed, "--" as "-", white space collapsed.
 */
std::string latexToText(const std::string& text);

/// Whether CoolProp's bibliography was embedded (it is not with a CoolProp without its sources): without it every
/// reference reads as its BibTeX key.
bool hasBibliography();

/**
 * @brief One reference as text (references.cite): "Lemmon et al. (2009)" in the "short" style, with the journal,
 *        volume and pages (or the title of a report or thesis) in the "full" style. The key itself when the
 *        bibliography does not have it or names no author.
 */
std::string cite(const std::string& key, const std::string& style = "short");

/**
 * @brief The name of the predefined blend of CoolProp with these components and mass fractions, e.g. "R410A", or ""
 *        (references.blend_name): every mass fraction within BLEND_TOLERANCE, the closest match winning.
 */
std::string blendName(const std::vector<std::string>& components, const std::vector<double>& massFractions);

/**
 * @brief The lines naming the sources of a fluid's properties (references.fluid_references).
 *
 * A pure fluid of the HEOS backend: "R290: Lemmon et al. (2009)". A mixture: its name and composition ("R410A: R32
 * 50.3 %, R125 49.7 % by mass", or "Mixture: ..." without a predefined blend of that composition), the equation of
 * state of every component, then the references of its mixing model. Humid air: ASHRAE RP-1485. Other backends:
 * the name of the backend, without a reference.
 * @param fluid a fluid record of a CoolPlot spec (FluidSpec.to_dict, or {"model": "humid_air", ...})
 * @param style "short" or "full"
 * @throws DiagramError for another style
 */
std::vector<std::string> fluidReferences(const json& fluid, const std::string& style = "short");

/// "CoolProp 8.1.0dev (8eb3b77), CoolPlot 1.5.0" (references.versions_line): the git revision for a development
/// version of CoolProp only; CoolPlot's version is that of coolplot-core (coolplot_version()).
std::string versionsLine();

/**
 * @brief The badge lines of a spec record as its settings.branding asks for them (Engine.source_lines): the
 *        references of its fluid ("references": "short" by default, "full", or "none") and the versions line
 *        ("versions": true by default). A record before CoolPlot 1.2 (settings with "calc" but no "branding") asks
 *        for neither, as Settings.from_dict reads it.
 * @param spec a spec record with its "fluid"
 * @throws DiagramError for a references style CoolPlot does not know
 */
std::vector<std::string> sourceLines(const json& spec);

}  // namespace diagram
}  // namespace coolsolve
