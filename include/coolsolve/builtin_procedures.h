#pragma once

#include "coolsolve/psychprops.h"

namespace coolsolve {

/**
 * @file builtin_procedures.h
 * @brief The built-in procedures whose outputs stand at fixed positions after the colon and may be skipped.
 *
 * EES files comment outputs out (`CALL NASA(X$, T : {cp}, h, {s})`, `CALL psychprops(... : , , h_1, , , w_1)`): the
 * output of a position that is not written is not computed. The parser records the position of each output written
 * (`ProcedureCall::outputSlots`, empty when they are the first ones in order) and refuses a position that does not exist;
 * the evaluator computes the output of the position. One description per procedure, read by both, so that every
 * such procedure maps its outputs the same way: CALL psychprops and CALL NASA.
 */
struct PositionedOutputs {
    const char* procedure;   ///< the name the messages give it
    int count;               ///< the number of outputs
    const char* list;        ///< their names in order, as the messages list them
};

/// CALL psychprops(... : T, v, h, s, u, w, R, WetB, DewP): see psychprops.h
inline constexpr PositionedOutputs kPsychpropsOutputs = {"psychprops", psychprops::kOutputCount, "T, v, h, s, u, w, R, WetB, DewP"};

/// CALL NASA(Species$, T : cp, h, s): the molar specific heat, enthalpy and entropy of a species of the NASA-Glenn data
inline constexpr PositionedOutputs kNasaOutputs = {"NASA", 3, "cp, h, s"};

}  // namespace coolsolve
