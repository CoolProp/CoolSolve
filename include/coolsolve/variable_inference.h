#pragma once
#include "coolsolve/ir.h"

namespace coolsolve {

// Infers variable properties (fluid, type) from usage in thermodynamic function calls
void inferVariables(IR& ir);

// Calculates and sets initial values for variables without guesses
// Uses inferred properties to calculate values at T=25C, x=0
// Defaults to 1.0 for others
void initializeVariables(IR& ir);

} // namespace coolsolve
