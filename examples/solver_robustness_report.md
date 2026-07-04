# CoolSolve Solver Robustness Report

**Generated:** Fri Jul  3 14:42:04 2026

Total example files tested: 42

**Legend:** OK = converged, PARSE = parse error, ANALYSIS = structural analysis error.  
Failure cells show: `ErrorCategory blkN |F|=residual` where N is the failed block size.

## Summary: Solve Success Rate by Configuration

| # | Configuration | Initials | Tearing | SymReduc | Solved | Total | Rate | Avg time (s) |
|---:|---|:---:|:---:|:---:|---:|---:|---:|---:|
| 1 | Default pipeline (with initials) | Yes | No | No | 38 | 41 | 92.7% | 0.027 |
| 2 | Newton only (with initials) | Yes | No | No | 34 | 41 | 82.9% | 0.011 |
| 3 | TrustRegion only (with initials) | Yes | No | No | 37 | 41 | 90.2% | 0.050 |
| 4 | TrustRegion + Hybrd K=5 (with initials) | Yes | No | No | 37 | 41 | 90.2% | 0.423 |
| 5 | LevenbergMarquardt only (with initials) | Yes | No | No | 33 | 41 | 80.5% | 0.366 |
| 6 | BisectionND only (with initials) | Yes | No | No | 17 | 41 | 41.5% | 0.003 |
| 7 | Homotopy only (with initials) | Yes | No | No | 35 | 41 | 85.4% | 0.036 |
| 8 | Partitioned only (with initials) | Yes | No | No | 16 | 41 | 39.0% | 0.011 |
| 9 | Default + Tearing (with initials) | Yes | Yes | No | 38 | 41 | 92.7% | 0.027 |
| 10 | Default + SymbolicReduction (with initials) | Yes | No | Yes | 38 | 41 | 92.7% | 0.027 |
| 11 | Default + Tearing + SymbolicReduction (with initials) | Yes | Yes | Yes | 38 | 41 | 92.7% | 0.027 |
| 12 | Default pipeline (NO initials) | No | No | No | 30 | 40 | 75.0% | 0.148 |
| 13 | Newton only (NO initials) | No | No | No | 24 | 40 | 60.0% | 0.278 |
| 14 | TrustRegion only (NO initials) | No | No | No | 26 | 40 | 65.0% | 1.295 |
| 15 | TrustRegion + Hybrd K=5 (NO initials) | No | No | No | 22 | 40 | 55.0% | 0.022 |
| 16 | LevenbergMarquardt only (NO initials) | No | No | No | 19 | 40 | 47.5% | 0.014 |
| 17 | BisectionND only (NO initials) | No | No | No | 14 | 40 | 35.0% | 0.004 |
| 18 | Homotopy only (NO initials) | No | No | No | 27 | 40 | 67.5% | 0.283 |
| 19 | Partitioned only (NO initials) | No | No | No | 14 | 40 | 35.0% | 0.015 |
| 20 | Default + Tearing (NO initials) | No | Yes | No | 30 | 40 | 75.0% | 0.146 |
| 21 | Default + SymbolicReduction (NO initials) | No | No | Yes | 30 | 40 | 75.0% | 0.146 |
| 22 | Default + Tearing + SymbolicReduction (NO initials) | No | Yes | Yes | 30 | 40 | 75.0% | 0.141 |
| 23 | Default + MultiStart parallel (NO initials) | No | No | No | 32 | 40 | 80.0% | 0.584 |
| 24 | Default + MultiStart parallel (with initials) | Yes | No | No | 38 | 41 | 92.7% | 0.028 |
| 25 | KINSOL linesearch (with initials) | Yes | No | No | 35 | 41 | 85.4% | 0.015 |
| 26 | KINSOL picard (with initials) | Yes | No | No | 15 | 41 | 36.6% | 0.003 |
| 27 | KINSOL fp/Anderson (with initials) | Yes | No | No | 26 | 41 | 63.4% | 0.010 |
| 28 | KINSOL linesearch (NO initials) | No | No | No | 22 | 40 | 55.0% | 0.036 |
| 29 | KINSOL picard (NO initials) | No | No | No | 10 | 40 | 25.0% | 0.004 |
| 30 | KINSOL fp/Anderson (NO initials) | No | No | No | 20 | 40 | 50.0% | 0.005 |

## Symbolic Reduction Impact

For configurations with symbolic reduction enabled, shows how multi-variable blocks were reduced.

| File | Config | Blocks | MultiVar Blocks (orig→eff) | Variables (orig→reduced) | Reductions Applied | Re-decomposition |
|---|---|---:|---|---|---:|---|
| advanced_features.eescode | Default + SymbolicReduction (with initials) | 48 | 1 → 1 | 2 → 2 | 0 | — |
| air_screw_compressor.eescode | Default + SymbolicReduction (with initials) | 28 | 1 → 1 | 13 → 13 | 0 | — |
| air_screw_compressor_simple.eescode | Default + SymbolicReduction (with initials) | 29 | 1 → 1 | 13 → 13 | 0 | — |
| condenser_3zones.eescode | Default + SymbolicReduction (with initials) | 50 | 1 → 1 | 62 → 62 | 0 | — |
| condenser_wet.eescode | Default + SymbolicReduction (with initials) | 20 | 1 → 1 | 2 → 2 | 0 | — |
| cooling_coil.eescode | Default + SymbolicReduction (with initials) | 49 | 1 → 1 | 12 → 12 | 0 | — |
| cooling_tower.eescode | Default + SymbolicReduction (with initials) | 34 | 2 → 2 | 22 → 22 | 0 | — |
| cooling_tower2.eescode | Default + SymbolicReduction (with initials) | 14 | 1 → 1 | 11 → 11 | 0 | — |
| cpbar.eescode | Default + SymbolicReduction (with initials) | 16 | 1 → 1 | 5 → 5 | 0 | — |
| evaporator.eescode | Default + SymbolicReduction (with initials) | 32 | 1 → 1 | 2 → 2 | 0 | — |
| exchangers2.eescode | Default + SymbolicReduction (with initials) | 26 | 1 → 1 | 4 → 4 | 0 | — |
| exchangers3.eescode | Default + SymbolicReduction (with initials) | 33 | 1 → 1 | 3 → 3 | 0 | — |
| heat_pump_MSTh_SB_R10.eescode | Default + SymbolicReduction (with initials) | 17 | 1 → 1 | 39 → 39 | 0 | — |
| humidair2.eescode | Default + SymbolicReduction (with initials) | 33 | 1 → 1 | 5 → 5 | 0 | — |
| internal_combustion_engine.eescode | Default + SymbolicReduction (with initials) | 41 | 2 → 2 | 7 → 7 | 0 | — |
| internal_combustion_engine_cpbar.eescode | Default + SymbolicReduction (with initials) | 49 | 5 → 5 | 25 → 25 | 0 | — |
| orc_co2.eescode | Default + SymbolicReduction (with initials) | 112 | 1 → 1 | 28 → 28 | 0 | — |
| orc_extraction.eescode | Default + SymbolicReduction (with initials) | 113 | 1 → 1 | 21 → 21 | 0 | — |
| orc_r245fa.eescode | Default + SymbolicReduction (with initials) | 152 | 8 → 8 | 38 → 38 | 0 | — |
| orc_simple.eescode | Default + SymbolicReduction (with initials) | 150 | 7 → 7 | 29 → 29 | 0 | — |
| piston_compressor.eescode | Default + SymbolicReduction (with initials) | 58 | 2 → 2 | 6 → 6 | 0 | — |
| rankine2.eescode | Default + SymbolicReduction (with initials) | 42 | 1 → 1 | 4 → 4 | 0 | — |
| refrigeration_compressor.eescode | Default + SymbolicReduction (with initials) | 53 | 2 → 2 | 6 → 6 | 0 | — |
| scroll_compressor.eescode | Default + SymbolicReduction (with initials) | 66 | 1 → 1 | 34 → 34 | 0 | — |
| turbocompressor.eescode | Default + SymbolicReduction (with initials) | 24 | 1 → 1 | 9 → 9 | 0 | — |
| zorlu_heat_pump.eescode | Default + SymbolicReduction (with initials) | 64 | 1 → 1 | 63 → 63 | 0 | — |
| advanced_features.eescode | Default + Tearing + SymbolicReduction (with initials) | 48 | 1 → 1 | 2 → 2 | 0 | — |
| air_screw_compressor.eescode | Default + Tearing + SymbolicReduction (with initials) | 28 | 1 → 1 | 13 → 13 | 0 | — |
| air_screw_compressor_simple.eescode | Default + Tearing + SymbolicReduction (with initials) | 29 | 1 → 1 | 13 → 13 | 0 | — |
| condenser_3zones.eescode | Default + Tearing + SymbolicReduction (with initials) | 50 | 1 → 1 | 62 → 62 | 0 | — |
| condenser_wet.eescode | Default + Tearing + SymbolicReduction (with initials) | 20 | 1 → 1 | 2 → 2 | 0 | — |
| cooling_coil.eescode | Default + Tearing + SymbolicReduction (with initials) | 49 | 1 → 1 | 12 → 12 | 0 | — |
| cooling_tower.eescode | Default + Tearing + SymbolicReduction (with initials) | 34 | 2 → 2 | 22 → 22 | 0 | — |
| cooling_tower2.eescode | Default + Tearing + SymbolicReduction (with initials) | 14 | 1 → 1 | 11 → 11 | 0 | — |
| cpbar.eescode | Default + Tearing + SymbolicReduction (with initials) | 16 | 1 → 1 | 5 → 5 | 0 | — |
| evaporator.eescode | Default + Tearing + SymbolicReduction (with initials) | 32 | 1 → 1 | 2 → 2 | 0 | — |
| exchangers2.eescode | Default + Tearing + SymbolicReduction (with initials) | 26 | 1 → 1 | 4 → 4 | 0 | — |
| exchangers3.eescode | Default + Tearing + SymbolicReduction (with initials) | 33 | 1 → 1 | 3 → 3 | 0 | — |
| heat_pump_MSTh_SB_R10.eescode | Default + Tearing + SymbolicReduction (with initials) | 17 | 1 → 1 | 39 → 39 | 0 | — |
| humidair2.eescode | Default + Tearing + SymbolicReduction (with initials) | 33 | 1 → 1 | 5 → 5 | 0 | — |
| internal_combustion_engine.eescode | Default + Tearing + SymbolicReduction (with initials) | 41 | 2 → 2 | 7 → 7 | 0 | — |
| internal_combustion_engine_cpbar.eescode | Default + Tearing + SymbolicReduction (with initials) | 49 | 5 → 5 | 25 → 25 | 0 | — |
| orc_co2.eescode | Default + Tearing + SymbolicReduction (with initials) | 112 | 1 → 1 | 28 → 28 | 0 | — |
| orc_extraction.eescode | Default + Tearing + SymbolicReduction (with initials) | 113 | 1 → 1 | 21 → 21 | 0 | — |
| orc_r245fa.eescode | Default + Tearing + SymbolicReduction (with initials) | 152 | 8 → 8 | 38 → 38 | 0 | — |
| orc_simple.eescode | Default + Tearing + SymbolicReduction (with initials) | 150 | 7 → 7 | 29 → 29 | 0 | — |
| piston_compressor.eescode | Default + Tearing + SymbolicReduction (with initials) | 58 | 2 → 2 | 6 → 6 | 0 | — |
| rankine2.eescode | Default + Tearing + SymbolicReduction (with initials) | 42 | 1 → 1 | 4 → 4 | 0 | — |
| refrigeration_compressor.eescode | Default + Tearing + SymbolicReduction (with initials) | 53 | 2 → 2 | 6 → 6 | 0 | — |
| scroll_compressor.eescode | Default + Tearing + SymbolicReduction (with initials) | 66 | 1 → 1 | 34 → 34 | 0 | — |
| turbocompressor.eescode | Default + Tearing + SymbolicReduction (with initials) | 24 | 1 → 1 | 9 → 9 | 0 | — |
| zorlu_heat_pump.eescode | Default + Tearing + SymbolicReduction (with initials) | 64 | 1 → 1 | 63 → 63 | 0 | — |
| advanced_features.eescode | Default + SymbolicReduction (NO initials) | 48 | 1 → 1 | 2 → 2 | 0 | — |
| air_screw_compressor.eescode | Default + SymbolicReduction (NO initials) | 28 | 1 → 1 | 13 → 13 | 0 | — |
| air_screw_compressor_simple.eescode | Default + SymbolicReduction (NO initials) | 29 | 1 → 1 | 13 → 13 | 0 | — |
| condenser_3zones.eescode | Default + SymbolicReduction (NO initials) | 50 | 1 → 1 | 62 → 62 | 0 | — |
| condenser_wet.eescode | Default + SymbolicReduction (NO initials) | 20 | 1 → 1 | 2 → 2 | 0 | — |
| cooling_coil.eescode | Default + SymbolicReduction (NO initials) | 49 | 1 → 1 | 12 → 12 | 0 | — |
| cooling_tower.eescode | Default + SymbolicReduction (NO initials) | 34 | 2 → 2 | 22 → 22 | 0 | — |
| cooling_tower2.eescode | Default + SymbolicReduction (NO initials) | 14 | 1 → 1 | 11 → 11 | 0 | — |
| cpbar.eescode | Default + SymbolicReduction (NO initials) | 16 | 1 → 1 | 5 → 5 | 0 | — |
| evaporator.eescode | Default + SymbolicReduction (NO initials) | 32 | 1 → 1 | 2 → 2 | 0 | — |
| exchangers2.eescode | Default + SymbolicReduction (NO initials) | 26 | 1 → 1 | 4 → 4 | 0 | — |
| exchangers3.eescode | Default + SymbolicReduction (NO initials) | 33 | 1 → 1 | 3 → 3 | 0 | — |
| heat_pump_MSTh_SB_R10.eescode | Default + SymbolicReduction (NO initials) | 17 | 1 → 1 | 39 → 39 | 0 | — |
| humidair2.eescode | Default + SymbolicReduction (NO initials) | 33 | 1 → 1 | 5 → 5 | 0 | — |
| internal_combustion_engine.eescode | Default + SymbolicReduction (NO initials) | 41 | 2 → 2 | 7 → 7 | 0 | — |
| internal_combustion_engine_cpbar.eescode | Default + SymbolicReduction (NO initials) | 49 | 2 → 2 | 12 → 12 | 0 | — |
| orc_co2.eescode | Default + SymbolicReduction (NO initials) | 112 | 1 → 1 | 28 → 28 | 0 | — |
| orc_extraction.eescode | Default + SymbolicReduction (NO initials) | 113 | 1 → 1 | 21 → 21 | 0 | — |
| orc_r245fa.eescode | Default + SymbolicReduction (NO initials) | 152 | 8 → 8 | 38 → 38 | 0 | — |
| orc_simple.eescode | Default + SymbolicReduction (NO initials) | 150 | 7 → 7 | 29 → 29 | 0 | — |
| piston_compressor.eescode | Default + SymbolicReduction (NO initials) | 58 | 2 → 2 | 6 → 6 | 0 | — |
| rankine2.eescode | Default + SymbolicReduction (NO initials) | 42 | 1 → 1 | 4 → 4 | 0 | — |
| refrigeration_compressor.eescode | Default + SymbolicReduction (NO initials) | 53 | 2 → 2 | 6 → 6 | 0 | — |
| scroll_compressor.eescode | Default + SymbolicReduction (NO initials) | 66 | 1 → 1 | 34 → 34 | 0 | — |
| turbocompressor.eescode | Default + SymbolicReduction (NO initials) | 24 | 1 → 1 | 9 → 9 | 0 | — |
| zorlu_heat_pump.eescode | Default + SymbolicReduction (NO initials) | 64 | 1 → 1 | 63 → 63 | 0 | — |
| advanced_features.eescode | Default + Tearing + SymbolicReduction (NO initials) | 48 | 1 → 1 | 2 → 2 | 0 | — |
| air_screw_compressor.eescode | Default + Tearing + SymbolicReduction (NO initials) | 28 | 1 → 1 | 13 → 13 | 0 | — |
| air_screw_compressor_simple.eescode | Default + Tearing + SymbolicReduction (NO initials) | 29 | 1 → 1 | 13 → 13 | 0 | — |
| condenser_3zones.eescode | Default + Tearing + SymbolicReduction (NO initials) | 50 | 1 → 1 | 62 → 62 | 0 | — |
| condenser_wet.eescode | Default + Tearing + SymbolicReduction (NO initials) | 20 | 1 → 1 | 2 → 2 | 0 | — |
| cooling_coil.eescode | Default + Tearing + SymbolicReduction (NO initials) | 49 | 1 → 1 | 12 → 12 | 0 | — |
| cooling_tower.eescode | Default + Tearing + SymbolicReduction (NO initials) | 34 | 2 → 2 | 22 → 22 | 0 | — |
| cooling_tower2.eescode | Default + Tearing + SymbolicReduction (NO initials) | 14 | 1 → 1 | 11 → 11 | 0 | — |
| cpbar.eescode | Default + Tearing + SymbolicReduction (NO initials) | 16 | 1 → 1 | 5 → 5 | 0 | — |
| evaporator.eescode | Default + Tearing + SymbolicReduction (NO initials) | 32 | 1 → 1 | 2 → 2 | 0 | — |
| exchangers2.eescode | Default + Tearing + SymbolicReduction (NO initials) | 26 | 1 → 1 | 4 → 4 | 0 | — |
| exchangers3.eescode | Default + Tearing + SymbolicReduction (NO initials) | 33 | 1 → 1 | 3 → 3 | 0 | — |
| heat_pump_MSTh_SB_R10.eescode | Default + Tearing + SymbolicReduction (NO initials) | 17 | 1 → 1 | 39 → 39 | 0 | — |
| humidair2.eescode | Default + Tearing + SymbolicReduction (NO initials) | 33 | 1 → 1 | 5 → 5 | 0 | — |
| internal_combustion_engine.eescode | Default + Tearing + SymbolicReduction (NO initials) | 41 | 2 → 2 | 7 → 7 | 0 | — |
| internal_combustion_engine_cpbar.eescode | Default + Tearing + SymbolicReduction (NO initials) | 49 | 2 → 2 | 12 → 12 | 0 | — |
| orc_co2.eescode | Default + Tearing + SymbolicReduction (NO initials) | 112 | 1 → 1 | 28 → 28 | 0 | — |
| orc_extraction.eescode | Default + Tearing + SymbolicReduction (NO initials) | 113 | 1 → 1 | 21 → 21 | 0 | — |
| orc_r245fa.eescode | Default + Tearing + SymbolicReduction (NO initials) | 152 | 8 → 8 | 38 → 38 | 0 | — |
| orc_simple.eescode | Default + Tearing + SymbolicReduction (NO initials) | 150 | 7 → 7 | 29 → 29 | 0 | — |
| piston_compressor.eescode | Default + Tearing + SymbolicReduction (NO initials) | 58 | 2 → 2 | 6 → 6 | 0 | — |
| rankine2.eescode | Default + Tearing + SymbolicReduction (NO initials) | 42 | 1 → 1 | 4 → 4 | 0 | — |
| refrigeration_compressor.eescode | Default + Tearing + SymbolicReduction (NO initials) | 53 | 2 → 2 | 6 → 6 | 0 | — |
| scroll_compressor.eescode | Default + Tearing + SymbolicReduction (NO initials) | 66 | 1 → 1 | 34 → 34 | 0 | — |
| turbocompressor.eescode | Default + Tearing + SymbolicReduction (NO initials) | 24 | 1 → 1 | 9 → 9 | 0 | — |
| zorlu_heat_pump.eescode | Default + Tearing + SymbolicReduction (NO initials) | 64 | 1 → 1 | 63 → 63 | 0 | — |

## Detailed Results: With Initials

| File | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt | TR | TR | LM | BisectionND | Homotopy | Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear | Nwt+TR+LM+BisectionND+Homotopy+Part | Kin(LS) | Kin(Pic) | Kin(FP) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| advanced_features.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 |F|=3e+90 | **OK** (0.00s) |
| air_screw_compressor.eescode | **OK** (0.03s) | Singular Jacobian blk13 |F|=322614.1 | **OK** (0.02s) | **OK** (0.02s) | Singular Jacobian blk13 |F|=365861.1 | Other blk13 | **OK** (0.63s) | Max iterations blk13 |F|=2e+08 | **OK** (0.02s) | **OK** (0.05s) | **OK** (0.02s) | **OK** (0.03s) | Other blk13 |F|=177866.5 | Other blk13 |F|=inf | Other blk13 |F|=inf |
| air_screw_compressor_simple.eescode | **OK** (0.03s) | Singular Jacobian blk13 |F|=304355.9 | **OK** (0.02s) | **OK** (0.02s) | Singular Jacobian blk13 |F|=365861.1 | Other blk13 | Other blk13 |F|=419347.5 | Max iterations blk13 |F|=2e+08 | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.03s) | Other blk13 |F|=22785.5 | Other blk13 |F|=inf | Other blk13 |F|=inf |
| boiler_cpbar.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| boiler_cpbar2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.02s) |
| compressor_refrigeration_simple.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| condenser_3zones.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.04s) | **OK** (0.22s) | Other blk62 | **OK** (0.03s) | Max iterations blk62 |F|=3335.5 | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.06s) | Other blk62 |F|=inf | Other blk62 |F|=inf |
| condenser_wet.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.07s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 |F|=inf | **OK** (0.00s) |
| cooling_coil.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk12 | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| cooling_tower.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (7.54s) | Other blk11 | **OK** (0.06s) | Max iterations blk11 |F|=27486.9 | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.03s) | Other blk11 | Other blk11 |
| cooling_tower2.eescode | **OK** (0.57s) | Singular Jacobian blk11 |F|=2e+07 | Max iterations blk11 |F|=2e+07 | Max iterations blk11 |F|=2e+07 | Singular Jacobian blk11 |F|=2e+07 | Other blk11 | Other blk11 |F|=2e+07 | **OK** (0.12s) | **OK** (0.58s) | **OK** (0.58s) | **OK** (0.59s) | **OK** (0.60s) | Singular Jacobian blk11 |F|=2e+07 | Other blk11 | Other blk11 |
| cpbar.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk5 |F|=4e+08 | **OK** (0.02s) | Max iterations blk5 |F|=5e+08 | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.01s) | Other blk5 |F|=inf | **OK** (0.01s) |
| evaporator.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.05s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 |F|=inf | **OK** (0.00s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 |F|=355.6 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.11s) | Other blk4 |F|=inf | **OK** (0.00s) |
| exchangers3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.01s) | Other blk3 |F|=123803.1 | **OK** (0.01s) | Other blk3 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk3 |F|=inf | **OK** (0.00s) |
| heat_pump_MSTh_SB_R10.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.74s) | Other blk39 | **OK** (0.01s) | Max iterations blk39 |F|=2e+158 | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.00s) | Other blk39 |F|=inf | Other blk39 |F|=inf |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.05s) | Other blk5 |F|=2.0 | **OK** (0.01s) | Max iterations blk5 |F|=1e+09 | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.01s) | Other blk5 |F|=inf | **OK** (0.01s) |
| internal_combustion_engine.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.01s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| internal_combustion_engine_cpbar.eescode | **OK** (0.04s) | **OK** (0.04s) | **OK** (1.45s) | **OK** (15.23s) | **OK** (0.04s) | Other blk2 |F|=0.0 | **OK** (0.12s) | Other blk2 | **OK** (0.04s) | **OK** (0.04s) | **OK** (0.04s) | **OK** (0.04s) | **OK** (0.03s) | Other blk2 |F|=inf | **OK** (0.13s) |
| lookup_demo.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| orc_co2.eescode | **OK** (0.03s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | Other blk28 | **OK** (0.03s) | Max iterations blk28 |F|=1e+14 | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.01s) | Other blk28 |F|=inf | Other blk28 |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Max iterations blk21 |F|=0.0 | Other blk21 | **OK** (0.02s) | Max iterations blk21 |F|=5e+47 | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) | Other blk21 |F|=inf | Other blk21 |F|=inf |
| orc_r245fa.eescode | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (1.98s) | Other blk8 |F|=904.1 | **OK** (0.06s) | Other blk2 | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | Other blk8 |F|=inf | Other blk8 |F|=6e+104 |
| orc_simple.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | Other blk7 |F|=0.0 | **OK** (0.03s) | Other blk2 | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.02s) | Other blk7 | **OK** (0.04s) |
| piston_compressor.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| pressuredrop.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| rankine1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| rankine2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 |F|=42594.2 | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 |F|=inf | **OK** (0.00s) |
| refrigeration1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration_compressor.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 |F|=0.1 | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 |F|=inf | **OK** (0.00s) |
| scroll_compressor.eescode | **OK** (0.11s) | **OK** (0.10s) | **OK** (0.09s) | **OK** (0.10s) | Max iterations blk34 |F|=0.8 | Other blk34 | **OK** (0.11s) | Max iterations blk34 |F|=inf | **OK** (0.11s) | **OK** (0.10s) | **OK** (0.10s) | **OK** (0.10s) | **OK** (0.10s) | Other blk34 |F|=inf | Other blk34 |
| simple_centrifugal_compressor.eescode | **OK** (0.00s) | Singular Jacobian blk1 |F|=0.3 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk1 |F|=0.3 | Other blk1 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk1 |F|=0.3 | Other blk1 |F|=6e+15 |
| storage_integraltable.eescode | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS |
| turbocompressor.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk9 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| turbocompressor_interpolate.eescode | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 |
| water_libr.eescode | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 |
| zorlu_heat_pump.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (1.25s) | Other blk63 | **OK** (0.05s) | Max iterations blk63 |F|=4576.0 | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.02s) | Other blk63 |F|=inf | Other blk63 |F|=inf |

## Detailed Results: Without Initials

| File | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt | TR | TR | LM | BisectionND | Homotopy | Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear | Nwt+TR+LM+BisectionND+Homotopy+Part | Kin(LS) | Kin(Pic) | Kin(FP) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| advanced_features.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 |F|=3e+90 | **OK** (0.00s) |
| air_screw_compressor.eescode | **OK** (0.03s) | Singular Jacobian blk13 |F|=322614.1 | **OK** (0.02s) | **OK** (0.02s) | Singular Jacobian blk13 |F|=365861.1 | Other blk13 | **OK** (0.75s) | Max iterations blk13 |F|=2e+08 | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.03s) | Other blk13 |F|=177866.5 | Other blk13 |F|=inf | Other blk13 |F|=inf |
| air_screw_compressor_simple.eescode | **OK** (0.03s) | Singular Jacobian blk13 |F|=304355.9 | **OK** (0.02s) | **OK** (0.02s) | Singular Jacobian blk13 |F|=365861.1 | Other blk13 | Other blk13 |F|=419347.5 | Max iterations blk13 |F|=2e+08 | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.02s) | **OK** (0.03s) | Other blk13 |F|=22785.5 | Other blk13 |F|=inf | Other blk13 |F|=inf |
| boiler_cpbar.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| boiler_cpbar2.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) |
| compressor_refrigeration_simple.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| condenser_3zones.eescode | **OK** (0.09s) | **OK** (0.08s) | Other blk62 |F|=19852.0 | Max iterations blk62 |F|=19288.0 | Max iterations blk62 |F|=921.2 | Other blk62 | **OK** (2.44s) | Max iterations blk62 |F|=1e+07 | **OK** (0.09s) | **OK** (0.09s) | **OK** (0.11s) | **OK** (0.09s) | Singular Jacobian blk62 |F|=254477.7 | Other blk62 |F|=inf | Other blk62 |F|=inf |
| condenser_wet.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.07s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 |F|=inf | **OK** (0.00s) |
| cooling_coil.eescode | Other blk12 |F|=468.5 | **OK** (1.09s) | **OK** (11.29s) | Max iterations blk12 |F|=1.6 | Max iterations blk12 |F|=58.0 | Other blk12 | **OK** (1.83s) | Max iterations blk12 |F|=4e+08 | Other blk12 |F|=468.5 | Other blk12 |F|=468.5 | Other blk12 |F|=468.5 | Other blk12 |F|=468.5 | **OK** (0.44s) | Other blk12 |F|=inf | Other blk12 |F|=8e+16 |
| cooling_tower.eescode | Max iterations blk11 |F|=8e+18 | Singular Jacobian blk11 |F|=2e+07 | Max iterations blk11 |F|=2e+07 | Max iterations blk11 |F|=2e+07 | Singular Jacobian blk11 |F|=2e+07 | Other blk11 | Other blk11 |F|=2e+07 | Max iterations blk11 |F|=8e+18 | Max iterations blk11 |F|=8e+18 | Max iterations blk11 |F|=8e+18 | Max iterations blk11 |F|=8e+18 | Max iterations blk11 |F|=8e+18 | Singular Jacobian blk11 |F|=2e+07 | Other blk11 | Other blk11 |
| cooling_tower2.eescode | **OK** (0.56s) | Singular Jacobian blk11 |F|=2e+07 | Max iterations blk11 |F|=2e+07 | Max iterations blk11 |F|=2e+07 | Singular Jacobian blk11 |F|=2e+07 | Other blk11 | Other blk11 |F|=2e+07 | **OK** (0.12s) | **OK** (0.57s) | **OK** (0.57s) | **OK** (0.58s) | **OK** (0.64s) | Singular Jacobian blk11 |F|=2e+07 | Other blk11 | Other blk11 |
| cpbar.eescode | **OK** (0.00s) | **OK** (0.00s) | Max iterations blk5 |F|=1e+08 | Max iterations blk5 |F|=1e+08 | **OK** (0.01s) | Other blk5 |F|=1e+08 | **OK** (0.02s) | Max iterations blk5 |F|=1e+08 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | Other blk5 |F|=3e+97 | **OK** (0.01s) |
| evaporator.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.05s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 |F|=inf | **OK** (0.00s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 |F|=4e+12 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Singular Jacobian blk4 |F|=272048.6 | Other blk4 |F|=inf | **OK** (0.00s) |
| exchangers3.eescode | **OK** (0.00s) | **OK** (0.00s) | Max iterations blk3 |F|=10003.6 | Max iterations blk3 |F|=5643.2 | Max iterations blk3 |F|=8.1 | Other blk3 |F|=4e+09 | **OK** (0.00s) | Other blk3 | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.00s) | Other blk3 |F|=inf | **OK** (0.00s) |
| heat_pump_MSTh_SB_R10.eescode | Max iterations blk39 |F|=5e+06 | Max iterations blk39 |F|=808738.9 | Max iterations blk39 |F|=68904.0 | Max iterations blk39 |F|=47705.3 | Max iterations blk39 |F|=1674.7 | Other blk39 | Other blk39 |F|=3e+06 | Max iterations blk39 |F|=3e+82 | Max iterations blk39 |F|=5e+06 | Max iterations blk39 |F|=5e+06 | Max iterations blk39 |F|=5e+06 | **OK** (12.16s) | Other blk39 |F|=808669.2 | Other blk39 |F|=inf | Other blk39 |F|=inf |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | **OK** (0.01s) | **OK** (0.01s) | Max iterations blk5 |F|=446.7 | **OK** (0.01s) | Max iterations blk5 |F|=0.0 | Other blk5 |F|=9930.0 | **OK** (0.02s) | Max iterations blk5 |F|=136222.2 | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.01s) | Other blk5 |F|=inf | Other blk5 |F|=inf |
| internal_combustion_engine.eescode | **OK** (0.18s) | Singular Jacobian blk5 |F|=399999.0 | Max iterations blk5 |F|=399999.0 | Max iterations blk5 |F|=399999.0 | Singular Jacobian blk5 |F|=399999.0 | Other blk5 |F|=4198.3 | Other blk5 |F|=399999.0 | Other blk2 | **OK** (0.01s) | **OK** (0.18s) | **OK** (0.01s) | **OK** (0.20s) | Singular Jacobian blk5 |F|=399999.0 | Other blk2 |F|=inf | **OK** (0.01s) |
| internal_combustion_engine_cpbar.eescode | Max iterations blk5 |F|=38.8 | Max iterations blk7 |F|=31452.4 | Max iterations blk7 |F|=0.3 | Max iterations blk7 |F|=0.3 | Max iterations blk7 |F|=0.3 | Other blk7 |F|=2078.6 | Other blk9 |F|=44089.1 | Max iterations blk7 |F|=2e+08 | Max iterations blk5 |F|=38.8 | Max iterations blk5 |F|=38.8 | Max iterations blk5 |F|=38.8 | Max iterations blk5 |F|=38.8 | Other blk7 |F|=24081.7 | Other blk7 |F|=inf | Other blk5 |F|=inf |
| lookup_demo.eescode | Other blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Other blk1 | Other blk1 | Other blk1 | Other blk1 | Other blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 |
| orc_co2.eescode | Max iterations blk28 |F|=1e+12 | Singular Jacobian blk28 |F|=491599.2 | Max iterations blk28 |F|=491599.2 | Max iterations blk28 |F|=491599.2 | Singular Jacobian blk28 |F|=491599.2 | Other blk28 | Other blk28 |F|=491599.2 | Max iterations blk28 |F|=1e+12 | Max iterations blk28 |F|=1e+12 | Max iterations blk28 |F|=1e+12 | Max iterations blk28 |F|=1e+12 | Max iterations blk28 |F|=1e+12 | Singular Jacobian blk28 |F|=491599.2 | Other blk28 | Other blk28 |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (0.45s) | Singular Jacobian blk21 |F|=inf | Max iterations blk21 |F|=inf | Max iterations blk21 |F|=inf | Singular Jacobian blk21 |F|=inf | Other blk21 | Other blk21 |F|=inf | **OK** (0.01s) | **OK** (0.43s) | **OK** (0.38s) | **OK** (0.38s) | **OK** (0.39s) | Singular Jacobian blk21 |F|=inf | Other blk21 |F|=inf | Other blk21 |F|=inf |
| orc_r245fa.eescode | **OK** (1.36s) | Singular Jacobian blk12 |F|=46791.8 | **OK** (0.76s) | Max iterations blk12 |F|=91775.1 | Max iterations blk6 |F|=5.1 | Other blk8 |F|=1930.2 | Other blk12 |F|=57429.4 | Other blk2 | **OK** (1.46s) | **OK** (1.36s) | **OK** (1.36s) | **OK** (0.64s) | Singular Jacobian blk12 |F|=64652.9 | Other blk8 |F|=inf | Other blk6 |F|=4120.2 |
| orc_simple.eescode | **OK** (0.27s) | **OK** (0.28s) | **OK** (0.28s) | **OK** (0.27s) | Max iterations blk6 |F|=2.8 | Other blk8 |F|=3417.9 | **OK** (1.38s) | Other blk2 | **OK** (0.28s) | **OK** (0.28s) | **OK** (0.27s) | **OK** (0.27s) | **OK** (0.27s) | Other blk8 | Other blk7 |F|=3e+26 |
| piston_compressor.eescode | **OK** (0.71s) | **OK** (0.04s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.04s) | Other blk4 |F|=0.2 | **OK** (0.01s) | Other blk2 | **OK** (0.74s) | **OK** (0.71s) | **OK** (0.74s) | **OK** (0.79s) | **OK** (0.00s) | Other blk4 |F|=7e+91 | **OK** (0.00s) |
| pressuredrop.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| rankine1.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| rankine2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.07s) | **OK** (0.07s) | **OK** (0.02s) | Other blk4 |F|=6e+07 | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.02s) | Other blk4 |F|=inf | **OK** (0.02s) |
| refrigeration1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration_compressor.eescode | **OK** (0.39s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.04s) | Other blk2 |F|=310.0 | **OK** (0.00s) | Other blk2 | **OK** (0.39s) | **OK** (0.36s) | **OK** (0.37s) | **OK** (0.46s) | **OK** (0.00s) | Other blk2 |F|=inf | **OK** (0.01s) |
| scroll_compressor.eescode | Max iterations blk34 |F|=1e+07 | **OK** (5.08s) | **OK** (5.98s) | Other blk34 |F|=7350.0 | Singular Jacobian blk34 |F|=1e+07 | Other blk34 | **OK** (1.00s) | Max iterations blk34 |F|=inf | Max iterations blk34 |F|=1e+07 | Max iterations blk34 |F|=1e+07 | Max iterations blk34 |F|=1e+07 | **OK** (2.60s) | Singular Jacobian blk34 |F|=5e+08 | Other blk34 |F|=inf | Other blk34 |F|=inf |
| simple_centrifugal_compressor.eescode | **OK** (0.00s) | Singular Jacobian blk1 |F|=0.3 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk1 |F|=0.3 | Other blk1 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk1 |F|=0.3 | Other blk1 |F|=6e+15 |
| storage_integraltable.eescode | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS |
| turbocompressor.eescode | **OK** (0.28s) | Singular Jacobian blk9 |F|=425994.2 | **OK** (0.12s) | Max iterations blk9 |F|=419375.3 | Singular Jacobian blk9 |F|=419155.1 | Other blk9 | **OK** (0.06s) | Max iterations blk9 |F|=6e+10 | **OK** (0.28s) | **OK** (0.31s) | **OK** (0.28s) | **OK** (0.28s) | Singular Jacobian blk9 |F|=424823.4 | Other blk9 |F|=inf | Other blk9 |F|=inf |
| turbocompressor_interpolate.eescode | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 |
| water_libr.eescode | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 |
| zorlu_heat_pump.eescode | Other blk63 |F|=6195.6 | Line search failed blk63 |F|=6196.8 | **OK** (15.04s) | Max iterations blk63 |F|=3e+06 | Max iterations blk63 |F|=65.5 | Other blk63 | **OK** (0.07s) | Max iterations blk63 |F|=2e+06 | Other blk63 |F|=6195.6 | Other blk63 |F|=6195.6 | Other blk63 |F|=6195.6 | Other blk63 |F|=6195.6 | Singular Jacobian blk63 |F|=6196.8 | Other blk63 |F|=inf | Other blk63 |F|=inf |

## Model Difficulty Ranking

Models ranked by number of configurations that failed to solve them.

| File | Failures / Configs | Failed Configurations |
|---|---:|---|
| water_libr.eescode | 30 / 30 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), TrustRegion + Hybrd K=5 (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials), Default + SymbolicReduction (with initials), Default + Tearing + SymbolicReduction (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), Default + MultiStart parallel (with initials), KINSOL linesearch (with initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| turbocompressor_interpolate.eescode | 30 / 30 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), TrustRegion + Hybrd K=5 (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials), Default + SymbolicReduction (with initials), Default + Tearing + SymbolicReduction (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), Default + MultiStart parallel (with initials), KINSOL linesearch (with initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| orc_co2.eescode | 19 / 30 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| cooling_tower.eescode | 19 / 30 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| internal_combustion_engine_cpbar.eescode | 18 / 30 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), KINSOL picard (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| cooling_tower2.eescode | 18 / 30 | Newton only (with initials), TrustRegion only (with initials), TrustRegion + Hybrd K=5 (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), KINSOL linesearch (with initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| heat_pump_MSTh_SB_R10.eescode | 18 / 30 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| zorlu_heat_pump.eescode | 17 / 30 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| air_screw_compressor_simple.eescode | 16 / 30 | Newton only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Newton only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), KINSOL linesearch (with initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| scroll_compressor.eescode | 16 / 30 | LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| lookup_demo.eescode | 15 / 30 | Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| orc_complex.eescode | 15 / 15 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), TrustRegion + Hybrd K=5 (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials), Default + SymbolicReduction (with initials), Default + Tearing + SymbolicReduction (with initials), Default + MultiStart parallel (with initials), KINSOL linesearch (with initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials) |
| air_screw_compressor.eescode | 14 / 30 | Newton only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Newton only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL linesearch (with initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| orc_extraction.eescode | 14 / 30 | LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| orc_r245fa.eescode | 13 / 30 | BisectionND only (with initials), Partitioned only (with initials), Newton only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| cooling_coil.eescode | 12 / 30 | BisectionND only (with initials), Default pipeline (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials), Default + SymbolicReduction (NO initials), Default + Tearing + SymbolicReduction (NO initials), Default + MultiStart parallel (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| condenser_3zones.eescode | 12 / 30 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| internal_combustion_engine.eescode | 10 / 30 | Partitioned only (with initials), Newton only (NO initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials) |
| simple_centrifugal_compressor.eescode | 10 / 30 | Newton only (with initials), Homotopy only (with initials), Partitioned only (with initials), Newton only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL fp/Anderson (with initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| exchangers3.eescode | 9 / 30 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials) |
| turbocompressor.eescode | 9 / 30 | BisectionND only (with initials), Newton only (NO initials), TrustRegion + Hybrd K=5 (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| humidair2.eescode | 9 / 30 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| cpbar.eescode | 8 / 30 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), TrustRegion + Hybrd K=5 (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials) |
| orc_simple.eescode | 8 / 30 | BisectionND only (with initials), Partitioned only (with initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials), KINSOL fp/Anderson (NO initials) |
| refrigeration_compressor.eescode | 6 / 30 | BisectionND only (with initials), Partitioned only (with initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials) |
| exchangers2.eescode | 5 / 30 | BisectionND only (with initials), BisectionND only (NO initials), KINSOL picard (with initials), KINSOL linesearch (NO initials), KINSOL picard (NO initials) |
| piston_compressor.eescode | 4 / 30 | Partitioned only (with initials), BisectionND only (NO initials), Partitioned only (NO initials), KINSOL picard (NO initials) |
| rankine2.eescode | 4 / 30 | BisectionND only (with initials), BisectionND only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials) |
| evaporator.eescode | 4 / 30 | Partitioned only (with initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials) |
| condenser_wet.eescode | 4 / 30 | Partitioned only (with initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials) |
| advanced_features.eescode | 4 / 30 | Partitioned only (with initials), Partitioned only (NO initials), KINSOL picard (with initials), KINSOL picard (NO initials) |

## Error Category Breakdown

Across all configurations and models:

| Error Category | Count | Fraction |
|---|---:|---:|
| Evaluation error | 39 | 10.0% |
| Line search failed | 1 | 0.3% |
| Max iterations | 86 | 22.1% |
| Other | 196 | 50.3% |
| Singular Jacobian | 38 | 9.7% |
| Unsupported function | 30 | 7.7% |

## Solver Pipeline Results

Shows which solver(s) were tried and their outcome for each model/config combination.

### Default pipeline (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| air_screw_compressor_simple.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | **OK** | Newton:Success |
| cooling_tower.eescode | **OK** | Newton:Success |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Newton:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Newton:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | Newton:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | Newton:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | Newton:Success |

### Newton only (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | FAIL | Newton:SingularJacobian |
| air_screw_compressor_simple.eescode | FAIL | Newton:SingularJacobian |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | **OK** | Newton:Success |
| cooling_tower.eescode | **OK** | Newton:Success |
| cooling_tower2.eescode | FAIL | Newton:SingularJacobian |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Newton:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Newton:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | Newton:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | Newton:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | FAIL | Newton:SingularJacobian |
| turbocompressor.eescode | **OK** | Newton:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | Newton:Success |

### TrustRegion only (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | TR:Success |
| air_screw_compressor.eescode | **OK** | TR:Success |
| air_screw_compressor_simple.eescode | **OK** | TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | TR:Success |
| condenser_wet.eescode | **OK** | TR:Success |
| cooling_coil.eescode | **OK** | TR:Success |
| cooling_tower.eescode | **OK** | TR:Success |
| cooling_tower2.eescode | FAIL | TR:MaxIterations |
| cpbar.eescode | **OK** | TR:Success |
| evaporator.eescode | **OK** | TR:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | TR:Success |
| exchangers3.eescode | **OK** | TR:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | TR:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | TR:Success |
| internal_combustion_engine.eescode | **OK** | TR:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | TR:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | TR:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | TR:Success |
| orc_r245fa.eescode | **OK** | TR:Success |
| orc_simple.eescode | **OK** | TR:Success |
| piston_compressor.eescode | **OK** | TR:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | TR:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | TR:Success |
| scroll_compressor.eescode | **OK** | TR:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | TR:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | TR:Success |

### TrustRegion + Hybrd K=5 (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | TR:Success |
| air_screw_compressor.eescode | **OK** | TR:Success |
| air_screw_compressor_simple.eescode | **OK** | TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | TR:Success |
| condenser_wet.eescode | **OK** | TR:Success |
| cooling_coil.eescode | **OK** | TR:Success |
| cooling_tower.eescode | **OK** | TR:Success |
| cooling_tower2.eescode | FAIL | TR:MaxIterations |
| cpbar.eescode | **OK** | TR:Success |
| evaporator.eescode | **OK** | TR:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | TR:Success |
| exchangers3.eescode | **OK** | TR:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | TR:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | TR:Success |
| internal_combustion_engine.eescode | **OK** | TR:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | TR:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | TR:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | TR:Success |
| orc_r245fa.eescode | **OK** | TR:Success |
| orc_simple.eescode | **OK** | TR:Success |
| piston_compressor.eescode | **OK** | TR:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | TR:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | TR:Success |
| scroll_compressor.eescode | **OK** | TR:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | TR:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | TR:Success |

### LevenbergMarquardt only (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | LM:Success |
| air_screw_compressor.eescode | FAIL | LM:SingularJacobian→LM:SingularJacobian |
| air_screw_compressor_simple.eescode | FAIL | LM:SingularJacobian→LM:SingularJacobian |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | LM:Success |
| condenser_wet.eescode | **OK** | LM:MaxIterations→LM:MaxIterations→LM:Success |
| cooling_coil.eescode | **OK** | LM:Success |
| cooling_tower.eescode | **OK** | LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:Success |
| cooling_tower2.eescode | FAIL | LM:SingularJacobian |
| cpbar.eescode | **OK** | LM:Success |
| evaporator.eescode | **OK** | LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | LM:Success |
| exchangers3.eescode | **OK** | LM:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | LM:MaxIterations→LM:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | LM:Success |
| internal_combustion_engine.eescode | **OK** | LM:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | LM:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | LM:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | FAIL | LM:MaxIterations→LM:MaxIterations→LM:MaxIterations |
| orc_r245fa.eescode | **OK** | LM:Success |
| orc_simple.eescode | **OK** | LM:Success |
| piston_compressor.eescode | **OK** | LM:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | LM:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | LM:Success |
| scroll_compressor.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | LM:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | LM:Success |

### BisectionND only (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | BisectionND:MaxIterations→BisectionND:Success |
| air_screw_compressor.eescode | FAIL | BisectionND:InvalidInput |
| air_screw_compressor_simple.eescode | FAIL | BisectionND:InvalidInput |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | BisectionND:InvalidInput |
| condenser_wet.eescode | **OK** | BisectionND:MaxIterations→BisectionND:Success |
| cooling_coil.eescode | FAIL | BisectionND:InvalidInput |
| cooling_tower.eescode | FAIL | BisectionND:InvalidInput |
| cooling_tower2.eescode | FAIL | BisectionND:InvalidInput |
| cpbar.eescode | FAIL | BisectionND:MaxIterations |
| evaporator.eescode | **OK** | BisectionND:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | FAIL | BisectionND:MaxIterations |
| exchangers3.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations |
| heat_pump_MSTh_SB_R10.eescode | FAIL | BisectionND:InvalidInput |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | BisectionND:MaxIterations |
| internal_combustion_engine.eescode | **OK** | BisectionND:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | FAIL | BisectionND:InvalidInput |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | FAIL | BisectionND:InvalidInput |
| orc_r245fa.eescode | FAIL | BisectionND:MaxIterations |
| orc_simple.eescode | FAIL | BisectionND:MaxIterations |
| piston_compressor.eescode | **OK** | BisectionND:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | FAIL | BisectionND:MaxIterations |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | FAIL | BisectionND:MaxIterations |
| scroll_compressor.eescode | FAIL | BisectionND:InvalidInput |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | FAIL | BisectionND:InvalidInput |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | BisectionND:InvalidInput |

### Homotopy only (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Homotopy:Success |
| air_screw_compressor.eescode | **OK** | Homotopy:Success |
| air_screw_compressor_simple.eescode | FAIL | Homotopy:MaxIterations |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Homotopy:Success |
| condenser_wet.eescode | **OK** | Homotopy:Success |
| cooling_coil.eescode | **OK** | Homotopy:Success |
| cooling_tower.eescode | **OK** | Homotopy:Success |
| cooling_tower2.eescode | FAIL | Homotopy:MaxIterations |
| cpbar.eescode | **OK** | Homotopy:Success |
| evaporator.eescode | **OK** | Homotopy:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Homotopy:Success |
| exchangers3.eescode | **OK** | Homotopy:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Homotopy:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Homotopy:Success |
| internal_combustion_engine.eescode | **OK** | Homotopy:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Homotopy:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | Homotopy:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | Homotopy:Success |
| orc_r245fa.eescode | **OK** | Homotopy:Success |
| orc_simple.eescode | **OK** | Homotopy:Success |
| piston_compressor.eescode | **OK** | Homotopy:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Homotopy:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Homotopy:Success |
| scroll_compressor.eescode | **OK** | Homotopy:Success |
| simple_centrifugal_compressor.eescode | FAIL | Homotopy:MaxIterations |
| turbocompressor.eescode | **OK** | Homotopy:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | Homotopy:Success |

### Partitioned only (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | FAIL | Part:MaxIterations |
| air_screw_compressor.eescode | FAIL | Part:MaxIterations |
| air_screw_compressor_simple.eescode | FAIL | Part:MaxIterations |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | Part:MaxIterations |
| condenser_wet.eescode | FAIL | Part:MaxIterations |
| cooling_coil.eescode | **OK** | Part:Success |
| cooling_tower.eescode | FAIL | Part:MaxIterations |
| cooling_tower2.eescode | **OK** | Part:Success |
| cpbar.eescode | FAIL | Part:MaxIterations |
| evaporator.eescode | FAIL | Part:MaxIterations |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Part:Success |
| exchangers3.eescode | FAIL | Part:MaxIterations |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Part:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | Part:MaxIterations |
| internal_combustion_engine.eescode | FAIL | Part:MaxIterations |
| internal_combustion_engine_cpbar.eescode | FAIL | Part:MaxIterations |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | FAIL | Part:MaxIterations |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | FAIL | Part:MaxIterations |
| orc_r245fa.eescode | FAIL | Part:MaxIterations |
| orc_simple.eescode | FAIL | Part:MaxIterations |
| piston_compressor.eescode | FAIL | Part:MaxIterations |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Part:MaxIterations→Part:MaxIterations→Part:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | FAIL | Part:MaxIterations |
| scroll_compressor.eescode | FAIL | Part:MaxIterations |
| simple_centrifugal_compressor.eescode | FAIL | Part:MaxIterations |
| turbocompressor.eescode | **OK** | Part:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Part:MaxIterations→Part:MaxIterations |

### Default + Tearing (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | — |
| air_screw_compressor_simple.eescode | **OK** | — |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | **OK** | — |
| cooling_tower.eescode | **OK** | — |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | — |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | — |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Newton:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | — |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Newton:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | — |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | — |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | — |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | — |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | — |

### Default + SymbolicReduction (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| air_screw_compressor_simple.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | **OK** | Newton:Success |
| cooling_tower.eescode | **OK** | Newton:Success |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Newton:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Newton:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | Newton:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | Newton:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | Newton:Success |

### Default + Tearing + SymbolicReduction (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | — |
| air_screw_compressor_simple.eescode | **OK** | — |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | **OK** | — |
| cooling_tower.eescode | **OK** | — |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | — |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | — |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Newton:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | — |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Newton:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | — |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | — |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | — |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | — |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | — |

### Default pipeline (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| air_screw_compressor_simple.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| cooling_tower.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:MaxIterations→Homotopy:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations |
| lookup_demo.eescode | FAIL | Newton:EvaluationError→TR:EvaluationError→LM:EvaluationError→BisectionND:EvaluationError→Homotopy:EvaluationError→Part:MaxIterations |
| orc_co2.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| orc_extraction.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |

### Newton only (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | FAIL | Newton:SingularJacobian |
| air_screw_compressor_simple.eescode | FAIL | Newton:SingularJacobian |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | **OK** | Newton:Success |
| cooling_tower.eescode | FAIL | Newton:SingularJacobian |
| cooling_tower2.eescode | FAIL | Newton:SingularJacobian |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Newton:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | FAIL | Newton:SingularJacobian |
| internal_combustion_engine_cpbar.eescode | FAIL | Newton:MaxIterations |
| lookup_demo.eescode | FAIL | Newton:EvaluationError |
| orc_co2.eescode | FAIL | Newton:SingularJacobian |
| orc_extraction.eescode | FAIL | Newton:SingularJacobian |
| orc_r245fa.eescode | FAIL | Newton:SingularJacobian |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | FAIL | Newton:SingularJacobian |
| turbocompressor.eescode | FAIL | Newton:SingularJacobian |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Newton:LineSearchFailed→Newton:LineSearchFailed |

### TrustRegion only (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | TR:Success |
| air_screw_compressor.eescode | **OK** | TR:Success |
| air_screw_compressor_simple.eescode | **OK** | TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | TR:MaxIterations→TR:LineSearchFailed |
| condenser_wet.eescode | **OK** | TR:Success |
| cooling_coil.eescode | **OK** | TR:Success |
| cooling_tower.eescode | FAIL | TR:MaxIterations |
| cooling_tower2.eescode | FAIL | TR:MaxIterations |
| cpbar.eescode | FAIL | TR:MaxIterations |
| evaporator.eescode | **OK** | TR:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | TR:Success |
| exchangers3.eescode | FAIL | TR:MaxIterations→TR:MaxIterations→TR:MaxIterations |
| heat_pump_MSTh_SB_R10.eescode | FAIL | TR:MaxIterations→TR:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | TR:MaxIterations→TR:MaxIterations |
| internal_combustion_engine.eescode | FAIL | TR:MaxIterations |
| internal_combustion_engine_cpbar.eescode | FAIL | TR:MaxIterations→TR:MaxIterations→TR:MaxIterations |
| lookup_demo.eescode | FAIL | TR:EvaluationError |
| orc_co2.eescode | FAIL | TR:MaxIterations |
| orc_extraction.eescode | FAIL | TR:MaxIterations |
| orc_r245fa.eescode | **OK** | TR:Success |
| orc_simple.eescode | **OK** | TR:Success |
| piston_compressor.eescode | **OK** | TR:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | TR:MaxIterations→TR:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | TR:Success |
| scroll_compressor.eescode | **OK** | TR:MaxIterations→TR:MaxIterations→TR:MaxIterations→TR:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | TR:MaxIterations→TR:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | TR:MaxIterations→TR:Success |

### TrustRegion + Hybrd K=5 (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | TR:Success |
| air_screw_compressor.eescode | **OK** | TR:Success |
| air_screw_compressor_simple.eescode | **OK** | TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | TR:MaxIterations→TR:MaxIterations→TR:MaxIterations |
| condenser_wet.eescode | **OK** | TR:Success |
| cooling_coil.eescode | FAIL | TR:LineSearchFailed→TR:MaxIterations |
| cooling_tower.eescode | FAIL | TR:MaxIterations |
| cooling_tower2.eescode | FAIL | TR:MaxIterations |
| cpbar.eescode | FAIL | TR:MaxIterations |
| evaporator.eescode | **OK** | TR:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | TR:Success |
| exchangers3.eescode | FAIL | TR:MaxIterations→TR:MaxIterations→TR:MaxIterations |
| heat_pump_MSTh_SB_R10.eescode | FAIL | TR:MaxIterations→TR:MaxIterations→TR:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | TR:Success |
| internal_combustion_engine.eescode | FAIL | TR:MaxIterations |
| internal_combustion_engine_cpbar.eescode | FAIL | TR:MaxIterations→TR:MaxIterations |
| lookup_demo.eescode | FAIL | TR:EvaluationError |
| orc_co2.eescode | FAIL | TR:MaxIterations |
| orc_extraction.eescode | FAIL | TR:MaxIterations |
| orc_r245fa.eescode | FAIL | TR:MaxIterations |
| orc_simple.eescode | **OK** | TR:Success |
| piston_compressor.eescode | **OK** | TR:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | TR:MaxIterations→TR:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | TR:MaxIterations→TR:Success |
| scroll_compressor.eescode | FAIL | TR:MaxIterations→TR:MaxIterations→TR:LineSearchFailed |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | FAIL | TR:MaxIterations |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | TR:MaxIterations |

### LevenbergMarquardt only (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | LM:Success |
| air_screw_compressor.eescode | FAIL | LM:SingularJacobian→LM:SingularJacobian |
| air_screw_compressor_simple.eescode | FAIL | LM:SingularJacobian→LM:SingularJacobian |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| condenser_wet.eescode | **OK** | LM:MaxIterations→LM:MaxIterations→LM:Success |
| cooling_coil.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| cooling_tower.eescode | FAIL | LM:SingularJacobian |
| cooling_tower2.eescode | FAIL | LM:SingularJacobian |
| cpbar.eescode | **OK** | LM:Success |
| evaporator.eescode | **OK** | LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | LM:Success |
| exchangers3.eescode | FAIL | LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations→LM:MaxIterations |
| heat_pump_MSTh_SB_R10.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| internal_combustion_engine.eescode | FAIL | LM:SingularJacobian |
| internal_combustion_engine_cpbar.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| lookup_demo.eescode | FAIL | LM:EvaluationError |
| orc_co2.eescode | FAIL | LM:SingularJacobian |
| orc_extraction.eescode | FAIL | LM:SingularJacobian |
| orc_r245fa.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| orc_simple.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |
| piston_compressor.eescode | **OK** | LM:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | LM:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | LM:Success |
| scroll_compressor.eescode | FAIL | LM:SingularJacobian→LM:SingularJacobian |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | FAIL | LM:SingularJacobian |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | LM:MaxIterations→LM:MaxIterations |

### BisectionND only (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | BisectionND:MaxIterations→BisectionND:Success |
| air_screw_compressor.eescode | FAIL | BisectionND:InvalidInput |
| air_screw_compressor_simple.eescode | FAIL | BisectionND:InvalidInput |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | BisectionND:InvalidInput |
| condenser_wet.eescode | **OK** | BisectionND:MaxIterations→BisectionND:Success |
| cooling_coil.eescode | FAIL | BisectionND:InvalidInput |
| cooling_tower.eescode | FAIL | BisectionND:InvalidInput |
| cooling_tower2.eescode | FAIL | BisectionND:InvalidInput |
| cpbar.eescode | FAIL | BisectionND:MaxIterations |
| evaporator.eescode | **OK** | BisectionND:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | FAIL | BisectionND:MaxIterations |
| exchangers3.eescode | FAIL | BisectionND:MaxIterations |
| heat_pump_MSTh_SB_R10.eescode | FAIL | BisectionND:InvalidInput |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations |
| internal_combustion_engine.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations |
| internal_combustion_engine_cpbar.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations |
| lookup_demo.eescode | FAIL | BisectionND:EvaluationError |
| orc_co2.eescode | FAIL | BisectionND:InvalidInput |
| orc_extraction.eescode | FAIL | BisectionND:InvalidInput |
| orc_r245fa.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations |
| orc_simple.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations |
| piston_compressor.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | FAIL | BisectionND:MaxIterations |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | FAIL | BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations→BisectionND:MaxIterations |
| scroll_compressor.eescode | FAIL | BisectionND:InvalidInput |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | FAIL | BisectionND:InvalidInput |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | BisectionND:InvalidInput |

### Homotopy only (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Homotopy:Success |
| air_screw_compressor.eescode | **OK** | Homotopy:Success |
| air_screw_compressor_simple.eescode | FAIL | Homotopy:MaxIterations |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Homotopy:Success |
| condenser_wet.eescode | **OK** | Homotopy:Success |
| cooling_coil.eescode | **OK** | Homotopy:Success |
| cooling_tower.eescode | FAIL | Homotopy:MaxIterations |
| cooling_tower2.eescode | FAIL | Homotopy:MaxIterations |
| cpbar.eescode | **OK** | Homotopy:Success |
| evaporator.eescode | **OK** | Homotopy:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Homotopy:Success |
| exchangers3.eescode | **OK** | Homotopy:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Homotopy:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Homotopy:Success |
| internal_combustion_engine.eescode | FAIL | Homotopy:MaxIterations |
| internal_combustion_engine_cpbar.eescode | FAIL | Homotopy:MaxIterations |
| lookup_demo.eescode | FAIL | Homotopy:EvaluationError |
| orc_co2.eescode | FAIL | Homotopy:MaxIterations |
| orc_extraction.eescode | FAIL | Homotopy:MaxIterations |
| orc_r245fa.eescode | FAIL | Homotopy:MaxIterations→Homotopy:MaxIterations |
| orc_simple.eescode | **OK** | Homotopy:Success |
| piston_compressor.eescode | **OK** | Homotopy:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Homotopy:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Homotopy:Success |
| scroll_compressor.eescode | **OK** | Homotopy:Success |
| simple_centrifugal_compressor.eescode | FAIL | Homotopy:MaxIterations |
| turbocompressor.eescode | **OK** | Homotopy:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | Homotopy:Success |

### Partitioned only (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | FAIL | Part:MaxIterations |
| air_screw_compressor.eescode | FAIL | Part:MaxIterations |
| air_screw_compressor_simple.eescode | FAIL | Part:MaxIterations |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | Part:MaxIterations |
| condenser_wet.eescode | FAIL | Part:MaxIterations |
| cooling_coil.eescode | FAIL | Part:MaxIterations |
| cooling_tower.eescode | FAIL | Part:MaxIterations |
| cooling_tower2.eescode | **OK** | Part:Success |
| cpbar.eescode | FAIL | Part:MaxIterations |
| evaporator.eescode | FAIL | Part:MaxIterations |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Part:Success |
| exchangers3.eescode | FAIL | Part:MaxIterations |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Part:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | Part:MaxIterations→Part:MaxIterations |
| internal_combustion_engine.eescode | FAIL | Part:MaxIterations |
| internal_combustion_engine_cpbar.eescode | FAIL | Part:MaxIterations |
| lookup_demo.eescode | FAIL | Part:MaxIterations |
| orc_co2.eescode | FAIL | Part:MaxIterations |
| orc_extraction.eescode | **OK** | Part:Success |
| orc_r245fa.eescode | FAIL | Part:MaxIterations |
| orc_simple.eescode | FAIL | Part:MaxIterations |
| piston_compressor.eescode | FAIL | Part:MaxIterations |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Part:MaxIterations→Part:MaxIterations→Part:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | FAIL | Part:MaxIterations |
| scroll_compressor.eescode | FAIL | Part:MaxIterations |
| simple_centrifugal_compressor.eescode | FAIL | Part:MaxIterations |
| turbocompressor.eescode | FAIL | Part:MaxIterations |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Part:MaxIterations→Part:MaxIterations |

### Default + Tearing (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | — |
| air_screw_compressor_simple.eescode | **OK** | — |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| cooling_tower.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | — |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations |
| lookup_demo.eescode | FAIL | Newton:EvaluationError→TR:EvaluationError→LM:EvaluationError→BisectionND:EvaluationError→Homotopy:EvaluationError→Part:MaxIterations |
| orc_co2.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| orc_extraction.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | — |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |

### Default + SymbolicReduction (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| air_screw_compressor_simple.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| cooling_tower.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:MaxIterations→Homotopy:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations |
| lookup_demo.eescode | FAIL | Newton:EvaluationError→TR:EvaluationError→LM:EvaluationError→BisectionND:EvaluationError→Homotopy:EvaluationError→Part:MaxIterations |
| orc_co2.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| orc_extraction.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |

### Default + Tearing + SymbolicReduction (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | — |
| air_screw_compressor_simple.eescode | **OK** | — |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| cooling_tower.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | — |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations |
| lookup_demo.eescode | FAIL | Newton:EvaluationError→TR:EvaluationError→LM:EvaluationError→BisectionND:EvaluationError→Homotopy:EvaluationError→Part:MaxIterations |
| orc_co2.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| orc_extraction.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | — |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |

### Default + MultiStart parallel (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| air_screw_compressor_simple.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |
| cooling_tower.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:MaxIterations→Homotopy:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations→Part:MaxIterations→Newton:MaxIterations→TR:MaxIterations→LM:MaxIterations→BisectionND:MaxIterations→Homotopy:MaxIterations |
| lookup_demo.eescode | FAIL | Newton:EvaluationError→TR:EvaluationError→LM:EvaluationError→BisectionND:EvaluationError→Homotopy:EvaluationError→Part:MaxIterations |
| orc_co2.eescode | FAIL | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:MaxIterations |
| orc_extraction.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Newton:LineSearchFailed→TR:MaxIterations→LM:MaxIterations→BisectionND:InvalidInput→Homotopy:MaxIterations |

### Default + MultiStart parallel (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Newton:Success |
| air_screw_compressor.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| air_screw_compressor_simple.eescode | **OK** | Newton:SingularJacobian→TR:Success |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Newton:Success |
| condenser_wet.eescode | **OK** | Newton:Success |
| cooling_coil.eescode | **OK** | Newton:Success |
| cooling_tower.eescode | **OK** | Newton:Success |
| cooling_tower2.eescode | **OK** | Newton:SingularJacobian→TR:MaxIterations→LM:SingularJacobian→BisectionND:InvalidInput→Homotopy:MaxIterations→Part:Success |
| cpbar.eescode | **OK** | Newton:Success |
| evaporator.eescode | **OK** | Newton:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Newton:Success |
| exchangers3.eescode | **OK** | Newton:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Newton:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Newton:Success |
| internal_combustion_engine.eescode | **OK** | Newton:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Newton:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | Newton:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | Newton:Success |
| orc_r245fa.eescode | **OK** | Newton:Success |
| orc_simple.eescode | **OK** | Newton:Success |
| piston_compressor.eescode | **OK** | Newton:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Newton:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Newton:Success |
| scroll_compressor.eescode | **OK** | Newton:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Newton:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | Newton:Success |

### KINSOL linesearch (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Kinsol:Success |
| air_screw_compressor.eescode | FAIL | Kinsol:LineSearchFailed→Kinsol:LineSearchFailed |
| air_screw_compressor_simple.eescode | FAIL | Kinsol:LineSearchFailed→Kinsol:LineSearchFailed |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | **OK** | Kinsol:Success |
| condenser_wet.eescode | **OK** | Kinsol:Success |
| cooling_coil.eescode | **OK** | Kinsol:Success |
| cooling_tower.eescode | **OK** | Kinsol:Success |
| cooling_tower2.eescode | FAIL | Kinsol:SingularJacobian |
| cpbar.eescode | **OK** | Kinsol:Success |
| evaporator.eescode | **OK** | Kinsol:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Kinsol:Success |
| exchangers3.eescode | **OK** | Kinsol:Success |
| heat_pump_MSTh_SB_R10.eescode | **OK** | Kinsol:Success |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Kinsol:Success |
| internal_combustion_engine.eescode | **OK** | Kinsol:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Kinsol:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | **OK** | Kinsol:Success |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | **OK** | Kinsol:Success |
| orc_r245fa.eescode | **OK** | Kinsol:Success |
| orc_simple.eescode | **OK** | Kinsol:Success |
| piston_compressor.eescode | **OK** | Kinsol:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Kinsol:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Kinsol:Success |
| scroll_compressor.eescode | **OK** | Kinsol:SingularJacobian→Kinsol:Success |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | **OK** | Kinsol:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | **OK** | Kinsol:Success |

### KINSOL picard (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | FAIL | Kinsol:MaxIterations |
| air_screw_compressor.eescode | FAIL | Kinsol:Diverged |
| air_screw_compressor_simple.eescode | FAIL | Kinsol:Diverged |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | Kinsol:Diverged |
| condenser_wet.eescode | FAIL | Kinsol:Diverged |
| cooling_coil.eescode | **OK** | Kinsol:Success |
| cooling_tower.eescode | FAIL | Kinsol:Diverged |
| cooling_tower2.eescode | FAIL | Kinsol:Diverged |
| cpbar.eescode | FAIL | Kinsol:Diverged |
| evaporator.eescode | FAIL | Kinsol:Diverged |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | FAIL | Kinsol:Diverged |
| exchangers3.eescode | FAIL | Kinsol:Diverged |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Kinsol:Diverged |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | Kinsol:Diverged |
| internal_combustion_engine.eescode | **OK** | Kinsol:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | Kinsol:Diverged |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | FAIL | Kinsol:Diverged |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | FAIL | Kinsol:Diverged |
| orc_r245fa.eescode | FAIL | Kinsol:Diverged |
| orc_simple.eescode | FAIL | Kinsol:Diverged |
| piston_compressor.eescode | **OK** | Kinsol:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | FAIL | Kinsol:Diverged |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | FAIL | Kinsol:Diverged |
| scroll_compressor.eescode | FAIL | Kinsol:Diverged |
| simple_centrifugal_compressor.eescode | FAIL | Kinsol:MaxIterations |
| turbocompressor.eescode | **OK** | Kinsol:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Kinsol:Diverged |

### KINSOL fp/Anderson (with initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Kinsol:Success |
| air_screw_compressor.eescode | FAIL | Kinsol:Diverged |
| air_screw_compressor_simple.eescode | FAIL | Kinsol:Diverged |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | Kinsol:Diverged |
| condenser_wet.eescode | **OK** | Kinsol:Success |
| cooling_coil.eescode | **OK** | Kinsol:Success |
| cooling_tower.eescode | FAIL | Kinsol:Diverged |
| cooling_tower2.eescode | FAIL | Kinsol:Diverged |
| cpbar.eescode | **OK** | Kinsol:Success |
| evaporator.eescode | **OK** | Kinsol:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Kinsol:Success |
| exchangers3.eescode | **OK** | Kinsol:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Kinsol:Diverged |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Kinsol:Success |
| internal_combustion_engine.eescode | **OK** | Kinsol:Success |
| internal_combustion_engine_cpbar.eescode | **OK** | Kinsol:Success |
| lookup_demo.eescode | **OK** | — |
| orc_co2.eescode | FAIL | Kinsol:Diverged |
| orc_complex.eescode | FAIL | — |
| orc_extraction.eescode | FAIL | Kinsol:Diverged |
| orc_r245fa.eescode | FAIL | Kinsol:MaxIterations |
| orc_simple.eescode | **OK** | Kinsol:Success |
| piston_compressor.eescode | **OK** | Kinsol:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Kinsol:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Kinsol:Success |
| scroll_compressor.eescode | FAIL | Kinsol:Diverged |
| simple_centrifugal_compressor.eescode | FAIL | Kinsol:MaxIterations |
| turbocompressor.eescode | **OK** | Kinsol:Success |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Kinsol:Diverged |

### KINSOL linesearch (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Kinsol:Success |
| air_screw_compressor.eescode | FAIL | Kinsol:LineSearchFailed→Kinsol:LineSearchFailed |
| air_screw_compressor_simple.eescode | FAIL | Kinsol:LineSearchFailed→Kinsol:LineSearchFailed |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | Kinsol:SingularJacobian |
| condenser_wet.eescode | **OK** | Kinsol:Success |
| cooling_coil.eescode | **OK** | Kinsol:Success |
| cooling_tower.eescode | FAIL | Kinsol:SingularJacobian |
| cooling_tower2.eescode | FAIL | Kinsol:SingularJacobian |
| cpbar.eescode | **OK** | Kinsol:Success |
| evaporator.eescode | **OK** | Kinsol:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | FAIL | Kinsol:SingularJacobian→Kinsol:SingularJacobian→Kinsol:SingularJacobian→Kinsol:SingularJacobian |
| exchangers3.eescode | **OK** | Kinsol:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Kinsol:MaxIterations |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | **OK** | Kinsol:Success |
| internal_combustion_engine.eescode | FAIL | Kinsol:SingularJacobian |
| internal_combustion_engine_cpbar.eescode | FAIL | Kinsol:MaxIterations→Kinsol:MaxIterations |
| lookup_demo.eescode | FAIL | Kinsol:EvaluationError |
| orc_co2.eescode | FAIL | Kinsol:SingularJacobian |
| orc_extraction.eescode | FAIL | Kinsol:SingularJacobian |
| orc_r245fa.eescode | FAIL | Kinsol:SingularJacobian→Kinsol:SingularJacobian |
| orc_simple.eescode | **OK** | Kinsol:Success |
| piston_compressor.eescode | **OK** | Kinsol:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Kinsol:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Kinsol:Success |
| scroll_compressor.eescode | FAIL | Kinsol:SingularJacobian |
| simple_centrifugal_compressor.eescode | **OK** | — |
| turbocompressor.eescode | FAIL | Kinsol:SingularJacobian |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Kinsol:SingularJacobian→Kinsol:SingularJacobian |

### KINSOL picard (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | FAIL | Kinsol:MaxIterations |
| air_screw_compressor.eescode | FAIL | Kinsol:Diverged |
| air_screw_compressor_simple.eescode | FAIL | Kinsol:Diverged |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | Kinsol:Diverged |
| condenser_wet.eescode | FAIL | Kinsol:Diverged |
| cooling_coil.eescode | FAIL | Kinsol:Diverged |
| cooling_tower.eescode | FAIL | Kinsol:Diverged |
| cooling_tower2.eescode | FAIL | Kinsol:Diverged |
| cpbar.eescode | FAIL | Kinsol:MaxIterations |
| evaporator.eescode | FAIL | Kinsol:Diverged |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | FAIL | Kinsol:Diverged |
| exchangers3.eescode | FAIL | Kinsol:Diverged |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Kinsol:Diverged |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | Kinsol:Diverged |
| internal_combustion_engine.eescode | FAIL | Kinsol:Diverged |
| internal_combustion_engine_cpbar.eescode | FAIL | Kinsol:Diverged |
| lookup_demo.eescode | FAIL | Kinsol:EvaluationError |
| orc_co2.eescode | FAIL | Kinsol:Diverged |
| orc_extraction.eescode | FAIL | Kinsol:Diverged |
| orc_r245fa.eescode | FAIL | Kinsol:Diverged |
| orc_simple.eescode | FAIL | Kinsol:Diverged |
| piston_compressor.eescode | FAIL | Kinsol:MaxIterations |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | FAIL | Kinsol:Diverged |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | FAIL | Kinsol:Diverged |
| scroll_compressor.eescode | FAIL | Kinsol:Diverged |
| simple_centrifugal_compressor.eescode | FAIL | Kinsol:MaxIterations |
| turbocompressor.eescode | FAIL | Kinsol:Diverged |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Kinsol:Diverged |

### KINSOL fp/Anderson (NO initials)

| File | Status | Pipeline (solver:result) |
|---|:---:|---|
| advanced_features.eescode | **OK** | Kinsol:Success |
| air_screw_compressor.eescode | FAIL | Kinsol:Diverged |
| air_screw_compressor_simple.eescode | FAIL | Kinsol:Diverged |
| boiler_cpbar.eescode | **OK** | — |
| boiler_cpbar2.eescode | **OK** | — |
| compressor_refrigeration_simple.eescode | **OK** | — |
| condenser_3zones.eescode | FAIL | Kinsol:Diverged |
| condenser_wet.eescode | **OK** | Kinsol:Success |
| cooling_coil.eescode | FAIL | Kinsol:MaxIterations |
| cooling_tower.eescode | FAIL | Kinsol:Diverged |
| cooling_tower2.eescode | FAIL | Kinsol:Diverged |
| cpbar.eescode | **OK** | Kinsol:Success |
| evaporator.eescode | **OK** | Kinsol:Success |
| exchangers1.eescode | **OK** | — |
| exchangers2.eescode | **OK** | Kinsol:Success |
| exchangers3.eescode | **OK** | Kinsol:Success |
| heat_pump_MSTh_SB_R10.eescode | FAIL | Kinsol:Diverged |
| humidair1.eescode | **OK** | — |
| humidair2.eescode | FAIL | Kinsol:MaxIterations→Kinsol:Diverged |
| internal_combustion_engine.eescode | **OK** | Kinsol:Success |
| internal_combustion_engine_cpbar.eescode | FAIL | Kinsol:Diverged |
| lookup_demo.eescode | FAIL | Kinsol:EvaluationError |
| orc_co2.eescode | FAIL | Kinsol:Diverged |
| orc_extraction.eescode | FAIL | Kinsol:Diverged |
| orc_r245fa.eescode | FAIL | Kinsol:MaxIterations |
| orc_simple.eescode | FAIL | Kinsol:MaxIterations |
| piston_compressor.eescode | **OK** | Kinsol:Success |
| pressuredrop.eescode | **OK** | — |
| rankine1.eescode | **OK** | — |
| rankine2.eescode | **OK** | Kinsol:Success |
| refrigeration1.eescode | **OK** | — |
| refrigeration2.eescode | **OK** | — |
| refrigeration3.eescode | **OK** | — |
| refrigeration_compressor.eescode | **OK** | Kinsol:Success |
| scroll_compressor.eescode | FAIL | Kinsol:Diverged |
| simple_centrifugal_compressor.eescode | FAIL | Kinsol:MaxIterations |
| turbocompressor.eescode | FAIL | Kinsol:Diverged |
| turbocompressor_interpolate.eescode | FAIL | — |
| water_libr.eescode | FAIL | — |
| zorlu_heat_pump.eescode | FAIL | Kinsol:Diverged |


## Detailed Error Messages

### Default pipeline (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Newton only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Singular Jacobian | 13 | 322614.05 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | 304355.91 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| cooling_tower2.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| simple_centrifugal_compressor.eescode | Singular Jacobian | 1 | 0.30 | Block 11 (size 1, vars: A) failed: SingularJacobian -  Initial //F//_inf = 0.295303, best achieved = |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### TrustRegion only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_tower2.eescode | Max iterations | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### TrustRegion + Hybrd K=5 (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_tower2.eescode | Max iterations | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### LevenbergMarquardt only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Singular Jacobian | 13 | 365861.14 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | 365861.14 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| cooling_tower2.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| orc_extraction.eescode | Max iterations | 21 | 0.00 | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Levenbe |
| scroll_compressor.eescode | Max iterations | 34 | 0.81 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [LevenbergMarquardt]  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### BisectionND only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: InvalidInput - [BisectionND] Bis |
| air_screw_compressor_simple.eescode | Other | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: InvalidInput - [BisectionND] Bis |
| condenser_3zones.eescode | Other | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: InvalidInput - [ |
| cooling_coil.eescode | Other | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: InvalidInput - [BisectionND]  |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: InvalidInput - [BisectionND] Bisec |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: InvalidInput - [BisectionND] Bisec |
| cpbar.eescode | Other | 5 | 4.5e+08 | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [BisectionND] BisectionND: max it |
| exchangers2.eescode | Other | 4 | 355.59 | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: MaxIterations - [BisectionND] BisectionND: max  |
| exchangers3.eescode | Other | 3 | 123803.12 | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [BisectionND] BisectionND: max |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: InvalidInput - [BisectionND] Bise |
| humidair2.eescode | Other | 5 | 1.97 | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| internal_combustion_engine_cpbar.eescode | Other | 2 | 0.00 | Block 39 (size 2, vars: t_7, c_p_g_67) failed: MaxIterations - [BisectionND] BisectionND: max iterat |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: InvalidInput - [BisectionND] BisectionND: bl |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: InvalidInput - [Bisectio |
| orc_r245fa.eescode | Other | 8 | 904.08 | Block 51 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| orc_simple.eescode | Other | 7 | 0.00 | Block 75 (size 7, vars: h_hf_su_tp, T_hf_su_tp, h_hf_ex_tp, ...) failed: MaxIterations - [BisectionN |
| rankine2.eescode | Other | 4 | 42594.19 | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [BisectionND |
| refrigeration_compressor.eescode | Other | 4 | 0.10 | Block 43 (size 4, vars: C, epsilon_v_2, V_dot_s, ...) failed: MaxIterations - [BisectionND] Bisectio |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: InvalidInput - [BisectionND] Bisectio |
| turbocompressor.eescode | Other | 9 | — | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: InvalidInput - [BisectionND] BisectionND: block |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | — | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: InvalidInput - [BisectionND]  |

### Homotopy only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor_simple.eescode | Other | 13 | 419347.47 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Homotopy] Homot |
| cooling_tower2.eescode | Other | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [Homotopy] Homotop |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| simple_centrifugal_compressor.eescode | Other | 1 | 0.30 | Block 11 (size 1, vars: A) failed: MaxIterations - [Homotopy] Homotopy: did not converge. Last t=1,  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Partitioned only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| advanced_features.eescode | Other | 2 | — | Block 43 (size 2, vars: st_total, st_diff) failed: MaxIterations -  Initial //F//_inf = 10, best ach |
| air_screw_compressor.eescode | Max iterations | 13 | 2.1e+08 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| air_screw_compressor_simple.eescode | Max iterations | 13 | 2.1e+08 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| condenser_3zones.eescode | Max iterations | 62 | 3335.53 | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| condenser_wet.eescode | Other | 2 | — | Block 19 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 178978, best |
| cooling_tower.eescode | Max iterations | 11 | 27486.88 | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [Partitioned] Pa |
| cpbar.eescode | Max iterations | 5 | 4.5e+08 | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [Partitioned] Partitioned solver: |
| evaporator.eescode | Other | 2 | — | Block 23 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 13869.6, bes |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations -  Initial //F//_inf = 310600, b |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 1.8e+158 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Partitioned] Par |
| humidair2.eescode | Max iterations | 5 | 9.8e+08 | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [Partitioned] Partitioned so |
| internal_combustion_engine.eescode | Other | 2 | — | Block 16 (size 2, vars: h_ex_2, h_su_2) failed: MaxIterations -  Initial //F//_inf = 8.9989e-08, bes |
| internal_combustion_engine_cpbar.eescode | Other | 2 | — | Block 38 (size 2, vars: t_6, c_p_g_6) failed: MaxIterations -  Initial //F//_inf = 8.24002e-10, best |
| orc_co2.eescode | Max iterations | 28 | 1.4e+14 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Partitioned] Partitioned so |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| orc_extraction.eescode | Max iterations | 21 | 5.3e+47 | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Partiti |
| orc_r245fa.eescode | Other | 2 | — | Block 55 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 1377.92, best a |
| orc_simple.eescode | Other | 2 | — | Block 53 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 1.19908e-08, be |
| piston_compressor.eescode | Other | 2 | — | Block 27 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations -  Initial //F//_inf = 3.31784e-0 |
| refrigeration_compressor.eescode | Other | 2 | — | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations -  Initial //F//_inf = 1.31242e-0 |
| scroll_compressor.eescode | Max iterations | 34 | inf | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Partitioned] Partiti |
| simple_centrifugal_compressor.eescode | Other | 1 | — | Block 11 (size 1, vars: A) failed: MaxIterations -  Initial //F//_inf = 0.295303, best achieved = 0. |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Max iterations | 63 | 4576.05 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations - [Partitioned] |

### Default + Tearing (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Default + SymbolicReduction (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Default + Tearing + SymbolicReduction (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Default pipeline (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_coil.eescode | Other | 12 | 468.47 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations -  Initial //F/ |
| cooling_tower.eescode | Max iterations | 11 | 7.7e+18 | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [TrustRegion] Tr |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 4.6e+06 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Newton] Max iter |
| internal_combustion_engine_cpbar.eescode | Max iterations | 5 | 38.80 | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: MaxIterations - [Newton] Max iterations (100) re |
| lookup_demo.eescode | Other | 1 | — | Block 1 (size 1, vars: h_interp) failed: MaxIterations - [Newton] INTERPOLATE(): lookup table 'data' |
| orc_co2.eescode | Max iterations | 28 | 1.0e+12 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| scroll_compressor.eescode | Max iterations | 34 | 1.3e+07 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [TrustRegion] Trust r |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | 6195.62 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations -  Initial //F/ |

### Newton only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Singular Jacobian | 13 | 322614.05 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | 304355.91 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| cooling_tower.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| cooling_tower2.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 808738.89 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Newton] Max iter |
| internal_combustion_engine.eescode | Singular Jacobian | 5 | 399999.00 | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: SingularJacobian -  Initial //F//_i |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | 31452.39 | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [Newton] Max iterations (100 |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [Newton] INTERPOLATE(): lookup table 'dat |
| orc_co2.eescode | Singular Jacobian | 28 | 491599.19 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: SingularJacobian -  Initial //F//_inf = 4915 |
| orc_extraction.eescode | Singular Jacobian | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: SingularJacobian |
| orc_r245fa.eescode | Singular Jacobian | 12 | 46791.80 | Block 82 (size 12, vars: rho_hf_su_ev, nu_hf_su_ev, k_hf_su_ev, ...) failed: SingularJacobian -  Ini |
| simple_centrifugal_compressor.eescode | Singular Jacobian | 1 | 0.30 | Block 11 (size 1, vars: A) failed: SingularJacobian -  Initial //F//_inf = 0.295303, best achieved = |
| turbocompressor.eescode | Singular Jacobian | 9 | 425994.21 | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: SingularJacobian -  Initial //F//_inf = 419405, |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Line search failed | 63 | 6196.76 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: LineSearchFailed - [Newton] L |

### TrustRegion only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Other | 62 | 19851.98 | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: LineSearchFailed |
| cooling_tower.eescode | Max iterations | 11 | 2.2e+07 | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| cooling_tower2.eescode | Max iterations | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| cpbar.eescode | Max iterations | 5 | 1.1e+08 | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [TrustRegion] Trust region: Max i |
| exchangers3.eescode | Max iterations | 3 | 10003.56 | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [TrustRegion] Trust region: Ma |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 68903.97 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [TrustRegion] Tru |
| humidair2.eescode | Max iterations | 5 | 446.65 | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| internal_combustion_engine.eescode | Max iterations | 5 | 399999.00 | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: MaxIterations - [TrustRegion] Trust |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | 0.33 | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [TrustRegion] INTERPOLATE(): lookup table |
| orc_co2.eescode | Max iterations | 28 | 491599.19 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| orc_extraction.eescode | Max iterations | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [TrustRe |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### TrustRegion + Hybrd K=5 (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Max iterations | 62 | 19287.98 | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| cooling_coil.eescode | Max iterations | 12 | 1.61 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [TrustRegion] |
| cooling_tower.eescode | Max iterations | 11 | 2.2e+07 | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| cooling_tower2.eescode | Max iterations | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| cpbar.eescode | Max iterations | 5 | 1.1e+08 | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [TrustRegion] Trust region: Max i |
| exchangers3.eescode | Max iterations | 3 | 5643.17 | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [TrustRegion] Trust region: Ma |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 47705.29 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [TrustRegion] Tru |
| internal_combustion_engine.eescode | Max iterations | 5 | 399999.00 | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: MaxIterations - [TrustRegion] Trust |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | 0.33 | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [TrustRegion] INTERPOLATE(): lookup table |
| orc_co2.eescode | Max iterations | 28 | 491599.19 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| orc_extraction.eescode | Max iterations | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [TrustRe |
| orc_r245fa.eescode | Max iterations | 12 | 91775.14 | Block 82 (size 12, vars: rho_hf_su_ev, nu_hf_su_ev, k_hf_su_ev, ...) failed: MaxIterations - [TrustR |
| scroll_compressor.eescode | Other | 34 | 7349.97 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: LineSearchFailed -  Initial //F//_inf |
| turbocompressor.eescode | Max iterations | 9 | 419375.30 | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: MaxIterations - [TrustRegion] Trust region: Max |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Max iterations | 63 | 2.7e+06 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations - [TrustRegion] |

### LevenbergMarquardt only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Singular Jacobian | 13 | 365861.14 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | 365861.14 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| condenser_3zones.eescode | Max iterations | 62 | 921.23 | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| cooling_coil.eescode | Max iterations | 12 | 57.96 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [LevenbergMar |
| cooling_tower.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| cooling_tower2.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| exchangers3.eescode | Max iterations | 3 | 8.13 | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [LevenbergMarquardt] Levenberg |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 1674.66 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [LevenbergMarquar |
| humidair2.eescode | Max iterations | 5 | 0.01 | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [LevenbergMarquardt] Levenbe |
| internal_combustion_engine.eescode | Singular Jacobian | 5 | 399999.00 | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: SingularJacobian -  Initial //F//_i |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | 0.33 | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [LevenbergMarquardt] Levenbe |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [LevenbergMarquardt] INTERPOLATE(): looku |
| orc_co2.eescode | Singular Jacobian | 28 | 491599.19 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: SingularJacobian -  Initial //F//_inf = 4915 |
| orc_extraction.eescode | Singular Jacobian | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: SingularJacobian |
| orc_r245fa.eescode | Max iterations | 6 | 5.11 | Block 59 (size 6, vars: h_cf_ex_tp, T_cf_ex_tp, h_cf_ex_cd, ...) failed: MaxIterations - [LevenbergM |
| orc_simple.eescode | Max iterations | 6 | 2.85 | Block 57 (size 6, vars: h_cf_ex_tp, T_cf_ex_tp, h_cf_ex_cd, ...) failed: MaxIterations - [LevenbergM |
| scroll_compressor.eescode | Singular Jacobian | 34 | 1.2e+07 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: SingularJacobian -  Initial //F//_inf |
| turbocompressor.eescode | Singular Jacobian | 9 | 419155.15 | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: SingularJacobian -  Initial //F//_inf = 419405, |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Max iterations | 63 | 65.47 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations - [LevenbergMar |

### BisectionND only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: InvalidInput - [BisectionND] Bis |
| air_screw_compressor_simple.eescode | Other | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: InvalidInput - [BisectionND] Bis |
| condenser_3zones.eescode | Other | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: InvalidInput - [ |
| cooling_coil.eescode | Other | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: InvalidInput - [BisectionND]  |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: InvalidInput - [BisectionND] Bisec |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: InvalidInput - [BisectionND] Bisec |
| cpbar.eescode | Other | 5 | 1.1e+08 | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [BisectionND] BisectionND: max it |
| exchangers2.eescode | Other | 4 | 4.4e+12 | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: MaxIterations - [BisectionND] BisectionND: max  |
| exchangers3.eescode | Other | 3 | 4.3e+09 | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [BisectionND] BisectionND: max |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: InvalidInput - [BisectionND] Bise |
| humidair2.eescode | Other | 5 | 9930.00 | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| internal_combustion_engine.eescode | Other | 5 | 4198.30 | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: MaxIterations - [BisectionND] Bisec |
| internal_combustion_engine_cpbar.eescode | Other | 7 | 2078.59 | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [BisectionND] BisectionND: no evaluable p |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: InvalidInput - [BisectionND] BisectionND: bl |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: InvalidInput - [Bisectio |
| orc_r245fa.eescode | Other | 8 | 1930.21 | Block 51 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| orc_simple.eescode | Other | 8 | 3417.88 | Block 49 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| piston_compressor.eescode | Other | 4 | 0.17 | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: MaxIterations - [BisectionND] Bise |
| rankine2.eescode | Other | 4 | 6.3e+07 | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [BisectionND |
| refrigeration_compressor.eescode | Other | 2 | 310.04 | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations - [BisectionND] BisectionND: max  |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: InvalidInput - [BisectionND] Bisectio |
| turbocompressor.eescode | Other | 9 | — | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: InvalidInput - [BisectionND] BisectionND: block |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | — | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: InvalidInput - [BisectionND]  |

### Homotopy only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor_simple.eescode | Other | 13 | 419347.47 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Homotopy] Homot |
| cooling_tower.eescode | Other | 11 | 2.2e+07 | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [Homotopy] Homotop |
| cooling_tower2.eescode | Other | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [Homotopy] Homotop |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | 3.4e+06 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Homotopy] Homoto |
| internal_combustion_engine.eescode | Other | 5 | 399999.00 | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: MaxIterations - [Homotopy] Homotopy |
| internal_combustion_engine_cpbar.eescode | Other | 9 | 44089.11 | Block 42 (size 9, vars: t_8, c_p_g_78, Q_dot_gw, ...) failed: MaxIterations - [Homotopy] Homotopy: d |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [Homotopy] INTERPOLATE(): lookup table 'd |
| orc_co2.eescode | Other | 28 | 491599.19 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Homotopy] Homotopy: did not |
| orc_extraction.eescode | Other | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Homotop |
| orc_r245fa.eescode | Other | 12 | 57429.39 | Block 82 (size 12, vars: rho_hf_su_ev, nu_hf_su_ev, k_hf_su_ev, ...) failed: MaxIterations - [Homoto |
| simple_centrifugal_compressor.eescode | Other | 1 | 0.30 | Block 11 (size 1, vars: A) failed: MaxIterations - [Homotopy] Homotopy: did not converge. Last t=1,  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Partitioned only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| advanced_features.eescode | Other | 2 | — | Block 43 (size 2, vars: st_total, st_diff) failed: MaxIterations -  Initial //F//_inf = 10, best ach |
| air_screw_compressor.eescode | Max iterations | 13 | 2.1e+08 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| air_screw_compressor_simple.eescode | Max iterations | 13 | 2.1e+08 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| condenser_3zones.eescode | Max iterations | 62 | 1.3e+07 | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| condenser_wet.eescode | Other | 2 | — | Block 19 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 178978, best |
| cooling_coil.eescode | Max iterations | 12 | 3.6e+08 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [Partitioned] |
| cooling_tower.eescode | Max iterations | 11 | 7.7e+18 | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [Partitioned] Pa |
| cpbar.eescode | Max iterations | 5 | 1.1e+08 | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [Partitioned] Partitioned solver: |
| evaporator.eescode | Other | 2 | — | Block 23 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 13869.6, bes |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations -  Initial //F//_inf = 4.38491e+ |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 2.5e+82 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Partitioned] Par |
| humidair2.eescode | Max iterations | 5 | 136222.20 | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [Partitioned] Partitioned so |
| internal_combustion_engine.eescode | Other | 2 | — | Block 16 (size 2, vars: h_ex_2, h_su_2) failed: MaxIterations -  Initial //F//_inf = 190817, best ac |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | 1.9e+08 | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [Partitioned] Partitioned so |
| lookup_demo.eescode | Other | 1 | — | Block 1 (size 1, vars: h_interp) failed: MaxIterations |
| orc_co2.eescode | Max iterations | 28 | 1.0e+12 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Partitioned] Partitioned so |
| orc_r245fa.eescode | Other | 2 | — | Block 55 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 76754.1, best a |
| orc_simple.eescode | Other | 2 | — | Block 53 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 76754.1, best a |
| piston_compressor.eescode | Other | 2 | — | Block 27 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations -  Initial //F//_inf = 11041.2, b |
| refrigeration_compressor.eescode | Other | 2 | — | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations -  Initial //F//_inf = 9483.31, b |
| scroll_compressor.eescode | Max iterations | 34 | inf | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Partitioned] Partiti |
| simple_centrifugal_compressor.eescode | Other | 1 | — | Block 11 (size 1, vars: A) failed: MaxIterations -  Initial //F//_inf = 0.295303, best achieved = 0. |
| turbocompressor.eescode | Max iterations | 9 | 6.3e+10 | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: MaxIterations - [Partitioned] Partitioned solve |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Max iterations | 63 | 1.7e+06 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations - [Partitioned] |

### Default + Tearing (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_coil.eescode | Other | 12 | 468.46 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - Tearing: sing |
| cooling_tower.eescode | Max iterations | 11 | 7.7e+18 | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - Tearing: singula |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 4.6e+06 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - Tearing: singular |
| internal_combustion_engine_cpbar.eescode | Max iterations | 5 | 38.80 | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: MaxIterations - Tearing: singular Schur compleme |
| lookup_demo.eescode | Other | 1 | — | Block 1 (size 1, vars: h_interp) failed: MaxIterations - [Newton] INTERPOLATE(): lookup table 'data' |
| orc_co2.eescode | Max iterations | 28 | 1.0e+12 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - Tearing: singular Schur comp |
| scroll_compressor.eescode | Max iterations | 34 | 1.3e+07 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - Tearing: singular Sch |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | 6195.62 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations - Tearing: sing |

### Default + SymbolicReduction (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_coil.eescode | Other | 12 | 468.46 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations -  Initial //F/ |
| cooling_tower.eescode | Max iterations | 11 | 7.7e+18 | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [TrustRegion] Tr |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 4.6e+06 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Newton] Max iter |
| internal_combustion_engine_cpbar.eescode | Max iterations | 5 | 38.80 | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: MaxIterations - [Newton] Max iterations (100) re |
| lookup_demo.eescode | Other | 1 | — | Block 1 (size 1, vars: h_interp) failed: MaxIterations - [Newton] INTERPOLATE(): lookup table 'data' |
| orc_co2.eescode | Max iterations | 28 | 1.0e+12 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| scroll_compressor.eescode | Max iterations | 34 | 1.3e+07 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [TrustRegion] Trust r |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | 6195.62 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations -  Initial //F/ |

### Default + Tearing + SymbolicReduction (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_coil.eescode | Other | 12 | 468.46 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - Tearing: sing |
| cooling_tower.eescode | Max iterations | 11 | 7.7e+18 | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - Tearing: singula |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | 4.6e+06 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - Tearing: singular |
| internal_combustion_engine_cpbar.eescode | Max iterations | 5 | 38.80 | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: MaxIterations - Tearing: singular Schur compleme |
| lookup_demo.eescode | Other | 1 | — | Block 1 (size 1, vars: h_interp) failed: MaxIterations - [Newton] INTERPOLATE(): lookup table 'data' |
| orc_co2.eescode | Max iterations | 28 | 1.0e+12 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - Tearing: singular Schur comp |
| scroll_compressor.eescode | Max iterations | 34 | 1.3e+07 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - Tearing: singular Sch |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | 6195.62 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations - Tearing: sing |

### Default + MultiStart parallel (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_coil.eescode | Other | 12 | 468.46 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations -  Initial //F/ |
| cooling_tower.eescode | Max iterations | 11 | 7.7e+18 | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [TrustRegion] Tr |
| internal_combustion_engine_cpbar.eescode | Max iterations | 5 | 38.80 | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: MaxIterations - [Newton] Max iterations (100) re |
| lookup_demo.eescode | Other | 1 | — | Block 1 (size 1, vars: h_interp) failed: MaxIterations - [Newton] INTERPOLATE(): lookup table 'data' |
| orc_co2.eescode | Max iterations | 28 | 1.0e+12 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | 6195.62 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: MaxIterations -  Initial //F/ |

### Default + MultiStart parallel (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### KINSOL linesearch (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | 177866.45 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: LineSearchFailed - [Kinsol] KINS |
| air_screw_compressor_simple.eescode | Other | 13 | 22785.53 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: LineSearchFailed - [Kinsol] KINS |
| cooling_tower2.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian - [Kinsol] KINSOL |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### KINSOL picard (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| advanced_features.eescode | Other | 2 | 3.1e+90 | Block 43 (size 2, vars: st_total, st_diff) failed: MaxIterations - [Kinsol] KINSOL(Picard): max iter |
| air_screw_compressor.eescode | Other | 13 | inf | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(Picar |
| air_screw_compressor_simple.eescode | Other | 13 | inf | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(Picar |
| condenser_3zones.eescode | Other | 62 | inf | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: Diverged - [Kins |
| condenser_wet.eescode | Other | 2 | inf | Block 19 (size 2, vars: M_dot_a, M_dot_cd) failed: Diverged - [Kinsol] KINSOL(Picard): residual beca |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(Picard) |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(Picard) |
| cpbar.eescode | Other | 5 | inf | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual beca |
| evaporator.eescode | Other | 2 | inf | Block 23 (size 2, vars: M_dot_a, M_dot_cd) failed: Diverged - [Kinsol] KINSOL(Picard): residual beca |
| exchangers2.eescode | Other | 4 | inf | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual be |
| exchangers3.eescode | Other | 3 | inf | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: Diverged - [Kinsol] KINSOL(Picard): residual b |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | inf | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: Diverged - [Kinsol] KINSOL(Picard |
| humidair2.eescode | Other | 5 | inf | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual |
| internal_combustion_engine_cpbar.eescode | Other | 2 | inf | Block 39 (size 2, vars: t_7, c_p_g_67) failed: Diverged - [Kinsol] KINSOL(Picard): residual became n |
| orc_co2.eescode | Other | 28 | inf | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| orc_extraction.eescode | Other | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: Diverged - [Kinsol] KINS |
| orc_r245fa.eescode | Other | 8 | inf | Block 51 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: Diverged - [Kinsol] KI |
| orc_simple.eescode | Other | 7 | — | Block 75 (size 7, vars: h_hf_su_tp, T_hf_su_tp, h_hf_ex_tp, ...) failed: Diverged - [Kinsol] KINSOL( |
| rankine2.eescode | Other | 4 | inf | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: Diverged - [Kinsol] KINSOL(P |
| refrigeration_compressor.eescode | Other | 2 | inf | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: Diverged - [Kinsol] KINSOL(Picard): residual be |
| scroll_compressor.eescode | Other | 34 | inf | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: Diverged - [Kinsol] KINSOL(Picard): r |
| simple_centrifugal_compressor.eescode | Other | 1 | 0.30 | Block 11 (size 1, vars: A) failed: MaxIterations - [Kinsol] KINSOL(Picard): max iterations (300) rea |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | inf | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: Diverged - [Kinsol] KINSOL(Pi |

### KINSOL fp/Anderson (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | inf | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(FP):  |
| air_screw_compressor_simple.eescode | Other | 13 | inf | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(FP):  |
| condenser_3zones.eescode | Other | 62 | inf | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: Diverged - [Kins |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(FP): re |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(FP): re |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | inf | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: Diverged - [Kinsol] KINSOL(FP): r |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: Diverged - [Kinsol] KINSOL(FP): residual bec |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 37: Construct 'module' is not yet handled by coolsolve   Line 196: Could not pa |
| orc_extraction.eescode | Other | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: Diverged - [Kinsol] KINS |
| orc_r245fa.eescode | Other | 8 | 5.6e+104 | Block 51 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Kinso |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: Diverged - [Kinsol] KINSOL(FP): resid |
| simple_centrifugal_compressor.eescode | Other | 1 | 5.5e+15 | Block 12 (size 1, vars: D) failed: MaxIterations - [Kinsol] KINSOL(FP): max iterations (300) reached |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | inf | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: Diverged - [Kinsol] KINSOL(FP |

### KINSOL linesearch (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | 177866.45 | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: LineSearchFailed - [Kinsol] KINS |
| air_screw_compressor_simple.eescode | Other | 13 | 22785.53 | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: LineSearchFailed - [Kinsol] KINS |
| condenser_3zones.eescode | Singular Jacobian | 62 | 254477.70 | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: SingularJacobian |
| cooling_tower.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian - [Kinsol] KINSOL |
| cooling_tower2.eescode | Singular Jacobian | 11 | 2.2e+07 | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian - [Kinsol] KINSOL |
| exchangers2.eescode | Singular Jacobian | 4 | 272048.58 | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: SingularJacobian - [Kinsol] KINSOL(LS): singula |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | 808669.19 | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Kinsol] KINSOL(L |
| internal_combustion_engine.eescode | Singular Jacobian | 5 | 399999.00 | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: SingularJacobian - [Kinsol] KINSOL( |
| internal_combustion_engine_cpbar.eescode | Other | 7 | 24081.66 | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [Kinsol] KINSOL(LS): max ite |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [Kinsol] KINSOL: INTERPOLATE(): lookup ta |
| orc_co2.eescode | Singular Jacobian | 28 | 491599.19 | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: SingularJacobian - [Kinsol] KINSOL(LS): sing |
| orc_extraction.eescode | Singular Jacobian | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: SingularJacobian - [Kins |
| orc_r245fa.eescode | Singular Jacobian | 12 | 64652.93 | Block 82 (size 12, vars: rho_hf_su_ev, nu_hf_su_ev, k_hf_su_ev, ...) failed: SingularJacobian - [Kin |
| scroll_compressor.eescode | Singular Jacobian | 34 | 5.2e+08 | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: SingularJacobian - [Kinsol] KINSOL(LS |
| turbocompressor.eescode | Singular Jacobian | 9 | 424823.39 | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: SingularJacobian - [Kinsol] KINSOL(LS): singula |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Singular Jacobian | 63 | 6196.80 | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: SingularJacobian - [Kinsol] K |

### KINSOL picard (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| advanced_features.eescode | Other | 2 | 3.1e+90 | Block 43 (size 2, vars: st_total, st_diff) failed: MaxIterations - [Kinsol] KINSOL(Picard): max iter |
| air_screw_compressor.eescode | Other | 13 | inf | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(Picar |
| air_screw_compressor_simple.eescode | Other | 13 | inf | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(Picar |
| condenser_3zones.eescode | Other | 62 | inf | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: Diverged - [Kins |
| condenser_wet.eescode | Other | 2 | inf | Block 19 (size 2, vars: M_dot_a, M_dot_cd) failed: Diverged - [Kinsol] KINSOL(Picard): residual beca |
| cooling_coil.eescode | Other | 12 | inf | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: Diverged - [Kinsol] KINSOL(Pi |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(Picard) |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(Picard) |
| cpbar.eescode | Other | 5 | 2.9e+97 | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [Kinsol] KINSOL(Picard): max iter |
| evaporator.eescode | Other | 2 | inf | Block 23 (size 2, vars: M_dot_a, M_dot_cd) failed: Diverged - [Kinsol] KINSOL(Picard): residual beca |
| exchangers2.eescode | Other | 4 | inf | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual be |
| exchangers3.eescode | Other | 3 | inf | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: Diverged - [Kinsol] KINSOL(Picard): residual b |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | inf | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: Diverged - [Kinsol] KINSOL(Picard |
| humidair2.eescode | Other | 5 | inf | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual |
| internal_combustion_engine.eescode | Other | 2 | inf | Block 16 (size 2, vars: h_ex_2, h_su_2) failed: Diverged - [Kinsol] KINSOL(Picard): residual became  |
| internal_combustion_engine_cpbar.eescode | Other | 7 | inf | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [Kinsol] KINSOL(Picard): INTERPOLATE(): l |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual |
| orc_extraction.eescode | Other | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: Diverged - [Kinsol] KINS |
| orc_r245fa.eescode | Other | 8 | inf | Block 51 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: Diverged - [Kinsol] KI |
| orc_simple.eescode | Other | 8 | — | Block 49 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: Diverged - [Kinsol] KI |
| piston_compressor.eescode | Other | 4 | 6.8e+91 | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: MaxIterations - [Kinsol] KINSOL(Pi |
| rankine2.eescode | Other | 4 | inf | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: Diverged - [Kinsol] KINSOL(P |
| refrigeration_compressor.eescode | Other | 2 | inf | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: Diverged - [Kinsol] KINSOL(Picard): residual be |
| scroll_compressor.eescode | Other | 34 | inf | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: Diverged - [Kinsol] KINSOL(Picard): r |
| simple_centrifugal_compressor.eescode | Other | 1 | 0.30 | Block 11 (size 1, vars: A) failed: MaxIterations - [Kinsol] KINSOL(Picard): max iterations (300) rea |
| turbocompressor.eescode | Other | 9 | inf | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: Diverged - [Kinsol] KINSOL(Picard): residual be |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | inf | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: Diverged - [Kinsol] KINSOL(Pi |

### KINSOL fp/Anderson (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | inf | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(FP):  |
| air_screw_compressor_simple.eescode | Other | 13 | inf | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: Diverged - [Kinsol] KINSOL(FP):  |
| condenser_3zones.eescode | Other | 62 | inf | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: Diverged - [Kins |
| cooling_coil.eescode | Other | 12 | 7.9e+16 | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [Kinsol] KINS |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(FP): re |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: Diverged - [Kinsol] KINSOL(FP): re |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | inf | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: Diverged - [Kinsol] KINSOL(FP): r |
| humidair2.eescode | Other | 5 | inf | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: Diverged - [Kinsol] KINSOL(FP): residual bec |
| internal_combustion_engine_cpbar.eescode | Other | 5 | inf | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: Diverged - [Kinsol] KINSOL(FP): residual became  |
| lookup_demo.eescode | Evaluation error | 1 | — | Block 1 (size 1, vars: h_interp) failed: EvaluationError - [Kinsol] KINSOL(FP): initial evaluation f |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: Diverged - [Kinsol] KINSOL(FP): residual bec |
| orc_extraction.eescode | Other | 21 | inf | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: Diverged - [Kinsol] KINS |
| orc_r245fa.eescode | Other | 6 | 4120.16 | Block 59 (size 6, vars: h_cf_ex_tp, T_cf_ex_tp, h_cf_ex_cd, ...) failed: MaxIterations - [Kinsol] KI |
| orc_simple.eescode | Other | 7 | 2.6e+26 | Block 75 (size 7, vars: h_hf_su_tp, T_hf_su_tp, h_hf_ex_tp, ...) failed: MaxIterations - [Kinsol] KI |
| scroll_compressor.eescode | Other | 34 | inf | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: Diverged - [Kinsol] KINSOL(FP): resid |
| simple_centrifugal_compressor.eescode | Other | 1 | 5.5e+15 | Block 12 (size 1, vars: D) failed: MaxIterations - [Kinsol] KINSOL(FP): max iterations (300) reached |
| turbocompressor.eescode | Other | 9 | inf | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: Diverged - [Kinsol] KINSOL(FP): residual became |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: 'lookup 1'. The first argument  |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |
| zorlu_heat_pump.eescode | Other | 63 | inf | Block 42 (size 63, vars: T_su_ev, T_sf_x1_ev, T_sf_x0_ev, ...) failed: Diverged - [Kinsol] KINSOL(FP |

