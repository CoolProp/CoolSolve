# CoolSolve Solver Robustness Report

**Generated:** Mon Mar  2 15:22:03 2026

Total example files tested: 39

**Legend:** OK = converged, PARSE = parse error, ANALYSIS = structural analysis error.  
Failure cells show: `ErrorCategory blkN |F|=residual` where N is the failed block size.

## Summary: Solve Success Rate by Configuration

| # | Configuration | Initials | Tearing | Solved | Total | Rate | Avg time (s) |
|---:|---|:---:|:---:|---:|---:|---:|---:|
| 1 | Default pipeline (with initials) | Yes | No | 35 | 38 | 92.1% | 0.171 |
| 2 | Newton only (with initials) | Yes | No | 31 | 38 | 81.6% | 0.020 |
| 3 | TrustRegion only (with initials) | Yes | No | 33 | 38 | 86.8% | 0.019 |
| 4 | LevenbergMarquardt only (with initials) | Yes | No | 29 | 38 | 76.3% | 0.032 |
| 5 | BisectionND only (with initials) | Yes | No | 15 | 38 | 39.5% | 0.005 |
| 6 | Homotopy only (with initials) | Yes | No | 31 | 38 | 81.6% | 0.043 |
| 7 | Partitioned only (with initials) | Yes | No | 17 | 38 | 44.7% | 0.074 |
| 8 | Default + Tearing (with initials) | Yes | Yes | 35 | 38 | 92.1% | 0.206 |
| 9 | Default pipeline (NO initials) | No | No | 27 | 37 | 73.0% | 0.280 |
| 10 | Newton only (NO initials) | No | No | 21 | 37 | 56.8% | 0.203 |
| 11 | TrustRegion only (NO initials) | No | No | 18 | 37 | 48.6% | 0.029 |
| 12 | LevenbergMarquardt only (NO initials) | No | No | 17 | 37 | 45.9% | 0.037 |
| 13 | BisectionND only (NO initials) | No | No | 13 | 37 | 35.1% | 0.005 |
| 14 | Homotopy only (NO initials) | No | No | 21 | 37 | 56.8% | 0.169 |
| 15 | Partitioned only (NO initials) | No | No | 14 | 37 | 37.8% | 0.088 |
| 16 | Default + Tearing (NO initials) | No | Yes | 28 | 37 | 75.7% | 0.246 |

## Detailed Results: With Initials

| File | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt | TR | LM | BisectionND | Homotopy | Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear |
|---|---|---|---|---|---|---|---|---|
| air_screw_compressor.eescode | **OK** (0.40s) | Max iterations blk13 | **OK** (0.03s) | Singular Jacobian blk13 | Other blk13 | Other blk13 | Max iterations blk13 | **OK** (0.17s) |
| air_screw_compressor_simple.eescode | **OK** (0.03s) | Singular Jacobian blk13 | **OK** (0.02s) | Singular Jacobian blk13 | Other blk13 | Other blk13 | Max iterations blk13 | **OK** (0.16s) |
| boiler_cpbar.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| boiler_cpbar2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) |
| compressor_refrigeration_simple.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| condenser_3zones.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk62 | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.05s) |
| condenser_wet.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk2 | **OK** (0.01s) |
| cooling_coil.eescode | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | Other blk12 | **OK** (0.08s) | **OK** (0.03s) | **OK** (0.13s) |
| cooling_tower.eescode | **OK** (0.18s) | **OK** (0.18s) | **OK** (0.22s) | Evaluation error blk11 | Other blk11 | **OK** (0.39s) | Max iterations blk11 | **OK** (0.83s) |
| cooling_tower2.eescode | **OK** (4.92s) | Singular Jacobian blk11 | Max iterations blk11 | Singular Jacobian blk11 | Other blk11 | Other blk11 | **OK** (1.15s) | **OK** (5.13s) |
| cpbar.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk5 | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) |
| evaporator.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk2 | **OK** (0.01s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.02s) | Other blk3 | **OK** (0.02s) | Other blk3 | **OK** (0.02s) |
| heat_pump_MSTh_SB_R10.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.42s) | Other blk39 | **OK** (0.01s) | Max iterations blk39 | **OK** (0.09s) |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.21s) | Other blk5 | **OK** (0.05s) | Max iterations blk5 | **OK** (0.06s) |
| internal_combustion_engine.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) |
| internal_combustion_engine_cpbar.eescode | **OK** (0.08s) | **OK** (0.08s) | Max iterations blk2 | **OK** (0.05s) | Other blk2 | **OK** (0.26s) | Other blk2 | **OK** (0.10s) |
| orc_co2.eescode | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | Other blk28 | **OK** (0.03s) | Max iterations blk28 | **OK** (0.06s) |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | Max iterations blk21 | Other blk21 | **OK** (0.03s) | Max iterations blk21 | **OK** (0.04s) |
| orc_r245fa.eescode | **OK** (0.04s) | **OK** (0.04s) | **OK** (0.04s) | **OK** (0.05s) | Other blk8 | **OK** (0.19s) | Other blk2 | **OK** (0.05s) |
| orc_simple.eescode | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | Other blk7 | **OK** (0.04s) | Other blk2 | **OK** (0.04s) |
| piston_compressor.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk2 | **OK** (0.00s) |
| pressuredrop.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| rankine1.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| rankine2.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk4 | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.01s) |
| refrigeration1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration_compressor.eescode | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk4 | **OK** (0.01s) | Other blk2 | **OK** (0.01s) |
| scroll_compressor.eescode | **OK** (0.11s) | **OK** (0.11s) | **OK** (0.10s) | Max iterations blk34 | Other blk34 | **OK** (0.12s) | Max iterations blk34 | **OK** (0.18s) |
| simple_centrifugal_compressor.eescode | **OK** (0.00s) | Singular Jacobian blk1 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk1 | Other blk1 | **OK** (0.00s) |
| storage_integraltable.eescode | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS |
| turbocompressor.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk9 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| turbocompressor_interpolate.eescode | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 |
| water_libr.eescode | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 |

## Detailed Results: Without Initials

| File | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt | TR | LM | BisectionND | Homotopy | Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear |
|---|---|---|---|---|---|---|---|---|
| air_screw_compressor.eescode | **OK** (0.42s) | Max iterations blk13 | **OK** (0.02s) | Singular Jacobian blk13 | Other blk13 | Other blk13 | Max iterations blk13 | **OK** (0.16s) |
| air_screw_compressor_simple.eescode | **OK** (0.03s) | Singular Jacobian blk13 | **OK** (0.02s) | Singular Jacobian blk13 | Other blk13 | Other blk13 | Max iterations blk13 | **OK** (0.16s) |
| boiler_cpbar.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| boiler_cpbar2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) |
| compressor_refrigeration_simple.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| condenser_3zones.eescode | **OK** (0.12s) | **OK** (0.12s) | Max iterations blk62 | Max iterations blk62 | Other blk62 | Other blk62 | Max iterations blk62 | **OK** (0.18s) |
| condenser_wet.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk2 | **OK** (0.01s) |
| cooling_coil.eescode | Other blk12 | Max iterations blk12 | Max iterations blk12 | Evaluation error blk12 | Other blk12 | Other blk12 | Max iterations blk12 | Other blk12 |
| cooling_tower.eescode | Max iterations blk11 | Singular Jacobian blk11 | Max iterations blk11 | Singular Jacobian blk11 | Other blk11 | Other blk11 | Max iterations blk11 | Max iterations blk11 |
| cooling_tower2.eescode | **OK** (4.99s) | Singular Jacobian blk11 | Max iterations blk11 | Singular Jacobian blk11 | Other blk11 | Other blk11 | **OK** (1.15s) | **OK** (5.10s) |
| cpbar.eescode | **OK** (0.01s) | **OK** (0.01s) | Max iterations blk5 | **OK** (0.01s) | Other blk5 | **OK** (0.03s) | Max iterations blk5 | **OK** (0.01s) |
| evaporator.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | Other blk2 | **OK** (0.01s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.03s) |
| exchangers3.eescode | **OK** (0.00s) | **OK** (0.00s) | Max iterations blk3 | Max iterations blk3 | Other blk3 | **OK** (0.01s) | Other blk3 | **OK** (0.02s) |
| heat_pump_MSTh_SB_R10.eescode | Max iterations blk39 | Max iterations blk39 | Max iterations blk39 | Max iterations blk39 | Other blk39 | Other blk39 | Max iterations blk39 | Max iterations blk39 |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | **OK** (0.01s) | **OK** (0.01s) | Max iterations blk5 | Max iterations blk5 | Other blk5 | Other blk5 | Max iterations blk5 | **OK** (0.03s) |
| internal_combustion_engine.eescode | **OK** (0.19s) | Singular Jacobian blk5 | Max iterations blk5 | Singular Jacobian blk5 | Other blk5 | Other blk5 | Other blk2 | **OK** (0.02s) |
| internal_combustion_engine_cpbar.eescode | Max iterations blk5 | Max iterations blk7 | Max iterations blk7 | Max iterations blk7 | Other blk7 | Other blk7 | Max iterations blk7 | Max iterations blk5 |
| orc_co2.eescode | Max iterations blk28 | Singular Jacobian blk28 | Max iterations blk28 | Singular Jacobian blk28 | Other blk28 | Other blk28 | Max iterations blk28 | Max iterations blk28 |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (0.81s) | Singular Jacobian blk21 | Max iterations blk21 | Singular Jacobian blk21 | Other blk21 | Other blk21 | **OK** (0.02s) | **OK** (0.15s) |
| orc_r245fa.eescode | **OK** (0.37s) | **OK** (0.36s) | Max iterations blk12 | **OK** (0.54s) | Other blk8 | **OK** (1.37s) | Other blk2 | **OK** (0.40s) |
| orc_simple.eescode | **OK** (0.30s) | **OK** (0.29s) | **OK** (0.30s) | Max iterations blk7 | Other blk8 | **OK** (1.85s) | Other blk2 | **OK** (0.31s) |
| piston_compressor.eescode | Max iterations blk4 | Max iterations blk4 | Singular Jacobian blk4 | Max iterations blk4 | Other blk4 | **OK** (0.00s) | Max iterations blk4 | **OK** (0.01s) |
| pressuredrop.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| rankine1.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| rankine2.eescode | **OK** (0.02s) | **OK** (0.02s) | Max iterations blk4 | **OK** (0.02s) | Other blk4 | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.02s) |
| refrigeration1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| refrigeration_compressor.eescode | Max iterations blk4 | Max iterations blk4 | Singular Jacobian blk4 | Max iterations blk4 | Other blk2 | **OK** (0.01s) | Other blk2 | Max iterations blk4 |
| scroll_compressor.eescode | Max iterations blk34 | **OK** (3.39s) | Max iterations blk34 | Singular Jacobian blk34 | Other blk34 | Other blk34 | Max iterations blk34 | Max iterations blk34 |
| simple_centrifugal_compressor.eescode | **OK** (0.00s) | Singular Jacobian blk1 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk1 | Other blk1 | **OK** (0.00s) |
| storage_integraltable.eescode | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS | ANALYSIS |
| turbocompressor.eescode | **OK** (0.24s) | Singular Jacobian blk9 | **OK** (0.14s) | Singular Jacobian blk9 | Other blk9 | **OK** (0.17s) | Max iterations blk9 | **OK** (0.25s) |
| turbocompressor_interpolate.eescode | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 | Evaluation error blk1 |
| water_libr.eescode | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 | Unsupported function blk1 |

## Model Difficulty Ranking

Models ranked by number of configurations that failed to solve them.

| File | Failures / Configs | Failed Configurations |
|---|---:|---|
| water_libr.eescode | 16 / 16 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| turbocompressor_interpolate.eescode | 16 / 16 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| internal_combustion_engine_cpbar.eescode | 11 / 16 | TrustRegion only (with initials), BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| cooling_tower.eescode | 11 / 16 | LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| air_screw_compressor_simple.eescode | 10 / 16 | Newton only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Newton only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| scroll_compressor.eescode | 10 / 16 | LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| cooling_tower2.eescode | 10 / 16 | Newton only (with initials), TrustRegion only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials) |
| orc_co2.eescode | 10 / 16 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| air_screw_compressor.eescode | 10 / 16 | Newton only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Newton only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| heat_pump_MSTh_SB_R10.eescode | 10 / 16 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| cooling_coil.eescode | 9 / 16 | BisectionND only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| refrigeration_compressor.eescode | 9 / 16 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| orc_complex.eescode | 8 / 8 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials) |
| orc_extraction.eescode | 8 / 16 | LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials) |
| humidair2.eescode | 7 / 16 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| internal_combustion_engine.eescode | 7 / 16 | Partitioned only (with initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| piston_compressor.eescode | 7 / 16 | Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| exchangers3.eescode | 6 / 16 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| condenser_3zones.eescode | 6 / 16 | BisectionND only (with initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| simple_centrifugal_compressor.eescode | 6 / 16 | Newton only (with initials), Homotopy only (with initials), Partitioned only (with initials), Newton only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| turbocompressor.eescode | 5 / 16 | BisectionND only (with initials), Newton only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| orc_r245fa.eescode | 5 / 16 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| orc_simple.eescode | 5 / 16 | BisectionND only (with initials), Partitioned only (with initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| cpbar.eescode | 4 / 16 | BisectionND only (with initials), TrustRegion only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| rankine2.eescode | 3 / 16 | BisectionND only (with initials), TrustRegion only (NO initials), BisectionND only (NO initials) |
| exchangers2.eescode | 2 / 16 | BisectionND only (with initials), BisectionND only (NO initials) |
| condenser_wet.eescode | 2 / 16 | Partitioned only (with initials), Partitioned only (NO initials) |
| evaporator.eescode | 2 / 16 | Partitioned only (with initials), Partitioned only (NO initials) |

## Error Category Breakdown

Across all configurations and models:

| Error Category | Count | Fraction |
|---|---:|---:|
| Evaluation error | 18 | 8.4% |
| Max iterations | 68 | 31.6% |
| Other | 88 | 40.9% |
| Singular Jacobian | 25 | 11.6% |
| Unsupported function | 16 | 7.4% |

## Detailed Error Messages

### Default pipeline (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Newton only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Max iterations | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Newton] Max ite |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| cooling_tower2.eescode | Singular Jacobian | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| simple_centrifugal_compressor.eescode | Singular Jacobian | 1 | — | Block 11 (size 1, vars: A) failed: SingularJacobian -  Initial //F//_inf = 0.295303, best achieved = |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### TrustRegion only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_tower2.eescode | Max iterations | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| internal_combustion_engine_cpbar.eescode | Max iterations | 2 | — | Block 39 (size 2, vars: t_7, c_p_g_67) failed: MaxIterations - [TrustRegion] Trust region: Max itera |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### LevenbergMarquardt only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Singular Jacobian | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| cooling_tower.eescode | Evaluation error | 11 | — | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: EvaluationError - [LevenbergMarq |
| cooling_tower2.eescode | Singular Jacobian | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| orc_extraction.eescode | Max iterations | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Levenbe |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [LevenbergMarquardt]  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
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
| cpbar.eescode | Other | 5 | — | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [BisectionND] BisectionND: max it |
| exchangers2.eescode | Other | 4 | — | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: MaxIterations - [BisectionND] BisectionND: max  |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [BisectionND] BisectionND: max |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: InvalidInput - [BisectionND] Bise |
| humidair2.eescode | Other | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| internal_combustion_engine_cpbar.eescode | Other | 2 | — | Block 39 (size 2, vars: t_7, c_p_g_67) failed: MaxIterations - [BisectionND] BisectionND: max iterat |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: InvalidInput - [BisectionND] BisectionND: bl |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: InvalidInput - [Bisectio |
| orc_r245fa.eescode | Other | 8 | — | Block 50 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| orc_simple.eescode | Other | 7 | — | Block 75 (size 7, vars: h_hf_su_tp, T_hf_su_tp, h_hf_ex_tp, ...) failed: MaxIterations - [BisectionN |
| rankine2.eescode | Other | 4 | — | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [BisectionND |
| refrigeration_compressor.eescode | Other | 4 | — | Block 43 (size 4, vars: C, epsilon_v_2, V_dot_s, ...) failed: MaxIterations - [BisectionND] Bisectio |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: InvalidInput - [BisectionND] Bisectio |
| turbocompressor.eescode | Other | 9 | — | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: InvalidInput - [BisectionND] BisectionND: block |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Homotopy only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Homotopy] Homot |
| air_screw_compressor_simple.eescode | Other | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Homotopy] Homot |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [Homotopy] Homotop |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| simple_centrifugal_compressor.eescode | Other | 1 | — | Block 11 (size 1, vars: A) failed: MaxIterations - [Homotopy] Homotopy: did not converge. Last t=1,  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Partitioned only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Max iterations | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| air_screw_compressor_simple.eescode | Max iterations | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| condenser_wet.eescode | Other | 2 | — | Block 19 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 178978, best |
| cooling_tower.eescode | Max iterations | 11 | — | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [Partitioned] Pa |
| evaporator.eescode | Other | 2 | — | Block 23 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 13869.6, bes |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations -  Initial //F//_inf = 310600, b |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Partitioned] Par |
| humidair2.eescode | Max iterations | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [Partitioned] Partitioned so |
| internal_combustion_engine.eescode | Other | 2 | — | Block 16 (size 2, vars: h_ex_2, h_su_2) failed: MaxIterations -  Initial //F//_inf = 8.9989e-08, bes |
| internal_combustion_engine_cpbar.eescode | Other | 2 | — | Block 38 (size 2, vars: t_6, c_p_g_6) failed: MaxIterations -  Initial //F//_inf = 8.24002e-10, best |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Partitioned] Partitioned so |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| orc_extraction.eescode | Max iterations | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Partiti |
| orc_r245fa.eescode | Other | 2 | — | Block 54 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 1377.92, best a |
| orc_simple.eescode | Other | 2 | — | Block 53 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 1.19908e-08, be |
| piston_compressor.eescode | Other | 2 | — | Block 27 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations -  Initial //F//_inf = 3.31784e-0 |
| refrigeration_compressor.eescode | Other | 2 | — | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations -  Initial //F//_inf = 2.50657e-0 |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Partitioned] Partiti |
| simple_centrifugal_compressor.eescode | Other | 1 | — | Block 11 (size 1, vars: A) failed: MaxIterations -  Initial //F//_inf = 0.295303, best achieved = 0. |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Default + Tearing (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Default pipeline (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_coil.eescode | Other | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [Newton] Solv |
| cooling_tower.eescode | Max iterations | 11 | — | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [TrustRegion] Tr |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Newton] Max iter |
| internal_combustion_engine_cpbar.eescode | Max iterations | 5 | — | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: MaxIterations - [Newton] Max iterations (100) re |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| piston_compressor.eescode | Max iterations | 4 | — | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: MaxIterations - [Newton] Max itera |
| refrigeration_compressor.eescode | Max iterations | 4 | — | Block 43 (size 4, vars: C, epsilon_v_2, V_dot_s, ...) failed: MaxIterations - [Newton] Max iteration |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Newton] Line search  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Newton only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Max iterations | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Newton] Max ite |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| cooling_coil.eescode | Max iterations | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [Newton] Max  |
| cooling_tower.eescode | Singular Jacobian | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| cooling_tower2.eescode | Singular Jacobian | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Newton] Max iter |
| internal_combustion_engine.eescode | Singular Jacobian | 5 | — | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: SingularJacobian -  Initial //F//_i |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | — | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [Newton] Max iterations (100 |
| orc_co2.eescode | Singular Jacobian | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: SingularJacobian -  Initial //F//_inf = 4573 |
| orc_extraction.eescode | Singular Jacobian | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: SingularJacobian |
| piston_compressor.eescode | Max iterations | 4 | — | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: MaxIterations - [Newton] Max itera |
| refrigeration_compressor.eescode | Max iterations | 4 | — | Block 43 (size 4, vars: C, epsilon_v_2, V_dot_s, ...) failed: MaxIterations - [Newton] Max iteration |
| simple_centrifugal_compressor.eescode | Singular Jacobian | 1 | — | Block 11 (size 1, vars: A) failed: SingularJacobian -  Initial //F//_inf = 0.295303, best achieved = |
| turbocompressor.eescode | Singular Jacobian | 9 | — | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: SingularJacobian -  Initial //F//_inf = 419405, |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### TrustRegion only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Max iterations | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| cooling_coil.eescode | Max iterations | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [TrustRegion] |
| cooling_tower.eescode | Max iterations | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| cooling_tower2.eescode | Max iterations | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [TrustRegion] Trus |
| cpbar.eescode | Max iterations | 5 | — | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [TrustRegion] Trust region: Max i |
| exchangers3.eescode | Max iterations | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [TrustRegion] Trust region: Ma |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [TrustRegion] Tru |
| humidair2.eescode | Max iterations | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| internal_combustion_engine.eescode | Max iterations | 5 | — | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: MaxIterations - [TrustRegion] Trust |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | — | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| orc_extraction.eescode | Max iterations | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [TrustRe |
| orc_r245fa.eescode | Max iterations | 12 | — | Block 81 (size 12, vars: rho_hf_su_ev, nu_hf_su_ev, k_hf_su_ev, ...) failed: MaxIterations - [TrustR |
| piston_compressor.eescode | Singular Jacobian | 4 | — | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: SingularJacobian -  Initial //F//_ |
| rankine2.eescode | Max iterations | 4 | — | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [TrustRegion |
| refrigeration_compressor.eescode | Singular Jacobian | 4 | — | Block 43 (size 4, vars: C, epsilon_v_2, V_dot_s, ...) failed: SingularJacobian -  Initial //F//_inf  |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [TrustRegion] Trust r |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### LevenbergMarquardt only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Singular Jacobian | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| air_screw_compressor_simple.eescode | Singular Jacobian | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: SingularJacobian -  Initial //F/ |
| condenser_3zones.eescode | Max iterations | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| cooling_coil.eescode | Evaluation error | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: EvaluationError - [LevenbergM |
| cooling_tower.eescode | Singular Jacobian | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| cooling_tower2.eescode | Singular Jacobian | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: SingularJacobian -  Initial //F//_ |
| exchangers3.eescode | Max iterations | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [LevenbergMarquardt] Levenberg |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [LevenbergMarquar |
| humidair2.eescode | Max iterations | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [LevenbergMarquardt] Levenbe |
| internal_combustion_engine.eescode | Singular Jacobian | 5 | — | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: SingularJacobian -  Initial //F//_i |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | — | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [LevenbergMarquardt] Levenbe |
| orc_co2.eescode | Singular Jacobian | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: SingularJacobian -  Initial //F//_inf = 4573 |
| orc_extraction.eescode | Singular Jacobian | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: SingularJacobian |
| orc_simple.eescode | Max iterations | 7 | — | Block 75 (size 7, vars: h_hf_su_tp, T_hf_su_tp, h_hf_ex_tp, ...) failed: MaxIterations - [LevenbergM |
| piston_compressor.eescode | Max iterations | 4 | — | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: MaxIterations - [LevenbergMarquard |
| refrigeration_compressor.eescode | Max iterations | 4 | — | Block 43 (size 4, vars: C, epsilon_v_2, V_dot_s, ...) failed: MaxIterations - [LevenbergMarquardt] L |
| scroll_compressor.eescode | Singular Jacobian | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: SingularJacobian -  Initial //F//_inf |
| turbocompressor.eescode | Singular Jacobian | 9 | — | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: SingularJacobian -  Initial //F//_inf = 419405, |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### BisectionND only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: InvalidInput - [BisectionND] Bis |
| air_screw_compressor_simple.eescode | Other | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: InvalidInput - [BisectionND] Bis |
| condenser_3zones.eescode | Other | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: InvalidInput - [ |
| cooling_coil.eescode | Other | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: InvalidInput - [BisectionND]  |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: InvalidInput - [BisectionND] Bisec |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: InvalidInput - [BisectionND] Bisec |
| cpbar.eescode | Other | 5 | — | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [BisectionND] BisectionND: max it |
| exchangers2.eescode | Other | 4 | — | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: MaxIterations - [BisectionND] BisectionND: max  |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [BisectionND] BisectionND: max |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: InvalidInput - [BisectionND] Bise |
| humidair2.eescode | Other | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| internal_combustion_engine.eescode | Other | 5 | — | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: MaxIterations - [BisectionND] Bisec |
| internal_combustion_engine_cpbar.eescode | Other | 7 | — | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: InvalidInput - [BisectionND] BisectionND: bl |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: InvalidInput - [Bisectio |
| orc_r245fa.eescode | Other | 8 | — | Block 50 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| orc_simple.eescode | Other | 8 | — | Block 49 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| piston_compressor.eescode | Other | 4 | — | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: MaxIterations - [BisectionND] Bise |
| rankine2.eescode | Other | 4 | — | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [BisectionND |
| refrigeration_compressor.eescode | Other | 2 | — | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations - [BisectionND] BisectionND: max  |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: InvalidInput - [BisectionND] Bisectio |
| turbocompressor.eescode | Other | 9 | — | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: InvalidInput - [BisectionND] BisectionND: block |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Homotopy only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Other | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Homotopy] Homot |
| air_screw_compressor_simple.eescode | Other | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Homotopy] Homot |
| condenser_3zones.eescode | Other | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| cooling_coil.eescode | Other | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [Homotopy] Ho |
| cooling_tower.eescode | Other | 11 | — | Block 14 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [Homotopy] Homotop |
| cooling_tower2.eescode | Other | 11 | — | Block 12 (size 11, vars: C_dot_max_f, omega_f, AU_f, ...) failed: MaxIterations - [Homotopy] Homotop |
| heat_pump_MSTh_SB_R10.eescode | Other | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Homotopy] Homoto |
| humidair2.eescode | Other | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [Homotopy] Homotopy: did not |
| internal_combustion_engine.eescode | Other | 5 | — | Block 35 (size 5, vars: T_thr_l, gamma_l, P_crit_l, ...) failed: MaxIterations - [Homotopy] Homotopy |
| internal_combustion_engine_cpbar.eescode | Other | 7 | — | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [Homotopy] Homotopy: did not |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Homotopy] Homotopy: did not |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Homotop |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Homotopy] Homotopy:  |
| simple_centrifugal_compressor.eescode | Other | 1 | — | Block 11 (size 1, vars: A) failed: MaxIterations - [Homotopy] Homotopy: did not converge. Last t=1,  |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Partitioned only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| air_screw_compressor.eescode | Max iterations | 13 | — | Block 13 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| air_screw_compressor_simple.eescode | Max iterations | 13 | — | Block 14 (size 13, vars: v_thr_leak, s_ex, h_thr_leak, ...) failed: MaxIterations - [Partitioned] Pa |
| condenser_3zones.eescode | Max iterations | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| condenser_wet.eescode | Other | 2 | — | Block 19 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 178978, best |
| cooling_coil.eescode | Max iterations | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - [Partitioned] |
| cooling_tower.eescode | Max iterations | 11 | — | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - [Partitioned] Pa |
| cpbar.eescode | Max iterations | 5 | — | Block 5 (size 5, vars: x, Q_4, e_min, ...) failed: MaxIterations - [Partitioned] Partitioned solver: |
| evaporator.eescode | Other | 2 | — | Block 23 (size 2, vars: M_dot_a, M_dot_cd) failed: MaxIterations -  Initial //F//_inf = 13869.6, bes |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations -  Initial //F//_inf = 4.38491e+ |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - [Partitioned] Par |
| humidair2.eescode | Max iterations | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [Partitioned] Partitioned so |
| internal_combustion_engine.eescode | Other | 2 | — | Block 16 (size 2, vars: h_ex_2, h_su_2) failed: MaxIterations -  Initial //F//_inf = 190817, best ac |
| internal_combustion_engine_cpbar.eescode | Max iterations | 7 | — | Block 22 (size 7, vars: v_4, p_4, W_dot_p, ...) failed: MaxIterations - [Partitioned] Partitioned so |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Partitioned] Partitioned so |
| orc_r245fa.eescode | Other | 2 | — | Block 54 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 76754.1, best a |
| orc_simple.eescode | Other | 2 | — | Block 53 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 76754.1, best a |
| piston_compressor.eescode | Max iterations | 4 | — | Block 19 (size 4, vars: epsilon_v_1, C, epsilon_v_2, ...) failed: MaxIterations - [Partitioned] Part |
| refrigeration_compressor.eescode | Other | 2 | — | Block 34 (size 2, vars: W_dot_loss_0, alpha) failed: MaxIterations -  Initial //F//_inf = 9483.31, b |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Partitioned] Partiti |
| simple_centrifugal_compressor.eescode | Other | 1 | — | Block 11 (size 1, vars: A) failed: MaxIterations -  Initial //F//_inf = 0.295303, best achieved = 0. |
| turbocompressor.eescode | Max iterations | 9 | — | Block 17 (size 9, vars: h_t_2, w, h_2s, ...) failed: MaxIterations - [Partitioned] Partitioned solve |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

### Default + Tearing (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| cooling_coil.eescode | Other | 12 | — | Block 35 (size 12, vars: T_cd, M_dot_cd, C_dot_min_wet_f, ...) failed: MaxIterations - Tearing: sing |
| cooling_tower.eescode | Max iterations | 11 | — | Block 15 (size 11, vars: t_wb_ex_r, h_a_ex_r, Q_dot_r, ...) failed: MaxIterations - Tearing: singula |
| heat_pump_MSTh_SB_R10.eescode | Max iterations | 39 | — | Block 13 (size 39, vars: Q_dot_rech, P_ex_cp, v_ex_1, ...) failed: MaxIterations - Tearing: singular |
| internal_combustion_engine_cpbar.eescode | Max iterations | 5 | — | Block 25 (size 5, vars: s_2, t_3, C_3, ...) failed: MaxIterations - Tearing: singular Schur compleme |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - Tearing: singular Schur comp |
| refrigeration_compressor.eescode | Max iterations | 4 | — | Block 43 (size 4, vars: C, epsilon_v_2, V_dot_s, ...) failed: MaxIterations - Tearing: max iteration |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - Tearing: singular Sch |
| turbocompressor_interpolate.eescode | Evaluation error | 1 | — | Block 4 (size 1, vars: M_r) failed: EvaluationError - Unknown fluid: lookup 1 |
| water_libr.eescode | Unsupported function | 1 | — | Block 16 (size 1, vars: x_gen) failed: EvaluationError - Unknown or unsupported function: X_LIBR wit |

