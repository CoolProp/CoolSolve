# CoolSolve Solver Robustness Report

**Generated:** Sun Mar  1 10:18:21 2026

Total example files tested: 18

**Legend:** OK = converged, PARSE = parse error, ANALYSIS = structural analysis error.  
Failure cells show: `ErrorCategory blkN |F|=residual` where N is the failed block size.

## Summary: Solve Success Rate by Configuration

| # | Configuration | Initials | Tearing | Solved | Total | Rate | Avg time (s) |
|---:|---|:---:|:---:|---:|---:|---:|---:|
| 1 | Default pipeline (with initials) | Yes | No | 17 | 18 | 94.4% | 0.064 |
| 2 | Newton only (with initials) | Yes | No | 17 | 18 | 94.4% | 0.064 |
| 3 | TrustRegion only (with initials) | Yes | No | 17 | 18 | 94.4% | 0.053 |
| 4 | LevenbergMarquardt only (with initials) | Yes | No | 15 | 18 | 83.3% | 0.143 |
| 5 | BisectionND only (with initials) | Yes | No | 7 | 18 | 38.9% | 0.009 |
| 6 | Homotopy only (with initials) | Yes | No | 17 | 18 | 94.4% | 0.248 |
| 7 | Partitioned only (with initials) | Yes | No | 10 | 18 | 55.6% | 0.015 |
| 8 | Default + Tearing (with initials) | Yes | Yes | 17 | 18 | 94.4% | 0.248 |
| 9 | Default pipeline (NO initials) | No | No | 14 | 17 | 82.4% | 0.570 |
| 10 | Newton only (NO initials) | No | No | 14 | 17 | 82.4% | 1.362 |
| 11 | TrustRegion only (NO initials) | No | No | 9 | 17 | 52.9% | 0.055 |
| 12 | LevenbergMarquardt only (NO initials) | No | No | 10 | 17 | 58.8% | 0.091 |
| 13 | BisectionND only (NO initials) | No | No | 7 | 17 | 41.2% | 0.008 |
| 14 | Homotopy only (NO initials) | No | No | 12 | 17 | 70.6% | 0.745 |
| 15 | Partitioned only (NO initials) | No | No | 10 | 17 | 58.8% | 0.018 |
| 16 | Default + Tearing (NO initials) | No | Yes | 14 | 17 | 82.4% | 0.658 |

## Detailed Results: With Initials

| File | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt | TR | LM | BisectionND | Homotopy | Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear |
|---|---|---|---|---|---|---|---|---|
| condenser_3zones.eescode | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | Other blk62 | **OK** (0.14s) | **OK** (0.05s) | **OK** (1.02s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers3.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (1.30s) | Other blk3 | **OK** (2.03s) | Other blk3 | **OK** (0.07s) |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.19s) | Other blk5 | **OK** (0.04s) | Evaluation error blk5 | **OK** (0.06s) |
| orc_co2.eescode | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | Other blk28 | **OK** (0.07s) | Max iterations blk28 | **OK** (0.18s) |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (0.09s) | **OK** (0.09s) | **OK** (0.09s) | Evaluation error blk21 | Other blk21 | **OK** (0.20s) | Max iterations blk21 | **OK** (0.29s) |
| orc_r245fa.eescode | **OK** (0.19s) | **OK** (0.19s) | **OK** (0.19s) | **OK** (0.34s) | Other blk8 | **OK** (0.59s) | Other blk2 | **OK** (0.26s) |
| orc_simple.eescode | **OK** (0.11s) | **OK** (0.11s) | **OK** (0.11s) | **OK** (0.11s) | Other blk7 | **OK** (0.17s) | Other blk2 | **OK** (0.13s) |
| pressuredrop.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| rankine1.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) |
| rankine2.eescode | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | Other blk4 | **OK** (0.03s) | **OK** (0.04s) | **OK** (0.03s) |
| refrigeration1.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration2.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| scroll_compressor.eescode | **OK** (0.47s) | **OK** (0.47s) | **OK** (0.28s) | Max iterations blk34 | Other blk34 | **OK** (0.87s) | Max iterations blk34 | **OK** (2.12s) |

## Detailed Results: Without Initials

| File | Nwt+TR+LM+BisectionND+Homotopy+Part | Nwt | TR | LM | BisectionND | Homotopy | Part | Nwt+TR+LM+BisectionND+Homotopy+Part+Tear |
|---|---|---|---|---|---|---|---|---|
| condenser_3zones.eescode | **OK** (3.41s) | **OK** (3.53s) | Max iterations blk62 | Evaluation error blk62 | Other blk62 | Other blk62 | Max iterations blk62 | **OK** (5.11s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.00s) | **OK** (0.00s) | Other blk4 | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.02s) |
| exchangers3.eescode | **OK** (0.03s) | **OK** (0.03s) | Max iterations blk3 | Max iterations blk3 | Other blk3 | **OK** (0.56s) | Other blk3 | **OK** (0.92s) |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | Other blk5 | Evaluation error blk5 | Evaluation error blk5 | Evaluation error blk5 | Other blk5 | Evaluation error blk5 | Evaluation error blk5 | Other blk5 |
| orc_co2.eescode | Max iterations blk28 | Singular Jacobian blk28 | Max iterations blk28 | Singular Jacobian blk28 | Other blk28 | Other blk28 | Max iterations blk28 | Max iterations blk28 |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (3.46s) | Singular Jacobian blk21 | Max iterations blk21 | Singular Jacobian blk21 | Other blk21 | Other blk21 | **OK** (0.07s) | **OK** (1.93s) |
| orc_r245fa.eescode | **OK** (0.53s) | **OK** (0.54s) | Max iterations blk12 | **OK** (0.82s) | Other blk8 | **OK** (1.75s) | Other blk2 | **OK** (0.66s) |
| orc_simple.eescode | **OK** (0.44s) | **OK** (0.44s) | **OK** (0.44s) | Max iterations blk7 | Other blk8 | **OK** (6.53s) | Other blk2 | **OK** (0.47s) |
| pressuredrop.eescode | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| rankine1.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) |
| rankine2.eescode | **OK** (0.04s) | **OK** (0.04s) | Max iterations blk4 | **OK** (0.04s) | Other blk4 | **OK** (0.04s) | **OK** (0.05s) | **OK** (0.04s) |
| refrigeration1.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration2.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| scroll_compressor.eescode | Other blk34 | **OK** (14.43s) | Max iterations blk34 | Singular Jacobian blk34 | Other blk34 | Other blk34 | Max iterations blk34 | Other blk34 |

## Model Difficulty Ranking

Models ranked by number of configurations that failed to solve them.

| File | Failures / Configs | Failed Configurations |
|---|---:|---|
| scroll_compressor.eescode | 10 / 16 | LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| humidair2.eescode | 10 / 16 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| orc_co2.eescode | 10 / 16 | BisectionND only (with initials), Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| orc_complex.eescode | 8 / 8 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), LevenbergMarquardt only (with initials), BisectionND only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials) |
| orc_extraction.eescode | 8 / 16 | LevenbergMarquardt only (with initials), BisectionND only (with initials), Partitioned only (with initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials) |
| exchangers3.eescode | 6 / 16 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| condenser_3zones.eescode | 6 / 16 | BisectionND only (with initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| orc_r245fa.eescode | 5 / 16 | BisectionND only (with initials), Partitioned only (with initials), TrustRegion only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| orc_simple.eescode | 5 / 16 | BisectionND only (with initials), Partitioned only (with initials), LevenbergMarquardt only (NO initials), BisectionND only (NO initials), Partitioned only (NO initials) |
| rankine2.eescode | 3 / 16 | BisectionND only (with initials), TrustRegion only (NO initials), BisectionND only (NO initials) |
| exchangers2.eescode | 2 / 16 | BisectionND only (with initials), BisectionND only (NO initials) |

## Error Category Breakdown

Across all configurations and models:

| Error Category | Count | Fraction |
|---|---:|---:|
| Evaluation error | 8 | 11.0% |
| Max iterations | 18 | 24.7% |
| Other | 42 | 57.5% |
| Singular Jacobian | 5 | 6.8% |

## Detailed Error Messages

### Default pipeline (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |

### Newton only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |

### TrustRegion only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |

### LevenbergMarquardt only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| orc_extraction.eescode | Evaluation error | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: EvaluationError - [Leven |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [LevenbergMarquardt]  |

### BisectionND only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Other | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: InvalidInput - [ |
| exchangers2.eescode | Other | 4 | — | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: MaxIterations - [BisectionND] BisectionND: max  |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [BisectionND] BisectionND: max |
| humidair2.eescode | Other | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: InvalidInput - [BisectionND] BisectionND: sy |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: InvalidInput - [Bisectio |
| orc_r245fa.eescode | Other | 8 | — | Block 50 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| orc_simple.eescode | Other | 7 | — | Block 75 (size 7, vars: h_hf_su_tp, T_hf_su_tp, h_hf_ex_tp, ...) failed: MaxIterations - [BisectionN |
| rankine2.eescode | Other | 4 | — | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [BisectionND |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: InvalidInput - [BisectionND] Bisectio |

### Homotopy only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |

### Partitioned only (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations -  Initial //F//_inf = 310600, b |
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - [Partitioned] Partitioned  |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Partitioned] Partitioned so |
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |
| orc_extraction.eescode | Max iterations | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Partiti |
| orc_r245fa.eescode | Other | 2 | — | Block 54 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 1377.92, best a |
| orc_simple.eescode | Other | 2 | — | Block 53 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 1.19908e-08, be |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Partitioned] Partiti |

### Default + Tearing (with initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| orc_complex.eescode | Other | ? | — | Parse failed:   Line 34: Construct 'module' is not yet handled by coolsolve   Line 193: Could not pa |

### Default pipeline (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| humidair2.eescode | Other | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [Newton] Solver timed out [T |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Newton] Solver timed |

### Newton only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - [Newton] HumidAir error: H |
| orc_co2.eescode | Singular Jacobian | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: SingularJacobian -  Initial //F//_inf = 4573 |
| orc_extraction.eescode | Singular Jacobian | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: SingularJacobian |

### TrustRegion only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Max iterations | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| exchangers3.eescode | Max iterations | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [TrustRegion] Trust region: Ma |
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - [TrustRegion] HumidAir err |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [TrustRegion] Trust region:  |
| orc_extraction.eescode | Max iterations | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [TrustRe |
| orc_r245fa.eescode | Max iterations | 12 | — | Block 81 (size 12, vars: rho_hf_su_ev, nu_hf_su_ev, k_hf_su_ev, ...) failed: MaxIterations - [TrustR |
| rankine2.eescode | Max iterations | 4 | — | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [TrustRegion |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [TrustRegion] Trust r |

### LevenbergMarquardt only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Evaluation error | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: EvaluationError  |
| exchangers3.eescode | Max iterations | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [LevenbergMarquardt] Levenberg |
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - [LevenbergMarquardt] Humid |
| orc_co2.eescode | Singular Jacobian | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: SingularJacobian -  Initial //F//_inf = 4573 |
| orc_extraction.eescode | Singular Jacobian | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: SingularJacobian |
| orc_simple.eescode | Max iterations | 7 | — | Block 75 (size 7, vars: h_hf_su_tp, T_hf_su_tp, h_hf_ex_tp, ...) failed: MaxIterations - [LevenbergM |
| scroll_compressor.eescode | Singular Jacobian | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: SingularJacobian -  Initial //F//_inf |

### BisectionND only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Other | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: InvalidInput - [ |
| exchangers2.eescode | Other | 4 | — | Block 25 (size 4, vars: U, h_w, T_wall, ...) failed: MaxIterations - [BisectionND] BisectionND: max  |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations - [BisectionND] BisectionND: max |
| humidair2.eescode | Other | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - [BisectionND] BisectionND: m |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: InvalidInput - [BisectionND] BisectionND: sy |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: InvalidInput - [Bisectio |
| orc_r245fa.eescode | Other | 8 | — | Block 50 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| orc_simple.eescode | Other | 8 | — | Block 49 (size 8, vars: DELTAp_vap_rec, p_ex_vap_rec, p_vap_rec, ...) failed: MaxIterations - [Bisec |
| rankine2.eescode | Other | 4 | — | Block 28 (size 4, vars: W_dot_t_2, W_dot_t_1, M_dot_steam, ...) failed: MaxIterations - [BisectionND |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: InvalidInput - [BisectionND] Bisectio |

### Homotopy only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Other | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - [Homotopy] HumidAir error: |
| orc_co2.eescode | Other | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Homotopy] Homotopy: did not |
| orc_extraction.eescode | Other | 21 | — | Block 40 (size 21, vars: P_f_in3_exp, w_in_4_exp, v_f_in3_exp, ...) failed: MaxIterations - [Homotop |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Homotopy] Homotopy:  |

### Partitioned only (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| condenser_3zones.eescode | Max iterations | 62 | — | Block 38 (size 62, vars: t_r_ex_cd_sh, M_dot_cf_cd_sh, H_dot_cf_cd_sh, ...) failed: MaxIterations -  |
| exchangers3.eescode | Other | 3 | — | Block 7 (size 3, vars: T_w_ex, T_w_bar, cp_w) failed: MaxIterations -  Initial //F//_inf = 4.38491e+ |
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - [Partitioned] Partitioned  |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - [Partitioned] Partitioned so |
| orc_r245fa.eescode | Other | 2 | — | Block 54 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 76754.1, best a |
| orc_simple.eescode | Other | 2 | — | Block 53 (size 2, vars: p_cd_v, h_cd_v) failed: MaxIterations -  Initial //F//_inf = 76754.1, best a |
| scroll_compressor.eescode | Max iterations | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - [Partitioned] Partiti |

### Default + Tearing (NO initials)

| File | Category | Block | Residual | Error (truncated) |
|---|---|---:|---:|---|
| humidair2.eescode | Other | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: MaxIterations - Tearing evaluation failed: H |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - Tearing: singular Schur comp |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - Tearing: singular Sch |

