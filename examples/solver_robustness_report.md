# CoolSolve Solver Robustness Report

**Generated:** Sat Feb 28 22:21:09 2026

Total example files tested: 18

**Legend:** OK = converged, PARSE = parse error, ANALYSIS = structural analysis error.  
Failure cells show: `ErrorCategory blkN |F|=residual` where N is the failed block size.

## Summary: Solve Success Rate by Configuration

| # | Configuration | Initials | Tearing | Solved | Total | Rate | Avg time (s) |
|---:|---|:---:|:---:|---:|---:|---:|---:|
| 1 | Default pipeline (with initials) | Yes | No | 17 | 18 | 94.4% | 0.066 |
| 2 | Newton only (with initials) | Yes | No | 17 | 18 | 94.4% | 0.067 |
| 3 | TrustRegion only (with initials) | Yes | No | 17 | 18 | 94.4% | 0.054 |
| 4 | LevenbergMarquardt only (with initials) | Yes | No | 15 | 18 | 83.3% | 0.145 |
| 5 | Homotopy only (with initials) | Yes | No | 17 | 18 | 94.4% | 0.254 |
| 6 | Partitioned only (with initials) | Yes | No | 10 | 18 | 55.6% | 0.015 |
| 7 | Default + Tearing (with initials) | Yes | Yes | 17 | 18 | 94.4% | 0.252 |
| 8 | Default pipeline (NO initials) | No | No | 14 | 17 | 82.4% | 0.575 |
| 9 | Newton only (NO initials) | No | No | 14 | 17 | 82.4% | 1.380 |
| 10 | TrustRegion only (NO initials) | No | No | 9 | 17 | 52.9% | 0.057 |
| 11 | LevenbergMarquardt only (NO initials) | No | No | 10 | 17 | 58.8% | 0.095 |
| 12 | Homotopy only (NO initials) | No | No | 12 | 17 | 70.6% | 0.758 |
| 13 | Partitioned only (NO initials) | No | No | 10 | 17 | 58.8% | 0.019 |
| 14 | Default + Tearing (NO initials) | No | Yes | 14 | 17 | 82.4% | 0.691 |

## Detailed Results: With Initials

| File | Nwt+TR+LM+Homotopy+Part | Nwt | TR | LM | Homotopy | Part | Nwt+TR+LM+Homotopy+Part+Tear |
|---|---|---|---|---|---|---|---|
| condenser_3zones.eescode | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.14s) | **OK** (0.05s) | **OK** (1.04s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers3.eescode | **OK** (0.02s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (1.32s) | **OK** (1.98s) | Other blk3 | **OK** (0.07s) |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.20s) | **OK** (0.04s) | Evaluation error blk5 | **OK** (0.06s) |
| orc_co2.eescode | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.05s) | **OK** (0.08s) | Max iterations blk28 | **OK** (0.19s) |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (0.09s) | **OK** (0.09s) | **OK** (0.09s) | Evaluation error blk21 | **OK** (0.21s) | Max iterations blk21 | **OK** (0.30s) |
| orc_r245fa.eescode | **OK** (0.20s) | **OK** (0.20s) | **OK** (0.19s) | **OK** (0.35s) | **OK** (0.61s) | Other blk2 | **OK** (0.26s) |
| orc_simple.eescode | **OK** (0.11s) | **OK** (0.13s) | **OK** (0.11s) | **OK** (0.11s) | **OK** (0.18s) | Other blk2 | **OK** (0.13s) |
| pressuredrop.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| rankine1.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) |
| rankine2.eescode | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.03s) | **OK** (0.04s) | **OK** (0.04s) | **OK** (0.03s) |
| refrigeration1.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration2.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| scroll_compressor.eescode | **OK** (0.48s) | **OK** (0.48s) | **OK** (0.29s) | Max iterations blk34 | **OK** (0.98s) | Max iterations blk34 | **OK** (2.14s) |

## Detailed Results: Without Initials

| File | Nwt+TR+LM+Homotopy+Part | Nwt | TR | LM | Homotopy | Part | Nwt+TR+LM+Homotopy+Part+Tear |
|---|---|---|---|---|---|---|---|
| condenser_3zones.eescode | **OK** (3.50s) | **OK** (3.54s) | Max iterations blk62 | Evaluation error blk62 | Other blk62 | Max iterations blk62 | **OK** (5.29s) |
| exchangers1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| exchangers2.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.01s) | **OK** (0.00s) | **OK** (0.03s) |
| exchangers3.eescode | **OK** (0.03s) | **OK** (0.03s) | Max iterations blk3 | Max iterations blk3 | **OK** (0.56s) | Other blk3 | **OK** (1.02s) |
| humidair1.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| humidair2.eescode | Evaluation error blk5 | Evaluation error blk5 | Evaluation error blk5 | Evaluation error blk5 | Evaluation error blk5 | Evaluation error blk5 | Evaluation error blk5 |
| orc_co2.eescode | Max iterations blk28 | Singular Jacobian blk28 | Max iterations blk28 | Singular Jacobian blk28 | Other blk28 | Max iterations blk28 | Max iterations blk28 |
| orc_complex.eescode | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE | PARSE |
| orc_extraction.eescode | **OK** (3.41s) | Singular Jacobian blk21 | Max iterations blk21 | Singular Jacobian blk21 | Other blk21 | **OK** (0.08s) | **OK** (2.04s) |
| orc_r245fa.eescode | **OK** (0.55s) | **OK** (0.55s) | Max iterations blk12 | **OK** (0.85s) | **OK** (1.79s) | Other blk2 | **OK** (0.70s) |
| orc_simple.eescode | **OK** (0.44s) | **OK** (0.45s) | **OK** (0.45s) | Max iterations blk7 | **OK** (6.63s) | Other blk2 | **OK** (0.50s) |
| pressuredrop.eescode | **OK** (0.02s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| rankine1.eescode | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) | **OK** (0.02s) |
| rankine2.eescode | **OK** (0.04s) | **OK** (0.04s) | Max iterations blk4 | **OK** (0.04s) | **OK** (0.04s) | **OK** (0.05s) | **OK** (0.04s) |
| refrigeration1.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration2.eescode | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) | **OK** (0.01s) |
| refrigeration3.eescode | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) | **OK** (0.00s) |
| scroll_compressor.eescode | Other blk34 | **OK** (14.62s) | Max iterations blk34 | Singular Jacobian blk34 | Other blk34 | Max iterations blk34 | Other blk34 |

## Model Difficulty Ranking

Models ranked by number of configurations that failed to solve them.

| File | Failures / Configs | Failed Configurations |
|---|---:|---|
| scroll_compressor.eescode | 8 / 14 | LevenbergMarquardt only (with initials), Partitioned only (with initials), Default pipeline (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| humidair2.eescode | 8 / 14 | Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| orc_co2.eescode | 8 / 14 | Partitioned only (with initials), Default pipeline (NO initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials), Default + Tearing (NO initials) |
| orc_complex.eescode | 7 / 7 | Default pipeline (with initials), Newton only (with initials), TrustRegion only (with initials), LevenbergMarquardt only (with initials), Homotopy only (with initials), Partitioned only (with initials), Default + Tearing (with initials) |
| orc_extraction.eescode | 6 / 14 | LevenbergMarquardt only (with initials), Partitioned only (with initials), Newton only (NO initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), Homotopy only (NO initials) |
| exchangers3.eescode | 4 / 14 | Partitioned only (with initials), TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), Partitioned only (NO initials) |
| condenser_3zones.eescode | 4 / 14 | TrustRegion only (NO initials), LevenbergMarquardt only (NO initials), Homotopy only (NO initials), Partitioned only (NO initials) |
| orc_r245fa.eescode | 3 / 14 | Partitioned only (with initials), TrustRegion only (NO initials), Partitioned only (NO initials) |
| orc_simple.eescode | 3 / 14 | Partitioned only (with initials), LevenbergMarquardt only (NO initials), Partitioned only (NO initials) |
| rankine2.eescode | 1 / 14 | TrustRegion only (NO initials) |

## Error Category Breakdown

Across all configurations and models:

| Error Category | Count | Fraction |
|---|---:|---:|
| Evaluation error | 10 | 19.2% |
| Max iterations | 18 | 34.6% |
| Other | 19 | 36.5% |
| Singular Jacobian | 5 | 9.6% |

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
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - [Newton] HumidAir error: H |
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
| humidair2.eescode | Evaluation error | 5 | — | Block 8 (size 5, vars: w_su, v_a_su, T_su, ...) failed: EvaluationError - Tearing evaluation failed: |
| orc_co2.eescode | Max iterations | 28 | — | Block 77 (size 28, vars: T[5], T[4], T[3], ...) failed: MaxIterations - Tearing: singular Schur comp |
| scroll_compressor.eescode | Other | 34 | — | Block 32 (size 34, vars: W_dot_loss, w_nad, w_ad, ...) failed: MaxIterations - Tearing: singular Sch |

