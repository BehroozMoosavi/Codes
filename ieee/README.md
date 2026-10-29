# Volatility-Adaptive Elastic Net for CAISO Load Forecasting

This repository contains the code, data, and verified results for the revised IEEE conference study:

**Volatility-Adaptive Elastic Net for CAISO Load Forecasting**

**Author:** Behrooz Moosavi Ramezanzadeh

## Overview

The project studies one-hour-ahead forecasting of CAISO hourly system demand using rolling penalized regression.

The proposed method, Volatility-Adaptive Rolling Elastic Net (VA-REN), allows the elastic-net mixing parameter to vary with recent load volatility.

The revised experiment uses a corrected causal forecasting design and CAISO OASIS hourly demand data with provider-supplied timestamps.

## Repository structure

### `matlab/`

MATLAB source code for:

- CAISO source validation and preprocessing
- corrected rolling elastic-net forecasting
- Ridge, LASSO, Static Elastic Net, Fixed Elastic Net, and VA-REN evaluation
- LSBoost nonlinear benchmarking
- checkpointed full-sample runs
- forecast-comparison inference
- numerical validation tests

Main scripts include:

- `prepare_caiso_oasis_inputs.m`
- `run_varen_step4_pilot.m`
- `run_varen_step5_full.m`
- `run_varen_step6_nonlinear_pilot.m`
- `run_varen_step6b_calibration.m`
- `run_varen_step7_full.m`
- `run_varen_step8_analysis.m`

### `data/processed/`

Processed CAISO OASIS data used in the revised empirical analysis.

Main files:

- `caiso_oasis_actual_2019_2023.csv`
- `caiso_oasis_inputs_2019_2023.mat`
- `source_manifest.csv`
- validation tables for yearly, daily, and archive checks

The dataset contains 43,824 unique contiguous hourly system-demand observations covering 2019 through 2023.

### `data/provenance/`

Source metadata and validation information for the CAISO OASIS archives.

This folder contains archive-level metadata, the source manifest, source-validation results, and the Step 3 audit summary.

### `results/`

Verified outputs from the revised forecasting experiment.

Important files include:

- `step5_full_results.zip` — full linear-model forecasting results
- `step7_full_results.zip` — full LSBoost benchmark results
- `step8_analysis_results.zip` — final statistical inference
- `step8_primary_comparisons.csv`
- `step8_dm_tests.csv`
- `step8_gw_joint_tests.csv`
- `step8_effect_intervals.csv`
- `step8_volatility_diagnostics.csv`
- `step8_model_timing.csv`

## Main empirical results

The common evaluation sample contains **34,343 one-hour-ahead forecasts** from January 2020 through December 2023.

| Model | RMSE (MW) | MAE (MW) | MAPE (%) |
|---|---:|---:|---:|
| Ridge | 236.365 | 172.656 | 0.6932 |
| LASSO / calibrated Static EN | 234.616 | 171.090 | 0.6876 |
| Fixed EN, alpha = 0.5 | 234.482 | 170.988 | 0.6873 |
| VA-REN | 234.483 | 170.996 | 0.6872 |
| LSBoost24h | 248.939 | 183.081 | 0.7362 |

VA-REN improves RMSE relative to calibrated Static EN by approximately **0.133 MW**, but the difference is not statistically significant.

VA-REN and Fixed EN with alpha = 0.5 perform essentially identically.

Relative to the tested LSBoost benchmark, VA-REN lowers RMSE by approximately **14.46 MW**, or **5.81%**.

## Statistical inference

The final analysis uses dependence-robust forecast-comparison procedures.

For the primary comparisons:

- VA-REN vs. Static EN: no statistically detectable improvement
- VA-REN vs. Fixed EN, alpha = 0.5: no statistically detectable difference
- VA-REN vs. LSBoost24h: materially and statistically lower squared forecast loss

See:

- `results/step8_primary_comparisons.csv`
- `results/step8_dm_tests.csv`
- `results/step8_gw_joint_tests.csv`
- `results/step8_effect_intervals.csv`

## Software

The revised analysis was run using:

- MATLAB R2025a
- Statistics and Machine Learning Toolbox

## Reproduction order

A full reproduction should follow:

1. `prepare_caiso_oasis_inputs.m`
2. `run_varen_step4_pilot.m`
3. `run_varen_step5_full.m`
4. `run_varen_step6_nonlinear_pilot.m`
5. `run_varen_step6b_calibration.m`
6. `run_varen_step7_full.m`
7. `run_varen_step8_analysis.m`

The full linear and nonlinear rolling experiments are computationally intensive.

## Data note

The forecasting target is CAISO **system demand**, not solar-subtracted net load.

The revised implementation uses only information available before each forecast target when constructing the volatility measure and training the forecasting models.

Published historical CAISO actuals may contain later revisions; the experiment assumes completed-hour demand is available at the next forecast origin.

## Code and results

Project repository:

https://github.com/BehroozMoosavi/Codes/tree/main/ieee

MATLAB code:

https://github.com/BehroozMoosavi/Codes/tree/main/ieee/matlab

Data:

https://github.com/BehroozMoosavi/Codes/tree/main/ieee/data

Results:

https://github.com/BehroozMoosavi/Codes/tree/main/ieee/results
