============================================================
  VA-REN PROJECT  —  Complete File Package
  30 files, flat folder structure
============================================================

WHAT THIS IS
------------
Complete MATLAB codebase for the VA-REN paper:
"Volatility-Adaptive Rolling Elastic Net for Short-Term
Electricity Demand Forecasting"

All files go in ONE folder alongside the data files.
No subfolders needed.

============================================================
RUN ORDER
============================================================

STEP 1 — Build the dataset  (run once)
   >> data_clean.m
   Reads 5 Excel files, builds all predictors, computes
   volatility index, standardizes, saves the .mat dataset.
   Required Excel files (must be in same folder):
     historicalemshourlyload-2019.xlsx
     historicalemshourlyload-2020.xlsx
     historicalemshourlyload-2021.xlsx
     historicalemshourlyload-2022.xlsx
     historicalemshourlyloadfor2023.xlsx
   Output: caiso_final_logchange_dataset_2019_2023.mat

STEP 2 — Run the models  (can run in parallel)
   >> naive.m          (runs in seconds)
   >> ridge.m          (runs in ~1-2 hours)
   >> lasso_method.m   (runs in ~2-3 hours)
   >> static_en.m      (runs in ~3-4 hours)
   >> va_ren.m         (runs in ~3-4 hours)
   Each saves a results_*.mat file.

STEP 3 — Run the full empirical analysis
   >> master_empirical.m
   Runs all tests, robustness checks, figures, and tables.
   Outputs go to output/tables/ and output/figures/

   Or run individual pieces:
   >> comparison_tables.m   quick LaTeX comparison
   >> tabulation.m          full 6-table paper output

============================================================
FILE INVENTORY  (30 files)
============================================================

CORE MODELS
  01 data_clean.m          build dataset
  02 naive.m               seasonal naive baseline
  03 ridge.m               rolling ridge regression
  04 lasso_method.m        rolling LASSO
  05 static_en.m           static elastic net (primary comparison)
  06 va_ren.m              VA-REN (proposed method)

TABLE WRITERS
  07 comparison_tables.m   quick Overleaf comparison (4 tables)
  08 tabulation.m          full paper tables (6 tables + CSVs)

SHARED UTILITIES
  09 newey_west_hac.m      HAC covariance estimator
  10 make_blocked_folds.m  temporally-ordered CV folds
  11 align_errors.m        timestamp-align two error vectors
  12 back_transform.m      z_std -> z_raw -> MW level
  13 significance_stars.m  * ** *** for p-values
  14 latex_escape.m        escape _ % & for LaTeX
  15 save_fig.m            save PDF + PNG

STATISTICAL TESTS
  16 dm_test.m             one-sided Diebold-Mariano test
  17 dm_regime.m           regime-specific DM (stable + ramp)
  18 mcs_test.m            model confidence set
  19 gw_test.m             Giacomini-White conditional test
  20 regime_selection.m    predictor selection frequency analysis

VISUALIZATION
  21 plot_all_results.m    5 publication figures
  22 plot_regime_selection.m  3 regime selection figures

LATEX WRITER
  23 write_selection_latex.m  table_7a + table_7b

ROBUSTNESS CHECKS
  24 run_varen_window.m    R1: window length sensitivity
  25 run_varen_volwindow.m R2: volatility window sensitivity
  26 run_varen_alphamap.m  R3: alpha mapping alternatives
  27 run_subperiod.m       R4: yearly subperiod stability
  28 run_adaptive_lasso.m  R5: adaptive LASSO baseline

MASTER + README
  29 master_empirical.m   runs everything after models are done
  30 README.txt           this file

============================================================
PAPER TABLES PRODUCED
============================================================

From tabulation.m and write_selection_latex.m:
  table_1_data_summary.tex/.csv
  table_2_main_accuracy.tex/.csv
  table_3_regime_accuracy.tex/.csv
  table_4_dm_tests.tex/.csv
  table_5_model_diagnostics.tex/.csv
  table_6_varen_mechanism.tex/.csv
  table_7a_selection_top10.tex
  table_7b_selection_groups.tex

Include in paper with:  \input{table_2_main_accuracy.tex}

============================================================
FIGURES PRODUCED
============================================================

From plot_all_results.m:
  fig1_volatility           alpha_t and V_t time series
  fig2_monthly_rmse         monthly RMSE comparison
  fig3_regime_bar           stable vs ramp RMSE bar chart
  fig4_window_robust        RMSE vs window W
  fig5_alpha_distribution   alpha_t histogram

From plot_regime_selection.m:
  figA_selection_diverging  all 70 predictors by Delta_j
  figB_selection_group      group-level bar chart
  figC_selection_heatmap    group x hour-of-day heatmap

============================================================
REQUIREMENTS
============================================================

  MATLAB R2020a or later
  Statistics and Machine Learning Toolbox
  (for lasso(), movvar(), exportgraphics)

  Estimated run time for full pipeline:
    data_clean.m   :  5-10 min
    5 model scripts:  10-15 hours total (can run in parallel)
    master_empirical:  4-6 hours (robustness checks dominate)

============================================================
HYPOTHESES TESTED
============================================================

  H1: VA-REN RMSE < Static EN RMSE  (full sample, DM test)
  H2: Ramp-hour DM > Stable-hour DM (regime concentration)
  H3a: mean_alpha_ramp < mean_alpha_stable
  H3b: mean_vol_ramp > mean_vol_stable
  H3c: Corr(alpha_t, active_set_size) > 0
  H4: VA-REN > Static EN > LASSO > Ridge > Seasonal Naive

  GW test: gamma_1 > 0  (advantage grows with volatility)
  MCS: VA-REN survives at alpha=0.10

============================================================
