%% ============================================================
%  master_empirical.m
%
%  Complete empirical analysis pipeline for VA-REN paper.
%  Fixed version — correct MCS array indexing.
%
%  RUN ORDER:
%    1. data_clean.m
%    2. run_all_models_parallel.m
%    3. run_all_robustness_parallel.m  (then fix_r4_subperiod.m)
%    4. master_empirical.m   <-- this script
%
%  OUTPUTS:
%    output/tables/   — all .tex and .csv table files
%    output/figures/  — all .pdf and .png figure files
%% ============================================================

clear; clc;

fprintf("============================================================\n");
fprintf("  VA-REN  Complete Empirical Analysis\n");
fprintf("  Start: %s\n", datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf("============================================================\n\n");

%% ----------------------------------------------------------
%  0. Setup directories
%% ----------------------------------------------------------

outDir_tab = "output/tables";
outDir_fig = "output/figures";
if ~exist(outDir_tab,"dir"), mkdir(outDir_tab); end
if ~exist(outDir_fig,"dir"), mkdir(outDir_fig); end

%% ----------------------------------------------------------
%  1. Load all result files
%% ----------------------------------------------------------

dataFile = "caiso_final_logchange_dataset_2019_2023.mat";

if ~isfile(dataFile)
    error("Missing: %s\nRun data_clean.m first.", dataFile);
end

resultFiles = struct( ...
    "Naive",    "results_baseline_seasonal_naive_logchange.mat", ...
    "Ridge",    "results_baseline_ridge_logchange.mat", ...
    "LASSO",    "results_baseline_lasso_logchange.mat", ...
    "StaticEN", "results_baseline_static_en_logchange.mat", ...
    "VAREN",    "results_varen_logchange.mat");

fnames = fieldnames(resultFiles);
R = struct();

fprintf("--- Loading result files ---\n");
for i = 1:numel(fnames)
    fn = resultFiles.(fnames{i});
    if isfile(fn)
        R.(fnames{i}) = load(fn);
        S = R.(fnames{i});
        fprintf("  %-10s  N=%-6d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
            fnames{i}, numel(S.error_level), S.RMSE_level, S.MAPE_level);
    else
        error("Missing result file: %s\nRun the model script first.", fn);
    end
end

D          = load(dataFile);
vol_series = D.evaluationData.NormalizedVolatility;
ts_eval    = D.evaluationData.Timestamp;
pred_names = D.standardizedNames(:);

fprintf("\nDataset loaded: %s\n\n", dataFile);

%% ----------------------------------------------------------
%  2. Pairwise DM Tests  (Layer 1)
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  LAYER 1: Pairwise DM Tests\n");
fprintf("============================================================\n\n");

dm_results = struct();

pairs = {
    "Naive",    "VAREN",    "Seasonal Naive", "VA-REN";
    "Ridge",    "VAREN",    "Ridge",          "VA-REN";
    "LASSO",    "VAREN",    "LASSO",          "VA-REN";
    "StaticEN", "VAREN",    "Static EN",      "VA-REN";
    "Naive",    "StaticEN", "Seasonal Naive", "Static EN"
};

for ci = 1:size(pairs,1)
    kA  = pairs{ci,1};  kB  = pairs{ci,2};
    lA  = pairs{ci,3};  lB  = pairs{ci,4};
    key = kA + "_vs_" + kB;
    dm_results.(key) = dm_test(R.(kA), R.(kB), lA, lB);
end

%% ----------------------------------------------------------
%  3. Regime DM Tests  (Layer 2)
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  LAYER 2: Regime DM Tests  (VA-REN vs Static EN)\n");
fprintf("============================================================\n\n");

[dm_stable, dm_ramp, dm_full] = ...
    dm_regime(R.StaticEN, R.VAREN, "Static EN", "VA-REN");

%% ----------------------------------------------------------
%  4. Model Confidence Set  (Layer 3)
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  LAYER 3: Model Confidence Set\n");
fprintf("============================================================\n\n");

mkeys = {"Naive","Ridge","LASSO","StaticEN","VAREN"};
mlabs = ["Seasonal Naive","Ridge","LASSO","Static EN","VA-REN"];

%% Find timestamps common to ALL five models
t_common = R.Naive.timestamps_valid(:);
for i = 2:5
    t_common = intersect(t_common, R.(mkeys{i}).timestamps_valid(:));
end

fprintf("Common timestamps for MCS: %d\n\n", numel(t_common));

%% Align each model's errors to the common timestamp grid
%% Guard against ia==0 (unmatched rows) before indexing
ec = cell(5,1);

for i = 1:5

    ts_model  = R.(mkeys{i}).timestamps_valid(:);
    err_model = R.(mkeys{i}).error_level(:);

    [~, ia] = ismember(t_common, ts_model);

    valid          = ia > 0;
    aligned        = NaN(numel(t_common), 1);
    aligned(valid) = err_model(ia(valid));
    ec{i}          = aligned;

    if any(~valid)
        warning("MCS: %d unmatched timestamps for %s", sum(~valid), mkeys{i});
    end

end

[mcs_set, elim_table] = mcs_test(ec, mlabs, 0.10, 1000);

fprintf("MCS (alpha=0.10) survives: %s\n\n", strjoin(mcs_set, ", "));

%% ----------------------------------------------------------
%  5. Giacomini-White Test  (Layer 4)
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  LAYER 4: Giacomini-White Test\n");
fprintf("============================================================\n\n");

gw_result = gw_test(R.StaticEN, R.VAREN, ...
                    vol_series, ts_eval, ...
                    "Static EN", "VA-REN");

%% ----------------------------------------------------------
%  6. Regime Selection Analysis
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  REGIME SELECTION ANALYSIS\n");
fprintf("============================================================\n\n");

if isfield(R.VAREN, "selection_matrix_valid")

    sel_result = regime_selection( ...
        "results_varen_logchange.mat", pred_names, 500);

    plot_regime_selection(sel_result, outDir_fig);
    write_selection_latex(sel_result, outDir_tab);
    save(fullfile(outDir_tab,"regime_selection_result.mat"), "sel_result");

else
    fprintf("  selection_matrix_valid not found in VA-REN results.\n");
    fprintf("  Skipping regime selection analysis.\n\n");
    sel_result = [];
end

%% ----------------------------------------------------------
%  7. Load robustness results (already computed)
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  LOADING ROBUSTNESS RESULTS\n");
fprintf("============================================================\n\n");

results_W   = [];
results_Wv  = [];
results_map = [];
T_sub       = [];
results_al  = [];

if isfile("robustness_W.mat")
    tmp = load("robustness_W.mat"); results_W = tmp.results_W;
    fprintf("  [OK] robustness_W.mat\n");
else
    fprintf("  [--] robustness_W.mat not found\n");
end

if isfile("robustness_Wv.mat")
    tmp = load("robustness_Wv.mat"); results_Wv = tmp.results_Wv;
    fprintf("  [OK] robustness_Wv.mat\n");
else
    fprintf("  [--] robustness_Wv.mat not found\n");
end

if isfile("robustness_map.mat")
    tmp = load("robustness_map.mat"); results_map = tmp.results_map;
    fprintf("  [OK] robustness_map.mat\n");
else
    fprintf("  [--] robustness_map.mat not found\n");
end

if isfile("subperiod.mat")
    tmp = load("subperiod.mat"); T_sub = tmp.T_sub;
    fprintf("  [OK] subperiod.mat\n");
else
    fprintf("  [--] subperiod.mat not found\n");
end

if isfile("results_adaptive_lasso.mat")
    results_al = load("results_adaptive_lasso.mat");
    fprintf("  [OK] results_adaptive_lasso.mat\n");
else
    fprintf("  [--] results_adaptive_lasso.mat not found\n");
end

fprintf("\n");

%% ----------------------------------------------------------
%  8. Figures
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  GENERATING FIGURES\n");
fprintf("============================================================\n\n");

plot_all_results( ...
    dataFile, ...
    "results_varen_logchange.mat", ...
    "results_baseline_static_en_logchange.mat", ...
    "results_baseline_seasonal_naive_logchange.mat", ...
    results_W, outDir_fig);

%% ----------------------------------------------------------
%  9. Paper tables
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  GENERATING PAPER TABLES\n");
fprintf("============================================================\n\n");

tabulation;

%% Re-declare everything after tabulation.m runs clear at top
outDir_tab = 'output/tables';
outDir_fig = 'output/figures';
if ~exist(outDir_tab,'dir'), mkdir(outDir_tab); end
if ~exist(outDir_fig,'dir'), mkdir(outDir_fig); end

%% Copy generated table files to output/tables/
for ext = {'.tex', '.csv'}
    ff = dir(['table_*' ext{1}]);
    for k = 1:numel(ff)
        copyfile(ff(k).name, fullfile(outDir_tab, ff(k).name));
    end
end

%% Reload results needed for final summary
resultFiles2 = struct( ...
    'Naive',    'results_baseline_seasonal_naive_logchange.mat', ...
    'Ridge',    'results_baseline_ridge_logchange.mat', ...
    'LASSO',    'results_baseline_lasso_logchange.mat', ...
    'StaticEN', 'results_baseline_static_en_logchange.mat', ...
    'VAREN',    'results_varen_logchange.mat');
fnames2 = fieldnames(resultFiles2);
R = struct();
for i = 1:numel(fnames2)
    R.(fnames2{i}) = load(resultFiles2.(fnames2{i}));
end

%% Reload statistical test results
tmp = load(fullfile(outDir_tab,'dm_regime_results.mat'));
dm_full   = tmp.dm_full;
dm_stable = tmp.dm_stable;
dm_ramp   = tmp.dm_ramp;

tmp2      = load(fullfile(outDir_tab,'gw_result.mat'));
gw_result = tmp2.gw_result;

tmp3    = load(fullfile(outDir_tab,'mcs_results.mat'));
mcs_set = tmp3.mcs_set;

%% ----------------------------------------------------------
%  10. Final summary
%% ----------------------------------------------------------

fprintf("============================================================\n");
fprintf("  FINAL RESULTS SUMMARY\n");
fprintf("============================================================\n\n");

fprintf("Main model accuracy:\n");
fprintf("  %-16s  %8s  %8s  %8s\n","Model","RMSE(MW)","MAE(MW)","MAPE(%%)");
fprintf("  %-16s  %8s  %8s  %8s\n","-----","--------","-------","-------");

model_list = {"Naive","Ridge","LASSO","StaticEN","VAREN"};
disp_names = {"Seasonal Naive","Ridge","LASSO","Static EN","VA-REN"};
staticRMSE = R.StaticEN.RMSE_level;

for i = 1:5
    S    = R.(model_list{i});
    dRMSE= 100*(staticRMSE - S.RMSE_level)/staticRMSE;
    fprintf("  %-16s  %8.2f  %8.2f  %8.4f  (Delta=%+.3f%%)\n", ...
        disp_names{i}, S.RMSE_level, S.MAE_level, S.MAPE_level, dRMSE);
end

fprintf("\nDM test: VA-REN vs Static EN\n");
fprintf("  Full sample  : DM=%+.4f  p=%.4f  %s\n", ...
    dm_full.DM,   dm_full.pVal,   dm_full.stars);
fprintf("  Stable hours : DM=%+.4f  p=%.4f  %s\n", ...
    dm_stable.DM, dm_stable.pVal, dm_stable.stars);
fprintf("  Ramp hours   : DM=%+.4f  p=%.4f  %s\n", ...
    dm_ramp.DM,   dm_ramp.pVal,   dm_ramp.stars);

fprintf("\nGW test:\n");
fprintf("  gamma_1=%+.4f  se=%.4f  t=%+.4f  p=%.4f  %s\n", ...
    gw_result.gamma1, gw_result.se_gamma1, ...
    gw_result.t_stat,  gw_result.p_val, gw_result.stars);

fprintf("\nMCS (alpha=0.10): %s\n", strjoin(mcs_set,", "));

fprintf("\nVA-REN mechanism:\n");
if isfield(R.VAREN,"mean_alpha")
    fprintf("  Mean alpha            : %.4f\n", R.VAREN.mean_alpha);
    fprintf("  Mean alpha (stable)   : %.4f\n", R.VAREN.mean_alpha_stable);
    fprintf("  Mean alpha (ramp)     : %.4f\n", R.VAREN.mean_alpha_ramp);
    fprintf("  Mean vol   (stable)   : %.4f\n", R.VAREN.mean_vol_stable);
    fprintf("  Mean vol   (ramp)     : %.4f\n", R.VAREN.mean_vol_ramp);
    fprintf("  Mean active set       : %.2f\n",  R.VAREN.mean_active_set);
    fprintf("  Mean active (stable)  : %.2f\n",  R.VAREN.mean_active_stable);
    fprintf("  Mean active (ramp)    : %.2f\n",  R.VAREN.mean_active_ramp);
    fprintf("  Corr(alpha,active set): %.4f\n",  R.VAREN.corr_alpha_active);
end

if isfield(R.VAREN,"RMSE_level_stable")
    fprintf("\nRegime accuracy:\n");
    fprintf("  %-16s  %8s  %8s\n","Model","Stable RMSE","Ramp RMSE");
    fprintf("  %-16s  %8s  %8s\n","-----","-----------","---------");
    for i = 1:5
        S = R.(model_list{i});
        if isfield(S,"RMSE_level_stable")
            fprintf("  %-16s  %8.2f  %8.2f\n", ...
                disp_names{i}, S.RMSE_level_stable, S.RMSE_level_ramp);
        else
            hv = hour(S.timestamps_valid);
            e  = S.error_level;
            sr = sqrt(mean(e(hv>=10&hv<=14).^2,"omitnan"));
            rr = sqrt(mean(e(hv>=16&hv<=21).^2,"omitnan"));
            fprintf("  %-16s  %8.2f  %8.2f\n", disp_names{i}, sr, rr);
        end
    end
end

%% ----------------------------------------------------------
%  11. Save all outputs
%% ----------------------------------------------------------

fprintf("\n============================================================\n");
fprintf("  SAVING OUTPUTS\n");
fprintf("============================================================\n\n");

save(fullfile(outDir_tab,"dm_results.mat"),       "dm_results");
save(fullfile(outDir_tab,"dm_regime_results.mat"), "dm_stable","dm_ramp","dm_full");
save(fullfile(outDir_tab,"mcs_results.mat"),       "mcs_set","elim_table");
save(fullfile(outDir_tab,"gw_result.mat"),         "gw_result");

if ~isempty(T_sub)
    writetable(T_sub, fullfile(outDir_tab,"subperiod_stability.csv"));
end

fprintf("Tables  -> %s\n", outDir_tab);
fprintf("Figures -> %s\n\n", outDir_fig);

fprintf("============================================================\n");
fprintf("  Pipeline complete: %s\n", datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf("  Copy the FINAL RESULTS SUMMARY above and send it.\n");
fprintf("  We will write the paper from those numbers.\n");
fprintf("============================================================\n");