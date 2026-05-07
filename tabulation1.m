%% ============================================================
%  make_full_paper_tables_logchange.m
%
%  Purpose:
%      Generate full LaTeX tables for the VA-REN paper.
%
%  This script loads:
%      caiso_final_logchange_dataset_2019_2023.mat
%
%  and the result files:
%      results_baseline_seasonal_naive_logchange.mat
%      results_baseline_ridge_logchange.mat
%      results_baseline_lasso_logchange.mat
%      results_baseline_static_en_logchange.mat
%      results_varen_logchange.mat
%
%  It writes complete LaTeX table code to:
%      table_1_data_summary.tex
%      table_2_main_accuracy.tex
%      table_3_regime_accuracy.tex
%      table_4_dm_tests.tex
%      table_5_model_diagnostics.tex
%      table_6_varen_mechanism.tex
%
%  Notes:
%      1. The main comparison is VA-REN versus Static EN.
%      2. Diebold-Mariano tests use squared level-demand errors.
%      3. Positive DM means the row model improves over Static EN.
%      4. Stable hours are 10:00--14:00 PT.
%      5. Ramp hours are 16:00--21:00 PT.
%% ============================================================

clear; clc;

fprintf("\nGenerating full paper-ready LaTeX tables...\n");

%% ------------------------------------------------------------
%  Required files
%% ------------------------------------------------------------

datasetFile = "caiso_final_logchange_dataset_2019_2023.mat";

resultFiles = {
    "results_baseline_seasonal_naive_logchange.mat", "Seasonal naive", "SeasonalNaive";
    "results_baseline_ridge_logchange.mat",          "Ridge",          "Ridge";
    "results_baseline_lasso_logchange.mat",          "LASSO",          "LASSO";
    "results_baseline_static_en_logchange.mat",      "Static EN",      "StaticEN";
    "results_varen_logchange.mat",                   "VA-REN",         "VAREN"
};

methodOrder = ["SeasonalNaive", "Ridge", "LASSO", "StaticEN", "VAREN"];

for i = 1:size(resultFiles, 1)
    if ~isfile(resultFiles{i, 1})
        warning("Missing result file: %s", resultFiles{i, 1});
    end
end

if ~isfile(datasetFile)
    error("Missing dataset file: %s", datasetFile);
end

%% ------------------------------------------------------------
%  Load dataset
%% ------------------------------------------------------------

D = load(datasetFile);

allData = D.allData;
modelData = D.modelData;
evaluationData = D.evaluationData;

fprintf("Loaded dataset: %s\n", datasetFile);

%% ------------------------------------------------------------
%  Load result structures
%% ------------------------------------------------------------

Results = struct();

for i = 1:size(resultFiles, 1)

    fileName = resultFiles{i, 1};
    displayName = resultFiles{i, 2};
    keyName = resultFiles{i, 3};

    if isfile(fileName)
        S = load(fileName);
        S.DisplayName = displayName;
        S.KeyName = keyName;
        Results.(keyName) = S;
        fprintf("Loaded result file: %s\n", fileName);
    else
        Results.(keyName) = struct();
        Results.(keyName).DisplayName = displayName;
        Results.(keyName).KeyName = keyName;
    end

end

%% ------------------------------------------------------------
%  Construct summary tables
%% ------------------------------------------------------------

dataSummary = buildDataSummaryTable(allData, modelData, evaluationData);

mainAccuracy = buildMainAccuracyTable(Results, methodOrder);

regimeAccuracy = buildRegimeAccuracyTable(Results, methodOrder);

dmTable = buildDMTable(Results, methodOrder);

diagnosticsTable = buildDiagnosticsTable(Results, methodOrder);

mechanismTable = buildMechanismTable(Results);

%% ------------------------------------------------------------
%  Display MATLAB tables
%% ------------------------------------------------------------

fprintf("\nTable 1: Data summary\n");
disp(dataSummary);

fprintf("\nTable 2: Main accuracy\n");
disp(mainAccuracy);

fprintf("\nTable 3: Regime accuracy\n");
disp(regimeAccuracy);

fprintf("\nTable 4: Diebold-Mariano tests\n");
disp(dmTable);

fprintf("\nTable 5: Model diagnostics\n");
disp(diagnosticsTable);

fprintf("\nTable 6: VA-REN mechanism\n");
disp(mechanismTable);

%% ------------------------------------------------------------
%  Write CSV versions
%% ------------------------------------------------------------

writetable(dataSummary,      "table_1_data_summary.csv");
writetable(mainAccuracy,     "table_2_main_accuracy.csv");
writetable(regimeAccuracy,   "table_3_regime_accuracy.csv");
writetable(dmTable,          "table_4_dm_tests.csv");
writetable(diagnosticsTable, "table_5_model_diagnostics.csv");
writetable(mechanismTable,   "table_6_varen_mechanism.csv");

%% ------------------------------------------------------------
%  Write LaTeX tables
%% ------------------------------------------------------------

writeDataSummaryLatex(dataSummary, "table_1_data_summary.tex");

writeMainAccuracyLatex(mainAccuracy, "table_2_main_accuracy.tex");

writeRegimeAccuracyLatex(regimeAccuracy, "table_3_regime_accuracy.tex");

writeDMLatex(dmTable, "table_4_dm_tests.tex");

writeDiagnosticsLatex(diagnosticsTable, "table_5_model_diagnostics.tex");

writeMechanismLatex(mechanismTable, "table_6_varen_mechanism.tex");

%% ------------------------------------------------------------
%  Completion message
%% ------------------------------------------------------------

fprintf("\nDone. Created LaTeX files:\n");
fprintf("  table_1_data_summary.tex\n");
fprintf("  table_2_main_accuracy.tex\n");
fprintf("  table_3_regime_accuracy.tex\n");
fprintf("  table_4_dm_tests.tex\n");
fprintf("  table_5_model_diagnostics.tex\n");
fprintf("  table_6_varen_mechanism.tex\n");

fprintf("\nYou can include them in your paper with:\n");
fprintf("  \\input{table_1_data_summary.tex}\n");
fprintf("  \\input{table_2_main_accuracy.tex}\n");
fprintf("  \\input{table_3_regime_accuracy.tex}\n");
fprintf("  \\input{table_4_dm_tests.tex}\n");
fprintf("  \\input{table_5_model_diagnostics.tex}\n");
fprintf("  \\input{table_6_varen_mechanism.tex}\n");

%% ============================================================
%  Local functions
%% ============================================================

function T = buildDataSummaryTable(allData, modelData, evaluationData)

    Period = strings(3, 1);
    StartTime = strings(3, 1);
    EndTime = strings(3, 1);
    Observations = NaN(3, 1);
    MeanLoadMW = NaN(3, 1);
    StdLoadMW = NaN(3, 1);
    MeanLogChange = NaN(3, 1);
    StdLogChange = NaN(3, 1);

    datasets = {allData, modelData, evaluationData};
    names = ["Full raw sample", "Model-ready sample", "Evaluation sample"];

    for i = 1:3

        X = datasets{i};

        Period(i) = names(i);
        StartTime(i) = string(X.Timestamp(1));
        EndTime(i) = string(X.Timestamp(end));
        Observations(i) = height(X);

        if ismember("CAISO_Load_MW", X.Properties.VariableNames)
            MeanLoadMW(i) = mean(X.CAISO_Load_MW, "omitnan");
            StdLoadMW(i) = std(X.CAISO_Load_MW, "omitnan");
        end

        if ismember("Delta_Log_Load", X.Properties.VariableNames)
            MeanLogChange(i) = mean(X.Delta_Log_Load, "omitnan");
            StdLogChange(i) = std(X.Delta_Log_Load, "omitnan");
        end

    end

    T = table(Period, StartTime, EndTime, Observations, ...
        MeanLoadMW, StdLoadMW, MeanLogChange, StdLogChange);

end

function T = buildMainAccuracyTable(Results, methodOrder)

    n = numel(methodOrder);

    Method = strings(n, 1);
    ValidObs = NaN(n, 1);

    MSFE_Level = NaN(n, 1);
    RMSE_Level = NaN(n, 1);
    MAE_Level = NaN(n, 1);
    MAPE_Level = NaN(n, 1);

    RMSE_Improvement_vs_StaticEN = NaN(n, 1);
    MAE_Improvement_vs_StaticEN = NaN(n, 1);

    for i = 1:n

        key = methodOrder(i);
        S = Results.(key);

        Method(i) = S.DisplayName;

        ValidObs(i) = getValidObs(S);

        MSFE_Level(i) = getFieldOrNaN(S, "MSFE_level");
        RMSE_Level(i) = getFieldOrNaN(S, "RMSE_level");
        MAE_Level(i) = getFieldOrNaN(S, "MAE_level");
        MAPE_Level(i) = getFieldOrNaN(S, "MAPE_level");

    end

    staticRMSE = RMSE_Level(methodOrder == "StaticEN");
    staticMAE  = MAE_Level(methodOrder == "StaticEN");

    for i = 1:n
        RMSE_Improvement_vs_StaticEN(i) = 100 * (staticRMSE - RMSE_Level(i)) / staticRMSE;
        MAE_Improvement_vs_StaticEN(i)  = 100 * (staticMAE - MAE_Level(i)) / staticMAE;
    end

    T = table(Method, ValidObs, MSFE_Level, RMSE_Level, MAE_Level, ...
        MAPE_Level, RMSE_Improvement_vs_StaticEN, MAE_Improvement_vs_StaticEN);

end

function T = buildRegimeAccuracyTable(Results, methodOrder)

    n = numel(methodOrder);

    Method = strings(n, 1);

    Stable_N = NaN(n, 1);
    Stable_RMSE = NaN(n, 1);
    Stable_MAE = NaN(n, 1);
    Stable_MAPE = NaN(n, 1);

    Ramp_N = NaN(n, 1);
    Ramp_RMSE = NaN(n, 1);
    Ramp_MAE = NaN(n, 1);
    Ramp_MAPE = NaN(n, 1);

    Ramp_RMSE_Improvement_vs_StaticEN = NaN(n, 1);
    Ramp_MAE_Improvement_vs_StaticEN = NaN(n, 1);

    for i = 1:n

        key = methodOrder(i);
        S = Results.(key);
        Method(i) = S.DisplayName;

        if ~hasFields(S, ["timestamps_valid", "error_level", "y_true_level_valid"])
            continue;
        end

        h = hour(S.timestamps_valid);
        e = S.error_level(:);
        y = S.y_true_level_valid(:);

        stableIdx = h >= 10 & h <= 14;
        rampIdx = h >= 16 & h <= 21;

        Stable_N(i) = sum(stableIdx);
        Ramp_N(i) = sum(rampIdx);

        Stable_RMSE(i) = sqrt(mean(e(stableIdx).^2, "omitnan"));
        Stable_MAE(i)  = mean(abs(e(stableIdx)), "omitnan");
        Stable_MAPE(i) = mean(abs(e(stableIdx) ./ y(stableIdx)), "omitnan") * 100;

        Ramp_RMSE(i) = sqrt(mean(e(rampIdx).^2, "omitnan"));
        Ramp_MAE(i)  = mean(abs(e(rampIdx)), "omitnan");
        Ramp_MAPE(i) = mean(abs(e(rampIdx) ./ y(rampIdx)), "omitnan") * 100;

    end

    staticRampRMSE = Ramp_RMSE(methodOrder == "StaticEN");
    staticRampMAE  = Ramp_MAE(methodOrder == "StaticEN");

    for i = 1:n
        Ramp_RMSE_Improvement_vs_StaticEN(i) = 100 * (staticRampRMSE - Ramp_RMSE(i)) / staticRampRMSE;
        Ramp_MAE_Improvement_vs_StaticEN(i)  = 100 * (staticRampMAE - Ramp_MAE(i)) / staticRampMAE;
    end

    T = table(Method, ...
        Stable_N, Stable_RMSE, Stable_MAE, Stable_MAPE, ...
        Ramp_N, Ramp_RMSE, Ramp_MAE, Ramp_MAPE, ...
        Ramp_RMSE_Improvement_vs_StaticEN, Ramp_MAE_Improvement_vs_StaticEN);

end

function T = buildDMTable(Results, methodOrder)

    baselineKey = "StaticEN";

    Method = strings(0, 1);
    Sample = strings(0, 1);
    N = [];
    MeanLossDiff = [];
    DM = [];
    PValue = [];
    Significance = strings(0, 1);

    base = Results.(baselineKey);

    if ~hasFields(base, ["timestamps_valid", "error_level"])
        warning("Static EN result does not contain required fields for DM tests.");
        T = table(Method, Sample, N, MeanLossDiff, DM, PValue, Significance);
        return;
    end

    for i = 1:numel(methodOrder)

        key = methodOrder(i);

        if key == baselineKey
            continue;
        end

        S = Results.(key);

        if ~hasFields(S, ["timestamps_valid", "error_level"])
            continue;
        end

        samples = ["Full sample", "Stable hours", "Ramp hours"];

        for q = 1:numel(samples)

            sampleName = samples(q);

            [eBase, eModel, tCommon] = alignErrorsByTimestamp(base, S);

            if isempty(eBase)
                continue;
            end

            h = hour(tCommon);

            if sampleName == "Stable hours"
                idx = h >= 10 & h <= 14;
            elseif sampleName == "Ramp hours"
                idx = h >= 16 & h <= 21;
            else
                idx = true(size(h));
            end

            eBaseSub = eBase(idx);
            eModelSub = eModel(idx);

            if numel(eBaseSub) < 30
                continue;
            end

            [dmStat, pVal, meanD] = dmTestOneSided(eBaseSub, eModelSub);

            Method(end + 1, 1) = S.DisplayName;
            Sample(end + 1, 1) = sampleName;
            N(end + 1, 1) = numel(eBaseSub);
            MeanLossDiff(end + 1, 1) = meanD;
            DM(end + 1, 1) = dmStat;
            PValue(end + 1, 1) = pVal;
            Significance(end + 1, 1) = significanceStars(pVal);

        end

    end

    T = table(Method, Sample, N, MeanLossDiff, DM, PValue, Significance);

end

function T = buildDiagnosticsTable(Results, methodOrder)

    n = numel(methodOrder);

    Method = strings(n, 1);
    Alpha = NaN(n, 1);
    MeanLambda = NaN(n, 1);
    MedianLambda = NaN(n, 1);
    MeanActiveSet = NaN(n, 1);
    MedianActiveSet = NaN(n, 1);

    for i = 1:n

        key = methodOrder(i);
        S = Results.(key);

        Method(i) = S.DisplayName;

        if key == "Ridge"
            Alpha(i) = 0;
        elseif key == "LASSO"
            Alpha(i) = 1;
        elseif key == "StaticEN"
            Alpha(i) = getFieldOrNaN(S, "alpha_static");
        elseif key == "VAREN"
            Alpha(i) = getFieldOrNaN(S, "mean_alpha");
        else
            Alpha(i) = NaN;
        end

        if isfield(S, "lambda_selected")
            MeanLambda(i) = mean(S.lambda_selected, "omitnan");
            MedianLambda(i) = median(S.lambda_selected, "omitnan");
        elseif isfield(S, "lambda_valid")
            MeanLambda(i) = mean(S.lambda_valid, "omitnan");
            MedianLambda(i) = median(S.lambda_valid, "omitnan");
        end

        if isfield(S, "active_set_size")
            MeanActiveSet(i) = mean(S.active_set_size, "omitnan");
            MedianActiveSet(i) = median(S.active_set_size, "omitnan");
        elseif isfield(S, "active_set_size_valid")
            MeanActiveSet(i) = mean(S.active_set_size_valid, "omitnan");
            MedianActiveSet(i) = median(S.active_set_size_valid, "omitnan");
        end

    end

    T = table(Method, Alpha, MeanLambda, MedianLambda, MeanActiveSet, MedianActiveSet);

end

function T = buildMechanismTable(Results)

    S = Results.VAREN;

    Quantity = [
        "Mean normalized volatility, stable hours";
        "Mean normalized volatility, ramp hours";
        "Mean alpha, stable hours";
        "Mean alpha, ramp hours";
        "Mean active-set size, stable hours";
        "Mean active-set size, ramp hours";
        "Correlation between alpha and active-set size"
    ];

    Value = [
        getFieldOrNaN(S, "mean_vol_stable");
        getFieldOrNaN(S, "mean_vol_ramp");
        getFieldOrNaN(S, "mean_alpha_stable");
        getFieldOrNaN(S, "mean_alpha_ramp");
        getFieldOrNaN(S, "mean_active_stable");
        getFieldOrNaN(S, "mean_active_ramp");
        getFieldOrNaN(S, "corr_alpha_active")
    ];

    Interpretation = [
        "Volatility should be low in the solar-plateau regime";
        "Volatility should be high in the evening ramp regime";
        "Stable hours should move VA-REN toward LASSO";
        "Ramp hours should move VA-REN toward Ridge";
        "Stable hours should produce a relatively sparse model";
        "Ramp hours should retain more correlated predictors";
        "Shows whether adaptive alpha changes model sparsity"
    ];

    T = table(Quantity, Value, Interpretation);

end

%% ============================================================
%  LaTeX writer functions
%% ============================================================

function writeDataSummaryLatex(T, outFile)

    fid = fopen(outFile, "w");

    fprintf(fid, "%% Auto-generated by make_full_paper_tables_logchange.m\n");
    fprintf(fid, "\\begin{table}[t]\n");
    fprintf(fid, "\\centering\n");
    fprintf(fid, "\\caption{Summary of the CAISO hourly demand dataset and transformed modeling samples. The response variable is the hourly log-change in demand, $z_t = 100[\\log(y_t)-\\log(y_{t-1})]$.}\n");
    fprintf(fid, "\\label{tab:data_summary}\n");
    fprintf(fid, "\\begin{tabular}{@{}lrrrrr@{}}\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Sample & Observations & Mean load & Std. load & Mean $z_t$ & Std. $z_t$ \\\\\n");
    fprintf(fid, "\\midrule\n");

    for i = 1:height(T)
        fprintf(fid, "%s & %s & %s & %s & %s & %s \\\\\n", ...
            latexEscape(T.Period(i)), ...
            fmtInt(T.Observations(i)), ...
            fmtNum(T.MeanLoadMW(i), 2), ...
            fmtNum(T.StdLoadMW(i), 2), ...
            fmtNum(T.MeanLogChange(i), 4), ...
            fmtNum(T.StdLogChange(i), 4));
    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{tabular}\n");
    fprintf(fid, "\\end{table}\n");

    fclose(fid);

end

function writeMainAccuracyLatex(T, outFile)

    fid = fopen(outFile, "w");

    bestRMSE = min(T.RMSE_Level, [], "omitnan");
    bestMAE  = min(T.MAE_Level, [], "omitnan");
    bestMAPE = min(T.MAPE_Level, [], "omitnan");

    fprintf(fid, "%% Auto-generated by make_full_paper_tables_logchange.m\n");
    fprintf(fid, "\\begin{table}[t]\n");
    fprintf(fid, "\\centering\n");
    fprintf(fid, "\\caption{Full-sample forecast accuracy over the 2020--2023 evaluation period. Lower values indicate better performance. The final two columns report percentage improvement relative to the static elastic net baseline; positive values indicate improvement.}\n");
    fprintf(fid, "\\label{tab:main_accuracy}\n");
    fprintf(fid, "\\begin{tabular}{@{}lrrrrrr@{}}\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Model & $N$ & MSFE & RMSE & MAE & MAPE & $\\Delta$RMSE vs. Static EN \\\\\n");
    fprintf(fid, "\\midrule\n");

    for i = 1:height(T)

        method = formatMethod(T.Method(i));

        rmseText = fmtNum(T.RMSE_Level(i), 2);
        maeText  = fmtNum(T.MAE_Level(i), 2);
        mapeText = fmtNum(T.MAPE_Level(i), 3);

        if abs(T.RMSE_Level(i) - bestRMSE) < 1e-10
            rmseText = "\\textbf{" + rmseText + "}";
        end

        if abs(T.MAE_Level(i) - bestMAE) < 1e-10
            maeText = "\\textbf{" + maeText + "}";
        end

        if abs(T.MAPE_Level(i) - bestMAPE) < 1e-10
            mapeText = "\\textbf{" + mapeText + "}";
        end

        fprintf(fid, "%s & %s & %s & %s & %s & %s & %s \\\\\n", ...
            method, ...
            fmtInt(T.ValidObs(i)), ...
            fmtNum(T.MSFE_Level(i), 2), ...
            rmseText, ...
            maeText, ...
            mapeText, ...
            fmtSigned(T.RMSE_Improvement_vs_StaticEN(i), 2));

    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{tabular}\n");
    fprintf(fid, "\\begin{tablenotes}\n");
    fprintf(fid, "\\footnotesize\n");
    fprintf(fid, "\\item Notes: MSFE is measured in MW$^2$; RMSE and MAE are measured in MW; MAPE is measured in percent. Forecasts are generated for the standardized log-change response and then transformed back to demand levels.\n");
    fprintf(fid, "\\end{tablenotes}\n");
    fprintf(fid, "\\end{table}\n");

    fclose(fid);

end

function writeRegimeAccuracyLatex(T, outFile)

    fid = fopen(outFile, "w");

    bestStableRMSE = min(T.Stable_RMSE, [], "omitnan");
    bestRampRMSE   = min(T.Ramp_RMSE, [], "omitnan");

    fprintf(fid, "%% Auto-generated by make_full_paper_tables_logchange.m\n");
    fprintf(fid, "\\begin{table}[t]\n");
    fprintf(fid, "\\centering\n");
    fprintf(fid, "\\caption{Regime-specific forecast accuracy. Stable hours correspond to the solar-plateau window, 10:00--14:00 PT. Ramp hours correspond to the evening Duck Curve ramp, 16:00--21:00 PT. Lower values indicate better performance.}\n");
    fprintf(fid, "\\label{tab:regime_accuracy}\n");
    fprintf(fid, "\\begin{tabular}{@{}lrrrrrrr@{}}\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "& \\multicolumn{3}{c}{Stable hours} & \\multicolumn{4}{c}{Ramp hours} \\\\\n");
    fprintf(fid, "\\cmidrule(lr){2-4}\\cmidrule(lr){5-8}\n");
    fprintf(fid, "Model & RMSE & MAE & MAPE & RMSE & MAE & MAPE & $\\Delta$RMSE vs. Static EN \\\\\n");
    fprintf(fid, "\\midrule\n");

    for i = 1:height(T)

        method = formatMethod(T.Method(i));

        stableRMSEText = fmtNum(T.Stable_RMSE(i), 2);
        rampRMSEText   = fmtNum(T.Ramp_RMSE(i), 2);

        if abs(T.Stable_RMSE(i) - bestStableRMSE) < 1e-10
            stableRMSEText = "\\textbf{" + stableRMSEText + "}";
        end

        if abs(T.Ramp_RMSE(i) - bestRampRMSE) < 1e-10
            rampRMSEText = "\\textbf{" + rampRMSEText + "}";
        end

        fprintf(fid, "%s & %s & %s & %s & %s & %s & %s & %s \\\\\n", ...
            method, ...
            stableRMSEText, ...
            fmtNum(T.Stable_MAE(i), 2), ...
            fmtNum(T.Stable_MAPE(i), 3), ...
            rampRMSEText, ...
            fmtNum(T.Ramp_MAE(i), 2), ...
            fmtNum(T.Ramp_MAPE(i), 3), ...
            fmtSigned(T.Ramp_RMSE_Improvement_vs_StaticEN(i), 2));

    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{tabular}\n");
    fprintf(fid, "\\begin{tablenotes}\n");
    fprintf(fid, "\\footnotesize\n");
    fprintf(fid, "\\item Notes: RMSE and MAE are measured in MW; MAPE and $\\Delta$RMSE are measured in percent. Positive improvement values indicate lower RMSE than Static EN.\n");
    fprintf(fid, "\\end{tablenotes}\n");
    fprintf(fid, "\\end{table}\n");

    fclose(fid);

end

function writeDMLatex(T, outFile)

    fid = fopen(outFile, "w");

    fprintf(fid, "%% Auto-generated by make_full_paper_tables_logchange.m\n");
    fprintf(fid, "\\begin{table}[t]\n");
    fprintf(fid, "\\centering\n");
    fprintf(fid, "\\caption{One-sided Diebold--Mariano tests against the static elastic net baseline. The loss differential is $d_t=e^2_{\\mathrm{Static\\ EN},t}-e^2_{m,t}$, so a positive statistic indicates that the row model improves upon Static EN.}\n");
    fprintf(fid, "\\label{tab:dm_tests}\n");
    fprintf(fid, "\\begin{tabular}{@{}llrrrr@{}}\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Model & Sample & $N$ & $\\bar d$ & DM & $p$-value \\\\\n");
    fprintf(fid, "\\midrule\n");

    for i = 1:height(T)

        method = formatMethod(T.Method(i));

        fprintf(fid, "%s & %s & %s & %s & %s%s & %s \\\\\n", ...
            method, ...
            latexEscape(T.Sample(i)), ...
            fmtInt(T.N(i)), ...
            fmtNum(T.MeanLossDiff(i), 2), ...
            fmtNum(T.DM(i), 3), ...
            T.Significance(i), ...
            fmtNum(T.PValue(i), 4));

    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{tabular}\n");
    fprintf(fid, "\\begin{tablenotes}\n");
    fprintf(fid, "\\footnotesize\n");
    fprintf(fid, "\\item Notes: $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$. The one-sided alternative is that the row model has lower expected squared forecast error than Static EN.\n");
    fprintf(fid, "\\end{tablenotes}\n");
    fprintf(fid, "\\end{table}\n");

    fclose(fid);

end

function writeDiagnosticsLatex(T, outFile)

    fid = fopen(outFile, "w");

    fprintf(fid, "%% Auto-generated by make_full_paper_tables_logchange.m\n");
    fprintf(fid, "\\begin{table}[t]\n");
    fprintf(fid, "\\centering\n");
    fprintf(fid, "\\caption{Model diagnostics for penalized regression methods. For VA-REN, the reported $\\alpha$ is the sample mean of the adaptive mixing parameter. For Ridge and LASSO, $\\alpha$ is fixed by definition.}\n");
    fprintf(fid, "\\label{tab:model_diagnostics}\n");
    fprintf(fid, "\\begin{tabular}{@{}lrrrrr@{}}\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Model & $\\alpha$ & Mean $\\lambda$ & Median $\\lambda$ & Mean active set & Median active set \\\\\n");
    fprintf(fid, "\\midrule\n");

    for i = 1:height(T)

        method = formatMethod(T.Method(i));

        fprintf(fid, "%s & %s & %s & %s & %s & %s \\\\\n", ...
            method, ...
            fmtNum(T.Alpha(i), 4), ...
            fmtNum(T.MeanLambda(i), 6), ...
            fmtNum(T.MedianLambda(i), 6), ...
            fmtNum(T.MeanActiveSet(i), 2), ...
            fmtNum(T.MedianActiveSet(i), 2));

    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{tabular}\n");
    fprintf(fid, "\\begin{tablenotes}\n");
    fprintf(fid, "\\footnotesize\n");
    fprintf(fid, "\\item Notes: Active-set size is the number of nonzero estimated coefficients in the rolling-window fitted model. Seasonal naive has no penalized-regression diagnostics.\n");
    fprintf(fid, "\\end{tablenotes}\n");
    fprintf(fid, "\\end{table}\n");

    fclose(fid);

end

function writeMechanismLatex(T, outFile)

    fid = fopen(outFile, "w");

    fprintf(fid, "%% Auto-generated by make_full_paper_tables_logchange.m\n");
    fprintf(fid, "\\begin{table}[t]\n");
    fprintf(fid, "\\centering\n");
    fprintf(fid, "\\caption{VA-REN mechanism diagnostics. Stable hours are 10:00--14:00 PT; ramp hours are 16:00--21:00 PT. The expected pattern is higher volatility, lower $\\alpha_t$, and a denser active set during ramp hours.}\n");
    fprintf(fid, "\\label{tab:varen_mechanism}\n");
    fprintf(fid, "\\begin{tabular}{@{}lrp{6.4cm}@{}}\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Quantity & Value & Interpretation \\\\\n");
    fprintf(fid, "\\midrule\n");

    for i = 1:height(T)

        fprintf(fid, "%s & %s & %s \\\\\n", ...
            latexEscape(T.Quantity(i)), ...
            fmtNum(T.Value(i), 4), ...
            latexEscape(T.Interpretation(i)));

    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{tabular}\n");
    fprintf(fid, "\\end{table}\n");

    fclose(fid);

end

%% ============================================================
%  Statistical helper functions
%% ============================================================

function [dmStat, pVal, meanD] = dmTestOneSided(eBase, eModel)

    eBase = eBase(:);
    eModel = eModel(:);

    valid = ~isnan(eBase) & ~isnan(eModel);

    eBase = eBase(valid);
    eModel = eModel(valid);

    d = eBase.^2 - eModel.^2;

    T = numel(d);

    meanD = mean(d, "omitnan");

    dCentered = d - meanD;

    h = floor(4 * (T / 100)^(2 / 9));

    gamma0 = mean(dCentered .* dCentered, "omitnan");

    longRunVar = gamma0;

    for j = 1:h
        gammaj = mean(dCentered((j + 1):end) .* dCentered(1:(end - j)), "omitnan");
        weight = 1 - j / (h + 1);
        longRunVar = longRunVar + 2 * weight * gammaj;
    end

    if longRunVar <= 0 || isnan(longRunVar)
        dmStat = NaN;
        pVal = NaN;
        return;
    end

    dmStat = meanD / sqrt(longRunVar / T);

    % One-sided p-value: P(Z > DM)
    pVal = 0.5 * erfc(dmStat / sqrt(2));

end

function [eBase, eModel, tCommon] = alignErrorsByTimestamp(base, model)

    tBase = base.timestamps_valid(:);
    tModel = model.timestamps_valid(:);

    eBaseAll = base.error_level(:);
    eModelAll = model.error_level(:);

    [tCommon, ia, ib] = intersect(tBase, tModel);

    eBase = eBaseAll(ia);
    eModel = eModelAll(ib);

end

function stars = significanceStars(p)

    if isnan(p)
        stars = "";
    elseif p < 0.01
        stars = "$^{***}$";
    elseif p < 0.05
        stars = "$^{**}$";
    elseif p < 0.10
        stars = "$^{*}$";
    else
        stars = "";
    end

end

%% ============================================================
%  Formatting helper functions
%% ============================================================

function value = getFieldOrNaN(S, fieldName)

    if isfield(S, fieldName)
        value = S.(fieldName);
        if numel(value) > 1
            value = value(1);
        end
    else
        value = NaN;
    end

end

function tf = hasFields(S, fields)

    tf = true;

    for i = 1:numel(fields)
        if ~isfield(S, fields(i))
            tf = false;
            return;
        end
    end

end

function n = getValidObs(S)

    if isfield(S, "error_level")
        n = numel(S.error_level);
    elseif isfield(S, "y_true_level_valid")
        n = numel(S.y_true_level_valid);
    else
        n = NaN;
    end

end

function s = fmtNum(x, decimals)

    if isempty(x) || isnan(x)
        s = "---";
    else
        s = sprintf("%.*f", decimals, x);
    end

end

function s = fmtSigned(x, decimals)

    if isempty(x) || isnan(x)
        s = "---";
    else
        if x > 0
            s = sprintf("+%.*f", decimals, x);
        else
            s = sprintf("%.*f", decimals, x);
        end
    end

end

function s = fmtInt(x)

    if isempty(x) || isnan(x)
        s = "---";
    else
        s = sprintf("%d", round(x));
    end

end

function out = latexEscape(in)

    out = string(in);

    out = replace(out, "_", "\\_");
    out = replace(out, "%", "\\%");
    out = replace(out, "&", "\\&");

end

function out = formatMethod(method)

    method = string(method);

    if method == "VA-REN"
        out = "\\textbf{VA-REN}";
    else
        out = latexEscape(method);
    end

end