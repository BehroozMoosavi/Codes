% Load all 5 results and print N for each
files = {
    "results_baseline_seasonal_naive_logchange.mat",  "Naive";
    "results_baseline_ridge_logchange.mat",           "Ridge";
    "results_baseline_lasso_logchange.mat",           "LASSO";
    "results_baseline_static_en_logchange.mat",       "Static EN";
    "results_varen_logchange.mat",                    "VA-REN"
};

for i = 1:5
    S = load(files{i,1});
    fprintf("%-12s  N=%d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
        files{i,2}, numel(S.error_level), S.RMSE_level, S.MAPE_level);
end