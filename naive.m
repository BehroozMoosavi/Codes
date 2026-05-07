%% ============================================================
%  naive.m
%  FILE 02 of 30
%
%  Seasonal Naive Baseline
%  yhat_t = y_{t-168}  (same hour, one week ago)
%
%  Input : caiso_final_logchange_dataset_2019_2023.mat
%  Output: results_baseline_seasonal_naive_logchange.mat
%% ============================================================

clear; clc;
load("caiso_final_logchange_dataset_2019_2023.mat");
fprintf("Running Seasonal Naive baseline...\n");

y_true_level = y_eval_level(:);
y_lag1_level = y_lag1_eval_level(:);
z_true_raw   = y_eval_z_raw(:);
z_true_std   = y_eval_z_std(:);

yhat_level = evaluationData.LoadLag_168(:);
yhat_z_raw = 100 * log(yhat_level ./ y_lag1_level);
yhat_z_std = (yhat_z_raw - z_mean) ./ z_std;

valid = ~isnan(y_true_level) & ~isnan(y_lag1_level) & ...
        ~isnan(z_true_raw)   & ~isnan(z_true_std)   & ...
        ~isnan(yhat_level)   & ~isnan(yhat_z_raw)   & ...
        ~isnan(yhat_z_std)   & y_true_level > 0      & ...
        y_lag1_level > 0     & yhat_level > 0;

y_true_level_valid = y_true_level(valid);
y_lag1_level_valid = y_lag1_level(valid);
z_true_raw_valid   = z_true_raw(valid);
z_true_std_valid   = z_true_std(valid);
yhat_level_valid   = yhat_level(valid);
yhat_z_raw_valid   = yhat_z_raw(valid);
yhat_z_std_valid   = yhat_z_std(valid);
timestamps_valid   = timestamps_eval(valid);

error_level = y_true_level_valid - yhat_level_valid;
error_z_raw = z_true_raw_valid   - yhat_z_raw_valid;
error_z_std = z_true_std_valid   - yhat_z_std_valid;

MSFE_level = mean(error_level.^2,"omitnan");
RMSE_level = sqrt(MSFE_level);
MAE_level  = mean(abs(error_level),"omitnan");
MAPE_level = mean(abs(error_level./y_true_level_valid),"omitnan")*100;

MSFE_z_raw = mean(error_z_raw.^2,"omitnan");
RMSE_z_raw = sqrt(MSFE_z_raw);
MAE_z_raw  = mean(abs(error_z_raw),"omitnan");

MSFE_z_std = mean(error_z_std.^2,"omitnan");
RMSE_z_std = sqrt(MSFE_z_std);
MAE_z_std  = mean(abs(error_z_std),"omitnan");

fprintf("\nSeasonal Naive Results:\n");
fprintf("  N          = %d\n",     numel(error_level));
fprintf("  RMSE level = %.4f MW\n", RMSE_level);
fprintf("  MAE  level = %.4f MW\n", MAE_level);
fprintf("  MAPE level = %.4f %%\n", MAPE_level);

save("results_baseline_seasonal_naive_logchange.mat", ...
    "y_true_level_valid","y_lag1_level_valid", ...
    "z_true_raw_valid","z_true_std_valid", ...
    "yhat_level_valid","yhat_z_raw_valid","yhat_z_std_valid", ...
    "error_level","error_z_raw","error_z_std", ...
    "MSFE_level","RMSE_level","MAE_level","MAPE_level", ...
    "MSFE_z_raw","RMSE_z_raw","MAE_z_raw", ...
    "MSFE_z_std","RMSE_z_std","MAE_z_std", ...
    "timestamps_valid");

fprintf("Saved: results_baseline_seasonal_naive_logchange.mat\n");
