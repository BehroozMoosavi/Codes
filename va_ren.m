%% ============================================================
%  va_ren.m
%  FILE 06 of 30
%
%  Volatility-Adaptive Rolling Elastic Net  (VA-REN)
%  alpha_t = 1 - NormalizedVolatility_t
%
%  High vol (ramp, 16-21): alpha_t small => Ridge-like (dense)
%  Low  vol (stable,10-14): alpha_t ~1   => LASSO-like (sparse)
%
%  Also stores selection_matrix_valid [T x p] for
%  regime_selection.m analysis.
%
%  Input : caiso_final_logchange_dataset_2019_2023.mat
%  Output: results_varen_logchange.mat
%% ============================================================

clear; clc;
load("caiso_final_logchange_dataset_2019_2023.mat");
fprintf("Running VA-REN...\n");

W         = 720;
K         = 5;
numLambda = 60;

X_all        = X;
y_all_std    = y_eval_z_std(:);
y_all_raw    = y_eval_z_raw(:);
y_true_level = y_eval_level(:);
y_lag1_level = y_lag1_eval_level(:);
alpha_all    = evaluationData.Alpha_VAREN(:);
vol_all      = evaluationData.NormalizedVolatility(:);
T            = length(y_all_std);
p            = size(X_all,2);

fprintf("T=%d  p=%d\n", T, p);

yhat_z_std          = NaN(T,1);
yhat_z_raw          = NaN(T,1);
yhat_level          = NaN(T,1);
lambda_selected     = NaN(T,1);
alpha_selected      = NaN(T,1);
vol_selected        = NaN(T,1);
active_set_size     = NaN(T,1);
selection_matrix    = false(T,p);

for t = W+1:T

    if mod(t,500)==0, fprintf("  VA-REN step %d/%d\n",t,T); end

    trainIdx = (t-W):(t-1);
    X_train  = X_all(trainIdx,:);
    y_train  = y_all_std(trainIdx);
    X_test   = X_all(t,:);
    alpha_t  = max(0.001, min(1.000, alpha_all(t)));
    vol_t    = vol_all(t);

    validTrain = all(~isnan(X_train),2) & ~isnan(y_train);
    X_train    = X_train(validTrain,:);
    y_train    = y_train(validTrain);

    if size(X_train,1) < K+10 || any(isnan(X_test)) || isnan(alpha_t)
        continue;
    end

    try
        [B,FI] = lasso(X_train,y_train,"Alpha",alpha_t,"CV",K, ...
                       "NumLambda",numLambda,"Standardize",false);
        idx   = FI.IndexMinMSE;
        beta  = B(:,idx);
        inter = FI.Intercept(idx);

        yhat_z_std(t) = inter + X_test*beta;
        yhat_z_raw(t) = z_mean + z_std*yhat_z_std(t);
        yhat_level(t) = y_lag1_level(t)*exp(yhat_z_raw(t)/100);

        lambda_selected(t)  = FI.Lambda(idx);
        alpha_selected(t)   = alpha_t;
        vol_selected(t)     = vol_t;
        active_set_size(t)  = sum(beta~=0);
        selection_matrix(t,:) = (beta~=0)';
    catch ME
        warning("VA-REN t=%d: %s",t,ME.message);
    end

end

valid = ~isnan(yhat_z_std) & ~isnan(yhat_level) & ~isnan(y_all_std) & ...
        ~isnan(y_true_level) & y_true_level>0 & yhat_level>0;

y_true_z_std_valid    = y_all_std(valid);
y_true_z_raw_valid    = y_all_raw(valid);
y_true_level_valid    = y_true_level(valid);
yhat_z_std_valid      = yhat_z_std(valid);
yhat_z_raw_valid      = yhat_z_raw(valid);
yhat_level_valid      = yhat_level(valid);
timestamps_valid      = timestamps_eval(valid);
alpha_valid           = alpha_selected(valid);
volatility_valid      = vol_selected(valid);
lambda_valid          = lambda_selected(valid);
active_set_size_valid = active_set_size(valid);
selection_matrix_valid= selection_matrix(valid,:);

error_z_std = y_true_z_std_valid - yhat_z_std_valid;
error_z_raw = y_true_z_raw_valid - yhat_z_raw_valid;
error_level = y_true_level_valid - yhat_level_valid;

MSFE_level = mean(error_level.^2,"omitnan"); RMSE_level = sqrt(MSFE_level);
MAE_level  = mean(abs(error_level),"omitnan");
MAPE_level = mean(abs(error_level./y_true_level_valid),"omitnan")*100;
MSFE_z_raw = mean(error_z_raw.^2,"omitnan"); RMSE_z_raw = sqrt(MSFE_z_raw);
MAE_z_raw  = mean(abs(error_z_raw),"omitnan");
MSFE_z_std = mean(error_z_std.^2,"omitnan"); RMSE_z_std = sqrt(MSFE_z_std);
MAE_z_std  = mean(abs(error_z_std),"omitnan");

h_v       = hour(timestamps_valid);
stableIdx = h_v>=10 & h_v<=14;
rampIdx   = h_v>=16 & h_v<=21;

RMSE_level_stable = sqrt(mean(error_level(stableIdx).^2,"omitnan"));
MAE_level_stable  = mean(abs(error_level(stableIdx)),"omitnan");
MAPE_level_stable = mean(abs(error_level(stableIdx)./y_true_level_valid(stableIdx)),"omitnan")*100;
RMSE_level_ramp   = sqrt(mean(error_level(rampIdx).^2,"omitnan"));
MAE_level_ramp    = mean(abs(error_level(rampIdx)),"omitnan");
MAPE_level_ramp   = mean(abs(error_level(rampIdx)./y_true_level_valid(rampIdx)),"omitnan")*100;

mean_alpha        = mean(alpha_valid,"omitnan");
median_alpha      = median(alpha_valid,"omitnan");
mean_alpha_stable = mean(alpha_valid(stableIdx),"omitnan");
mean_alpha_ramp   = mean(alpha_valid(rampIdx),"omitnan");
mean_vol_stable   = mean(volatility_valid(stableIdx),"omitnan");
mean_vol_ramp     = mean(volatility_valid(rampIdx),"omitnan");
mean_active_set   = mean(active_set_size_valid,"omitnan");
median_active_set = median(active_set_size_valid,"omitnan");
mean_active_stable= mean(active_set_size_valid(stableIdx),"omitnan");
mean_active_ramp  = mean(active_set_size_valid(rampIdx),"omitnan");
corr_alpha_active = corr(alpha_valid,active_set_size_valid,"Rows","complete");

fprintf("\nVA-REN Results: N=%d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
        numel(error_level), RMSE_level, MAPE_level);
fprintf("  Stable RMSE=%.2f  Ramp RMSE=%.2f\n", RMSE_level_stable, RMSE_level_ramp);
fprintf("  Mean alpha=%.4f  stable=%.4f  ramp=%.4f\n", ...
        mean_alpha, mean_alpha_stable, mean_alpha_ramp);
fprintf("  Mean active set=%.2f  stable=%.2f  ramp=%.2f\n", ...
        mean_active_set, mean_active_stable, mean_active_ramp);
fprintf("  Corr(alpha,active set)=%.4f\n", corr_alpha_active);
fprintf("  Selection matrix: [%d x %d]\n", ...
        size(selection_matrix_valid,1), size(selection_matrix_valid,2));

save("results_varen_logchange.mat", ...
    "y_true_z_std_valid","y_true_z_raw_valid","y_true_level_valid", ...
    "yhat_z_std_valid","yhat_z_raw_valid","yhat_level_valid", ...
    "error_z_std","error_z_raw","error_level", ...
    "MSFE_z_std","RMSE_z_std","MAE_z_std", ...
    "MSFE_z_raw","RMSE_z_raw","MAE_z_raw", ...
    "MSFE_level","RMSE_level","MAE_level","MAPE_level", ...
    "RMSE_level_stable","MAE_level_stable","MAPE_level_stable", ...
    "RMSE_level_ramp","MAE_level_ramp","MAPE_level_ramp", ...
    "timestamps_valid","alpha_valid","volatility_valid", ...
    "lambda_valid","active_set_size_valid","selection_matrix_valid", ...
    "mean_alpha","median_alpha","mean_alpha_stable","mean_alpha_ramp", ...
    "mean_vol_stable","mean_vol_ramp","mean_active_set","median_active_set", ...
    "mean_active_stable","mean_active_ramp","corr_alpha_active", ...
    "stableIdx","rampIdx","W","K","numLambda");

fprintf("Saved: results_varen_logchange.mat\n");
