%% ============================================================
%  lasso_method.m
%  FILE 04 of 30
%
%  Rolling LASSO Baseline  (alpha=1, pure L1)
%  W=720, K=5, numLambda=60.
%
%  Input : caiso_final_logchange_dataset_2019_2023.mat
%  Output: results_baseline_lasso_logchange.mat
%% ============================================================

clear; clc;
load("caiso_final_logchange_dataset_2019_2023.mat");
fprintf("Running Rolling LASSO baseline...\n");

W         = 720;
K         = 5;
numLambda = 60;

X_all        = X;
y_all_std    = y_eval_z_std(:);
y_all_raw    = y_eval_z_raw(:);
y_true_level = y_eval_level(:);
y_lag1_level = y_lag1_eval_level(:);
T            = length(y_all_std);
p            = size(X_all,2);

fprintf("T=%d  p=%d\n", T, p);

yhat_z_std      = NaN(T,1);
yhat_z_raw      = NaN(T,1);
yhat_level      = NaN(T,1);
lambda_selected = NaN(T,1);
active_set_size = NaN(T,1);

for t = W+1:T

    if mod(t,500)==0, fprintf("  LASSO step %d/%d\n",t,T); end

    trainIdx   = (t-W):(t-1);
    X_train    = X_all(trainIdx,:);
    y_train    = y_all_std(trainIdx);
    X_test     = X_all(t,:);

    validTrain = all(~isnan(X_train),2) & ~isnan(y_train);
    X_train    = X_train(validTrain,:);
    y_train    = y_train(validTrain);

    if size(X_train,1) < K+10 || any(isnan(X_test)), continue; end

    try
        [B,FI] = lasso(X_train,y_train,"Alpha",1,"CV",K, ...
                       "NumLambda",numLambda,"Standardize",false);
        idx   = FI.IndexMinMSE;
        beta  = B(:,idx);
        inter = FI.Intercept(idx);

        yhat_z_std(t) = inter + X_test*beta;
        yhat_z_raw(t) = z_mean + z_std*yhat_z_std(t);
        yhat_level(t) = y_lag1_level(t)*exp(yhat_z_raw(t)/100);

        lambda_selected(t) = FI.Lambda(idx);
        active_set_size(t) = sum(beta~=0);
    catch ME
        warning("LASSO t=%d: %s",t,ME.message);
    end

end

valid = ~isnan(yhat_z_std) & ~isnan(yhat_level) & ~isnan(y_all_std) & ...
        ~isnan(y_true_level) & y_true_level>0 & yhat_level>0;

y_true_z_std_valid = y_all_std(valid);
y_true_z_raw_valid = y_all_raw(valid);
y_true_level_valid = y_true_level(valid);
yhat_z_std_valid   = yhat_z_std(valid);
yhat_z_raw_valid   = yhat_z_raw(valid);
yhat_level_valid   = yhat_level(valid);
timestamps_valid   = timestamps_eval(valid);

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

fprintf("\nLASSO Results: N=%d  RMSE=%.2f MW  MAE=%.2f MW  MAPE=%.4f%%\n", ...
        numel(error_level), RMSE_level, MAE_level, MAPE_level);
fprintf("  Mean lambda=%.6f  Mean active set=%.2f\n", ...
        mean(lambda_selected,"omitnan"), mean(active_set_size,"omitnan"));

save("results_baseline_lasso_logchange.mat", ...
    "y_true_z_std_valid","y_true_z_raw_valid","y_true_level_valid", ...
    "yhat_z_std_valid","yhat_z_raw_valid","yhat_level_valid", ...
    "error_z_std","error_z_raw","error_level", ...
    "MSFE_z_std","RMSE_z_std","MAE_z_std", ...
    "MSFE_z_raw","RMSE_z_raw","MAE_z_raw", ...
    "MSFE_level","RMSE_level","MAE_level","MAPE_level", ...
    "timestamps_valid","lambda_selected","active_set_size","W","K","numLambda");

fprintf("Saved: results_baseline_lasso_logchange.mat\n");
