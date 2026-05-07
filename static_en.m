%% ============================================================
%  static_en.m
%  FILE 05 of 30
%
%  Rolling Static Elastic Net Baseline
%  alpha_static selected ONCE on obs 1:W before any forecast.
%  Lambda re-selected by CV at every rolling step.
%  Primary comparison model for VA-REN.
%
%  Input : caiso_final_logchange_dataset_2019_2023.mat
%  Output: results_baseline_static_en_logchange.mat
%% ============================================================

clear; clc;
load("caiso_final_logchange_dataset_2019_2023.mat");
fprintf("Running Static Elastic Net baseline...\n");

W         = 720;
K         = 5;
numLambda = 60;
alphaGrid = [0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];

X_all        = X;
y_all_std    = y_eval_z_std(:);
y_all_raw    = y_eval_z_raw(:);
y_true_level = y_eval_level(:);
y_lag1_level = y_lag1_eval_level(:);
T            = length(y_all_std);

%% --- Step 1: select alpha_static on first W obs ---
fprintf("Selecting alpha_static on obs 1:%d...\n", W);

X_init = X_all(1:W,:); y_init = y_all_std(1:W);
vI     = all(~isnan(X_init),2) & ~isnan(y_init);
X_init = X_init(vI,:); y_init = y_init(vI);

if size(X_init,1) < K+10
    error("Not enough valid obs in initial window.");
end

alphaCV_MSE = NaN(numel(alphaGrid),1);
alphaCV_SE  = NaN(numel(alphaGrid),1);

for a = 1:numel(alphaGrid)
    try
        [~,FI] = lasso(X_init,y_init,"Alpha",alphaGrid(a),"CV",K, ...
                       "NumLambda",numLambda,"Standardize",false);
        alphaCV_MSE(a) = FI.MSE(FI.IndexMinMSE);
        alphaCV_SE(a)  = FI.SE(FI.IndexMinMSE);
        fprintf("  alpha=%.2f  CV MSE=%.6f\n", alphaGrid(a), alphaCV_MSE(a));
    catch ME
        warning("alpha %.2f failed: %s", alphaGrid(a), ME.message);
        alphaCV_MSE(a) = Inf;
    end
end

[~,bestIdx]  = min(alphaCV_MSE);
alpha_static = alphaGrid(bestIdx);
fprintf("Selected alpha_static = %.2f\n", alpha_static);

%% --- Step 2: rolling forecasts with fixed alpha_static ---
yhat_z_std      = NaN(T,1);
yhat_z_raw      = NaN(T,1);
yhat_level      = NaN(T,1);
lambda_selected = NaN(T,1);
active_set_size = NaN(T,1);

for t = W+1:T

    if mod(t,500)==0, fprintf("  Static EN step %d/%d\n",t,T); end

    trainIdx   = (t-W):(t-1);
    X_train    = X_all(trainIdx,:);
    y_train    = y_all_std(trainIdx);
    X_test     = X_all(t,:);

    validTrain = all(~isnan(X_train),2) & ~isnan(y_train);
    X_train    = X_train(validTrain,:);
    y_train    = y_train(validTrain);

    if size(X_train,1) < K+10 || any(isnan(X_test)), continue; end

    try
        [B,FI] = lasso(X_train,y_train,"Alpha",alpha_static,"CV",K, ...
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
        warning("Static EN t=%d: %s",t,ME.message);
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

fprintf("\nStatic EN Results: alpha=%.2f  N=%d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
        alpha_static, numel(error_level), RMSE_level, MAPE_level);

save("results_baseline_static_en_logchange.mat", ...
    "y_true_z_std_valid","y_true_z_raw_valid","y_true_level_valid", ...
    "yhat_z_std_valid","yhat_z_raw_valid","yhat_level_valid", ...
    "error_z_std","error_z_raw","error_level", ...
    "MSFE_z_std","RMSE_z_std","MAE_z_std", ...
    "MSFE_z_raw","RMSE_z_raw","MAE_z_raw", ...
    "MSFE_level","RMSE_level","MAE_level","MAPE_level", ...
    "timestamps_valid","lambda_selected","active_set_size", ...
    "alpha_static","alphaGrid","alphaCV_MSE","alphaCV_SE","W","K","numLambda");

fprintf("Saved: results_baseline_static_en_logchange.mat\n");
