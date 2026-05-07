%% ============================================================
%  ridge.m
%  FILE 03 of 30
%
%  Rolling Ridge Regression Baseline
%  Objective: min_beta ||y-Xb||^2 + lambda*||b||^2
%  W=720, K=5 blocked CV, lambdaGrid 60 log-spaced values.
%
%  Input : caiso_final_logchange_dataset_2019_2023.mat
%  Output: results_baseline_ridge_logchange.mat
%% ============================================================

clear; clc;
load("caiso_final_logchange_dataset_2019_2023.mat");
fprintf("Running Rolling Ridge baseline...\n");

W          = 720;
K          = 5;
lambdaGrid = logspace(-4,4,60);

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

for t = W+1:T

    if mod(t,500)==0, fprintf("  Ridge step %d/%d\n",t,T); end

    trainIdx   = (t-W):(t-1);
    X_train    = X_all(trainIdx,:);
    y_train    = y_all_std(trainIdx);
    X_test     = X_all(t,:);

    validTrain = all(~isnan(X_train),2) & ~isnan(y_train);
    X_train    = X_train(validTrain,:);
    y_train    = y_train(validTrain);

    if size(X_train,1) < K+10 || any(isnan(X_test)), continue; end

    X_train_i  = [ones(size(X_train,1),1), X_train];
    X_test_i   = [1, X_test];

    foldID     = makeBlockedFolds(size(X_train_i,1), K);
    bestLambda = lambdaGrid(1);
    bestCV     = Inf;

    for l = 1:numel(lambdaGrid)
        lambda = lambdaGrid(l);
        cvErr  = NaN(K,1);
        for k = 1:K
            valIdx = foldID==k; trIdx = ~valIdx;
            Xtr=X_train_i(trIdx,:); ytr=y_train(trIdx);
            Xval=X_train_i(valIdx,:); yval=y_train(valIdx);
            P=eye(size(Xtr,2)); P(1,1)=0;
            beta = (Xtr'*Xtr + lambda*P)\(Xtr'*ytr);
            cvErr(k) = mean((yval - Xval*beta).^2,"omitnan");
        end
        mCV = mean(cvErr,"omitnan");
        if mCV < bestCV, bestCV=mCV; bestLambda=lambda; end
    end

    lambda_selected(t) = bestLambda;
    P = eye(size(X_train_i,2)); P(1,1)=0;
    betaFinal = (X_train_i'*X_train_i + bestLambda*P)\(X_train_i'*y_train);

    yhat_z_std(t) = X_test_i * betaFinal;
    yhat_z_raw(t) = z_mean + z_std * yhat_z_std(t);
    yhat_level(t) = y_lag1_level(t) * exp(yhat_z_raw(t)/100);

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

fprintf("\nRidge Results: N=%d  RMSE=%.2f MW  MAE=%.2f MW  MAPE=%.4f%%\n", ...
        numel(error_level), RMSE_level, MAE_level, MAPE_level);

save("results_baseline_ridge_logchange.mat", ...
    "y_true_z_std_valid","y_true_z_raw_valid","y_true_level_valid", ...
    "yhat_z_std_valid","yhat_z_raw_valid","yhat_level_valid", ...
    "error_z_std","error_z_raw","error_level", ...
    "MSFE_z_std","RMSE_z_std","MAE_z_std", ...
    "MSFE_z_raw","RMSE_z_raw","MAE_z_raw", ...
    "MSFE_level","RMSE_level","MAE_level","MAPE_level", ...
    "timestamps_valid","lambda_selected","W","K","lambdaGrid");

fprintf("Saved: results_baseline_ridge_logchange.mat\n");

function foldID = makeBlockedFolds(n,K)
    foldID=zeros(n,1); edges=round(linspace(1,n+1,K+1));
    for k=1:K, foldID(edges(k):edges(k+1)-1)=k; end
end
