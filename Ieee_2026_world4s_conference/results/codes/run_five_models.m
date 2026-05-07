%% ============================================================
%  run_all_models_parallel.m
%
%  Runs all 5 forecasting models in PARALLEL on your Mac.
%
%  Uses MATLAB parfor with a local parallel pool.
%  Each model runs on its own CPU core simultaneously.
%  Total time = time of the SLOWEST model (not sum of all).
%
%  PREREQUISITES
%  -------------
%  1. data_clean.m must have already been run once.
%     The file caiso_final_logchange_dataset_2019_2023.mat
%     must exist in the current folder.
%
%  2. Parallel Computing Toolbox must be installed.
%     Check: >> ver   (look for Parallel Computing Toolbox)
%
%  RUN ORDER
%  ---------
%  Step 1:  >> data_clean          (run once, ~10 min)
%  Step 2:  >> run_all_models_parallel   (runs all 5 at once)
%  Step 3:  >> master_empirical    (analysis, tables, figures)
%
%  WHAT GETS PRODUCED
%  ------------------
%  results_baseline_seasonal_naive_logchange.mat
%  results_baseline_ridge_logchange.mat
%  results_baseline_lasso_logchange.mat
%  results_baseline_static_en_logchange.mat
%  results_varen_logchange.mat
%
%  ESTIMATED TIME ON MAC
%  ----------------------
%  With 5 cores (one per model):
%    naive      :  < 1 minute
%    ridge      :  8-15 hours
%    lasso      :  15-25 hours
%    static_en  :  15-25 hours
%    va_ren     :  15-25 hours
%  Total wall time ~ 15-25 hours (limited by slowest model)
%  vs ~70-90 hours sequential
%% ============================================================

clear; clc;

fprintf("============================================================\n");
fprintf("  VA-REN — Running all 5 models in parallel\n");
fprintf("  Start: %s\n", datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf("============================================================\n\n");

%% ----------------------------------------------------------
%  Check dataset exists
%% ----------------------------------------------------------

dataFile = "caiso_final_logchange_dataset_2019_2023.mat";

if ~isfile(dataFile)
    error("Dataset not found: %s\n\nRun data_clean.m first.", dataFile);
end

fprintf("Dataset found: %s\n\n", dataFile);

%% ----------------------------------------------------------
%  Set up parallel pool
%  Uses as many workers as you have physical CPU cores.
%  Leave 1 core free for the OS.
%% ----------------------------------------------------------

nCores     = feature('numcores');
nWorkers   = max(1, min(5, nCores - 1));   % max 5 (one per model), leave 1 free

fprintf("CPU cores available : %d\n", nCores);
fprintf("Workers to use      : %d\n\n", nWorkers);

%% Start or reuse existing pool
pool = gcp('nocreate');

if isempty(pool)
    fprintf("Starting parallel pool with %d workers...\n", nWorkers);
    pool = parpool('local', nWorkers);
    fprintf("Pool started.\n\n");
else
    fprintf("Using existing pool with %d workers.\n\n", pool.NumWorkers);
end

%% ----------------------------------------------------------
%  Define all 5 models
%% ----------------------------------------------------------

models = {
    'naive',         'run_naive';
    'ridge',         'run_ridge';
    'lasso_method',  'run_lasso';
    'static_en',     'run_static_en';
    'va_ren',        'run_va_ren'
};

nModels  = size(models, 1);
workDir  = pwd;

fprintf("Models to run:\n");
for i = 1:nModels
    fprintf("  %d. %s\n", i, models{i,1});
end
fprintf("\n");

%% ----------------------------------------------------------
%  Storage for results and timing
%% ----------------------------------------------------------

results   = cell(nModels, 1);
timings   = zeros(nModels, 1);
errors    = cell(nModels, 1);
success   = false(nModels, 1);

%% ----------------------------------------------------------
%  Run all models in parallel via parfor
%  Each iteration runs one complete model independently.
%% ----------------------------------------------------------

fprintf("Launching parallel jobs...\n\n");

t_total_start = tic;

parfor i = 1:nModels

    modelName = models{i, 1};
    t_start   = tic;

    fprintf("[Worker %d] Starting: %s  (%s)\n", ...
            get_worker_id(), modelName, datestr(now,'HH:MM:SS'));

    try
        %% Each worker runs the full model function
        switch modelName

            case 'naive'
                results{i} = run_naive(workDir);

            case 'ridge'
                results{i} = run_ridge(workDir);

            case 'lasso_method'
                results{i} = run_lasso(workDir);

            case 'static_en'
                results{i} = run_static_en(workDir);

            case 'va_ren'
                results{i} = run_va_ren(workDir);

        end

        timings(i) = toc(t_start);
        success(i) = true;

        fprintf("[Worker %d] DONE: %s  (%.1f min)\n", ...
                get_worker_id(), modelName, timings(i)/60);

    catch ME

        timings(i) = toc(t_start);
        success(i) = false;
        errors{i}  = ME;

        fprintf("[Worker %d] FAILED: %s — %s\n", ...
                get_worker_id(), modelName, ME.message);

    end

end

t_total = toc(t_total_start);

%% ----------------------------------------------------------
%  Save all results to .mat files
%% ----------------------------------------------------------

fprintf("\nSaving result files...\n");

outFiles = {
    'results_baseline_seasonal_naive_logchange.mat';
    'results_baseline_ridge_logchange.mat';
    'results_baseline_lasso_logchange.mat';
    'results_baseline_static_en_logchange.mat';
    'results_varen_logchange.mat'
};

for i = 1:nModels
    if success(i) && ~isempty(results{i})
        S = results{i};
        save(outFiles{i}, '-struct', 'S');
        fprintf("  Saved: %s\n", outFiles{i});
    end
end

%% ----------------------------------------------------------
%  Print summary
%% ----------------------------------------------------------

fprintf("\n============================================================\n");
fprintf("  PARALLEL RUN SUMMARY\n");
fprintf("  End: %s\n", datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf("  Total wall time: %.1f hours\n", t_total/3600);
fprintf("============================================================\n\n");

fprintf("  %-20s  %-10s  %-12s  %s\n", "Model","Status","Time","RMSE (MW)");
fprintf("  %-20s  %-10s  %-12s  %s\n", "-----","------","----","---------");

for i = 1:nModels
    modelName = models{i,1};
    if success(i)
        rmse_str = "---";
        if ~isempty(results{i}) && isfield(results{i},'RMSE_level')
            rmse_str = sprintf("%.2f", results{i}.RMSE_level);
        end
        fprintf("  %-20s  %-10s  %8.1f min  %s\n", ...
                modelName, "OK", timings(i)/60, rmse_str);
    else
        fprintf("  %-20s  %-10s  %8.1f min  FAILED: %s\n", ...
                modelName, "FAILED", timings(i)/60, errors{i}.message);
    end
end

fprintf("\n");

n_ok     = sum(success);
n_failed = sum(~success);

if n_failed == 0
    fprintf("  All %d models completed successfully.\n", nModels);
    fprintf("  Run master_empirical.m to generate all tests and tables.\n\n");
else
    fprintf("  %d model(s) completed.  %d model(s) failed.\n", n_ok, n_failed);
    fprintf("  Check error messages above and re-run failed models individually.\n\n");
end

%% ============================================================
%  MODEL FUNCTIONS
%  Each function contains the complete model logic.
%  They return a struct that gets saved as the result .mat file.
%% ============================================================

function S = run_naive(workDir)
%  Seasonal Naive: yhat_t = y_{t-168}

    addpath(workDir);
    load(fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat"));

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

    S.y_true_level_valid = y_true_level(valid);
    S.y_lag1_level_valid = y_lag1_level(valid);
    S.z_true_raw_valid   = z_true_raw(valid);
    S.z_true_std_valid   = z_true_std(valid);
    S.yhat_level_valid   = yhat_level(valid);
    S.yhat_z_raw_valid   = yhat_z_raw(valid);
    S.yhat_z_std_valid   = yhat_z_std(valid);
    S.timestamps_valid   = timestamps_eval(valid);

    S.error_level = S.y_true_level_valid - S.yhat_level_valid;
    S.error_z_raw = S.z_true_raw_valid   - S.yhat_z_raw_valid;
    S.error_z_std = S.z_true_std_valid   - S.yhat_z_std_valid;

    S.MSFE_level = mean(S.error_level.^2,"omitnan");
    S.RMSE_level = sqrt(S.MSFE_level);
    S.MAE_level  = mean(abs(S.error_level),"omitnan");
    S.MAPE_level = mean(abs(S.error_level./S.y_true_level_valid),"omitnan")*100;
    S.MSFE_z_raw = mean(S.error_z_raw.^2,"omitnan");
    S.RMSE_z_raw = sqrt(S.MSFE_z_raw);
    S.MAE_z_raw  = mean(abs(S.error_z_raw),"omitnan");
    S.MSFE_z_std = mean(S.error_z_std.^2,"omitnan");
    S.RMSE_z_std = sqrt(S.MSFE_z_std);
    S.MAE_z_std  = mean(abs(S.error_z_std),"omitnan");

    fprintf("  Naive: N=%d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
            numel(S.error_level), S.RMSE_level, S.MAPE_level);
end

function S = run_ridge(workDir)
%  Rolling Ridge: W=720, K=5 blocked CV, lambdaGrid 60 values

    addpath(workDir);
    load(fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat"));

    W          = 720;
    K          = 5;
    lambdaGrid = logspace(-4, 4, 60);

    X_all        = X;
    y_all_std    = y_eval_z_std(:);
    y_all_raw    = y_eval_z_raw(:);
    y_true_level = y_eval_level(:);
    y_lag1_level = y_lag1_eval_level(:);
    T            = length(y_all_std);

    yhat_z_std      = NaN(T,1);
    yhat_z_raw      = NaN(T,1);
    yhat_level      = NaN(T,1);
    lambda_selected = NaN(T,1);

    for t = W+1:T

        if mod(t,1000)==0
            fprintf("  Ridge: step %d/%d\n", t, T);
        end

        trainIdx = (t-W):(t-1);
        X_train  = X_all(trainIdx,:);
        y_train  = y_all_std(trainIdx);
        X_test   = X_all(t,:);

        vt = all(~isnan(X_train),2) & ~isnan(y_train);
        X_train = X_train(vt,:); y_train = y_train(vt);

        if size(X_train,1) < K+10 || any(isnan(X_test)), continue; end

        X_i = [ones(size(X_train,1),1), X_train];
        X_t = [1, X_test];

        foldID = make_blocked_folds_local(size(X_i,1), K);
        bestL  = lambdaGrid(1);
        bestCV = Inf;

        for l = 1:numel(lambdaGrid)
            lam = lambdaGrid(l); cv = NaN(K,1);
            for k = 1:K
                vi = foldID==k; ti = ~vi;
                P  = eye(size(X_i(ti,:),2)); P(1,1)=0;
                b  = (X_i(ti,:)'*X_i(ti,:)+lam*P)\(X_i(ti,:)'*y_train(ti));
                cv(k) = mean((y_train(vi)-X_i(vi,:)*b).^2,"omitnan");
            end
            mc = mean(cv,"omitnan");
            if mc < bestCV, bestCV=mc; bestL=lam; end
        end

        P = eye(size(X_i,2)); P(1,1)=0;
        bf = (X_i'*X_i + bestL*P)\(X_i'*y_train);

        yhat_z_std(t) = X_t * bf;
        yhat_z_raw(t) = z_mean + z_std * yhat_z_std(t);
        yhat_level(t) = y_lag1_level(t) * exp(yhat_z_raw(t)/100);
        lambda_selected(t) = bestL;

    end

    S = build_result_struct(y_all_std, y_all_raw, y_true_level, ...
                            yhat_z_std, yhat_z_raw, yhat_level, ...
                            timestamps_eval, z_mean, z_std);
    S.lambda_selected = lambda_selected;
    S.W = W; S.K = K; S.lambdaGrid = lambdaGrid;

    fprintf("  Ridge: N=%d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
            numel(S.error_level), S.RMSE_level, S.MAPE_level);
end

function S = run_lasso(workDir)
%  Rolling LASSO: alpha=1, W=720, K=5, numLambda=60

    addpath(workDir);
    load(fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat"));

    W         = 720;
    K         = 5;
    numLambda = 60;

    X_all        = X;
    y_all_std    = y_eval_z_std(:);
    y_all_raw    = y_eval_z_raw(:);
    y_true_level = y_eval_level(:);
    y_lag1_level = y_lag1_eval_level(:);
    T            = length(y_all_std);

    yhat_z_std      = NaN(T,1);
    yhat_z_raw      = NaN(T,1);
    yhat_level      = NaN(T,1);
    lambda_selected = NaN(T,1);
    active_set_size = NaN(T,1);

    for t = W+1:T

        if mod(t,1000)==0
            fprintf("  LASSO: step %d/%d\n", t, T);
        end

        trainIdx = (t-W):(t-1);
        X_train  = X_all(trainIdx,:);
        y_train  = y_all_std(trainIdx);
        X_test   = X_all(t,:);

        vt = all(~isnan(X_train),2) & ~isnan(y_train);
        X_train = X_train(vt,:); y_train = y_train(vt);

        if size(X_train,1) < K+10 || any(isnan(X_test)), continue; end

        try
            [B,FI] = lasso(X_train, y_train, "Alpha",1, "CV",K, ...
                           "NumLambda",numLambda, "Standardize",false);
            idx   = FI.IndexMinMSE;
            beta  = B(:,idx);

            yhat_z_std(t) = FI.Intercept(idx) + X_test*beta;
            yhat_z_raw(t) = z_mean + z_std*yhat_z_std(t);
            yhat_level(t) = y_lag1_level(t)*exp(yhat_z_raw(t)/100);

            lambda_selected(t) = FI.Lambda(idx);
            active_set_size(t) = sum(beta~=0);
        catch
        end

    end

    S = build_result_struct(y_all_std, y_all_raw, y_true_level, ...
                            yhat_z_std, yhat_z_raw, yhat_level, ...
                            timestamps_eval, z_mean, z_std);
    S.lambda_selected = lambda_selected;
    S.active_set_size = active_set_size;
    S.W = W; S.K = K; S.numLambda = numLambda;

    fprintf("  LASSO: N=%d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
            numel(S.error_level), S.RMSE_level, S.MAPE_level);
end

function S = run_static_en(workDir)
%  Static Elastic Net: alpha selected once on first W obs, then fixed

    addpath(workDir);
    load(fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat"));

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

    %% Select alpha_static on first W obs
    Xi = X_all(1:W,:); yi = y_all_std(1:W);
    vi = all(~isnan(Xi),2) & ~isnan(yi);
    Xi = Xi(vi,:); yi = yi(vi);

    bestMSE = Inf; alpha_static = 0.50;
    for a = 1:numel(alphaGrid)
        try
            [~,FI] = lasso(Xi, yi, "Alpha",alphaGrid(a), "CV",K, ...
                           "NumLambda",numLambda, "Standardize",false);
            if FI.MSE(FI.IndexMinMSE) < bestMSE
                bestMSE      = FI.MSE(FI.IndexMinMSE);
                alpha_static = alphaGrid(a);
            end
        catch
        end
    end
    fprintf("  Static EN: alpha_static=%.2f\n", alpha_static);

    yhat_z_std      = NaN(T,1);
    yhat_z_raw      = NaN(T,1);
    yhat_level      = NaN(T,1);
    lambda_selected = NaN(T,1);
    active_set_size = NaN(T,1);

    for t = W+1:T

        if mod(t,1000)==0
            fprintf("  Static EN: step %d/%d\n", t, T);
        end

        trainIdx = (t-W):(t-1);
        X_train  = X_all(trainIdx,:);
        y_train  = y_all_std(trainIdx);
        X_test   = X_all(t,:);

        vt = all(~isnan(X_train),2) & ~isnan(y_train);
        X_train = X_train(vt,:); y_train = y_train(vt);

        if size(X_train,1) < K+10 || any(isnan(X_test)), continue; end

        try
            [B,FI] = lasso(X_train, y_train, "Alpha",alpha_static, "CV",K, ...
                           "NumLambda",numLambda, "Standardize",false);
            idx   = FI.IndexMinMSE;
            beta  = B(:,idx);

            yhat_z_std(t) = FI.Intercept(idx) + X_test*beta;
            yhat_z_raw(t) = z_mean + z_std*yhat_z_std(t);
            yhat_level(t) = y_lag1_level(t)*exp(yhat_z_raw(t)/100);

            lambda_selected(t) = FI.Lambda(idx);
            active_set_size(t) = sum(beta~=0);
        catch
        end

    end

    S = build_result_struct(y_all_std, y_all_raw, y_true_level, ...
                            yhat_z_std, yhat_z_raw, yhat_level, ...
                            timestamps_eval, z_mean, z_std);
    S.lambda_selected = lambda_selected;
    S.active_set_size = active_set_size;
    S.alpha_static    = alpha_static;
    S.alphaGrid       = alphaGrid;
    S.W = W; S.K = K; S.numLambda = numLambda;

    fprintf("  Static EN: N=%d  RMSE=%.2f MW  MAPE=%.4f%%\n", ...
            numel(S.error_level), S.RMSE_level, S.MAPE_level);
end

function S = run_va_ren(workDir)
%  VA-REN: alpha_t = 1 - NormalizedVolatility_t

    addpath(workDir);
    load(fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat"));

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

    yhat_z_std          = NaN(T,1);
    yhat_z_raw          = NaN(T,1);
    yhat_level          = NaN(T,1);
    lambda_selected     = NaN(T,1);
    alpha_selected      = NaN(T,1);
    vol_selected        = NaN(T,1);
    active_set_size     = NaN(T,1);
    selection_matrix    = false(T,p);

    for t = W+1:T

        if mod(t,1000)==0
            fprintf("  VA-REN: step %d/%d\n", t, T);
        end

        trainIdx = (t-W):(t-1);
        X_train  = X_all(trainIdx,:);
        y_train  = y_all_std(trainIdx);
        X_test   = X_all(t,:);
        alpha_t  = max(0.001, min(1.000, alpha_all(t)));
        vol_t    = vol_all(t);

        vt = all(~isnan(X_train),2) & ~isnan(y_train);
        X_train = X_train(vt,:); y_train = y_train(vt);

        if size(X_train,1) < K+10 || any(isnan(X_test)) || isnan(alpha_t)
            continue;
        end

        try
            [B,FI] = lasso(X_train, y_train, "Alpha",alpha_t, "CV",K, ...
                           "NumLambda",numLambda, "Standardize",false);
            idx  = FI.IndexMinMSE;
            beta = B(:,idx);

            yhat_z_std(t) = FI.Intercept(idx) + X_test*beta;
            yhat_z_raw(t) = z_mean + z_std*yhat_z_std(t);
            yhat_level(t) = y_lag1_level(t)*exp(yhat_z_raw(t)/100);

            lambda_selected(t)   = FI.Lambda(idx);
            alpha_selected(t)    = alpha_t;
            vol_selected(t)      = vol_t;
            active_set_size(t)   = sum(beta~=0);
            selection_matrix(t,:)= (beta~=0)';
        catch
        end

    end

    S = build_result_struct(y_all_std, y_all_raw, y_true_level, ...
                            yhat_z_std, yhat_z_raw, yhat_level, ...
                            timestamps_eval, z_mean, z_std);

    valid = S.valid_idx;
    S.alpha_valid           = alpha_selected(valid);
    S.volatility_valid      = vol_selected(valid);
    S.lambda_valid          = lambda_selected(valid);
    S.active_set_size_valid = active_set_size(valid);
    S.selection_matrix_valid= selection_matrix(valid,:);

    %% Regime metrics
    hv        = hour(S.timestamps_valid);
    stableIdx = hv>=10 & hv<=14;
    rampIdx   = hv>=16 & hv<=21;
    e         = S.error_level;
    yt        = S.y_true_level_valid;

    S.RMSE_level_stable = sqrt(mean(e(stableIdx).^2,"omitnan"));
    S.MAE_level_stable  = mean(abs(e(stableIdx)),"omitnan");
    S.MAPE_level_stable = mean(abs(e(stableIdx)./yt(stableIdx)),"omitnan")*100;
    S.RMSE_level_ramp   = sqrt(mean(e(rampIdx).^2,"omitnan"));
    S.MAE_level_ramp    = mean(abs(e(rampIdx)),"omitnan");
    S.MAPE_level_ramp   = mean(abs(e(rampIdx)./yt(rampIdx)),"omitnan")*100;

    av = S.alpha_valid; vv = S.volatility_valid; asv = S.active_set_size_valid;
    S.mean_alpha        = mean(av,"omitnan");
    S.median_alpha      = median(av,"omitnan");
    S.mean_alpha_stable = mean(av(stableIdx),"omitnan");
    S.mean_alpha_ramp   = mean(av(rampIdx),"omitnan");
    S.mean_vol_stable   = mean(vv(stableIdx),"omitnan");
    S.mean_vol_ramp     = mean(vv(rampIdx),"omitnan");
    S.mean_active_set   = mean(asv,"omitnan");
    S.median_active_set = median(asv,"omitnan");
    S.mean_active_stable= mean(asv(stableIdx),"omitnan");
    S.mean_active_ramp  = mean(asv(rampIdx),"omitnan");
    S.corr_alpha_active = corr(av, asv, "Rows","complete");
    S.stableIdx         = stableIdx;
    S.rampIdx           = rampIdx;
    S.W = W; S.K = K; S.numLambda = numLambda;

    fprintf("  VA-REN: N=%d  RMSE=%.2f MW  stable=%.2f  ramp=%.2f\n", ...
            numel(S.error_level), S.RMSE_level, S.RMSE_level_stable, S.RMSE_level_ramp);
end

%% ============================================================
%  SHARED HELPER FUNCTIONS
%% ============================================================

function S = build_result_struct(y_all_std, y_all_raw, y_true_level, ...
                                  yhat_z_std, yhat_z_raw, yhat_level, ...
                                  timestamps_eval, z_mean, z_std)
%  Filters to valid rows and computes standard accuracy metrics.
%  Returns a struct ready to be saved as a .mat file.

    valid = ~isnan(yhat_z_std)   & ~isnan(yhat_z_raw)   & ...
            ~isnan(yhat_level)   & ~isnan(y_all_std)     & ...
            ~isnan(y_all_raw)    & ~isnan(y_true_level)  & ...
            y_true_level > 0     & yhat_level > 0;

    S.valid_idx             = valid;
    S.y_true_z_std_valid    = y_all_std(valid);
    S.y_true_z_raw_valid    = y_all_raw(valid);
    S.y_true_level_valid    = y_true_level(valid);
    S.yhat_z_std_valid      = yhat_z_std(valid);
    S.yhat_z_raw_valid      = yhat_z_raw(valid);
    S.yhat_level_valid      = yhat_level(valid);
    S.timestamps_valid      = timestamps_eval(valid);

    S.error_z_std = S.y_true_z_std_valid - S.yhat_z_std_valid;
    S.error_z_raw = S.y_true_z_raw_valid - S.yhat_z_raw_valid;
    S.error_level = S.y_true_level_valid - S.yhat_level_valid;

    S.MSFE_level = mean(S.error_level.^2,"omitnan");
    S.RMSE_level = sqrt(S.MSFE_level);
    S.MAE_level  = mean(abs(S.error_level),"omitnan");
    S.MAPE_level = mean(abs(S.error_level./S.y_true_level_valid),"omitnan")*100;

    S.MSFE_z_raw = mean(S.error_z_raw.^2,"omitnan");
    S.RMSE_z_raw = sqrt(S.MSFE_z_raw);
    S.MAE_z_raw  = mean(abs(S.error_z_raw),"omitnan");

    S.MSFE_z_std = mean(S.error_z_std.^2,"omitnan");
    S.RMSE_z_std = sqrt(S.MSFE_z_std);
    S.MAE_z_std  = mean(abs(S.error_z_std),"omitnan");
end

function foldID = make_blocked_folds_local(n, K)
%  Contiguous (blocked) CV folds — no future leakage.
    foldID = zeros(n,1);
    edges  = round(linspace(1, n+1, K+1));
    for k = 1:K
        foldID(edges(k):edges(k+1)-1) = k;
    end
end

function id = get_worker_id()
%  Returns current parallel worker ID (or 0 for main thread).
    try
        w  = getCurrentWorker();
        id = w.ProcessId;
    catch
        id = 0;
    end
end