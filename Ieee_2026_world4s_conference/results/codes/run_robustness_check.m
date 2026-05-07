%% ============================================================
%  run_all_robustness_parallel.m
%
%  Runs all 5 robustness checks in PARALLEL on your Mac.
%
%  PREREQUISITES
%  -------------
%  All 5 model result files must exist:
%    results_baseline_seasonal_naive_logchange.mat
%    results_baseline_ridge_logchange.mat
%    results_baseline_lasso_logchange.mat
%    results_baseline_static_en_logchange.mat
%    results_varen_logchange.mat
%
%  HOW TO USE
%  ----------
%  Step 1:  >> data_clean
%  Step 2:  >> run_all_models_parallel
%  Step 3:  >> run_all_robustness_parallel    <-- this script
%  Step 4:  >> master_empirical
%
%  WHAT RUNS IN PARALLEL
%  ---------------------
%  R1: run_varen_window      window sizes [336 504 720 1080 1440]
%  R2: run_varen_volwindow   vol windows  [24 48 72 168]
%  R3: run_varen_alphamap    3 alpha mappings
%  R4: run_subperiod         yearly stability (fast, no re-estimation)
%  R5: run_adaptive_lasso    adaptive LASSO baseline
%
%  ESTIMATED TIME
%  --------------
%  R1 (heaviest): ~50-80 hours  (5 window sizes x 2 models)
%  R2           : ~20-30 hours
%  R3           : ~20-30 hours
%  R4           : < 5 minutes
%  R5           : ~20-25 hours
%  Total wall time ~ 50-80 hours (limited by R1)
%  vs ~130-165 hours sequential
%
%  OUTPUTS
%  -------
%  robustness_W.mat
%  robustness_Wv.mat
%  robustness_map.mat
%  subperiod.mat
%  subperiod_stability.csv
%  results_adaptive_lasso.mat
%% ============================================================

clear; clc;

fprintf("============================================================\n");
fprintf("  VA-REN — Running all 5 robustness checks in parallel\n");
fprintf("  Start: %s\n", datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf("============================================================\n\n");

%% ----------------------------------------------------------
%  Check all required files exist
%% ----------------------------------------------------------

dataFile = "caiso_final_logchange_dataset_2019_2023.mat";

required = {
    dataFile;
    "results_baseline_seasonal_naive_logchange.mat";
    "results_baseline_ridge_logchange.mat";
    "results_baseline_lasso_logchange.mat";
    "results_baseline_static_en_logchange.mat";
    "results_varen_logchange.mat"
};

fprintf("Checking required files...\n");
missing = false;
for i = 1:numel(required)
    if isfile(required{i})
        fprintf("  [OK]      %s\n", required{i});
    else
        fprintf("  [MISSING] %s\n", required{i});
        missing = true;
    end
end

if missing
    error("\nSome required files are missing.\nRun data_clean.m and run_all_models_parallel.m first.");
end
fprintf("\n");

%% ----------------------------------------------------------
%  Start parallel pool
%% ----------------------------------------------------------

nCores   = feature('numcores');
nWorkers = max(1, min(5, nCores - 1));

fprintf("CPU cores  : %d\n", nCores);
fprintf("Workers    : %d\n\n", nWorkers);

pool = gcp('nocreate');
if isempty(pool)
    fprintf("Starting parallel pool (%d workers)...\n\n", nWorkers);
    pool = parpool('local', nWorkers);
else
    fprintf("Using existing pool (%d workers).\n\n", pool.NumWorkers);
end

%% ----------------------------------------------------------
%  Define robustness checks
%% ----------------------------------------------------------

checks = {
    'R1_window',       'Window length sensitivity    W=[336 504 720 1080 1440]';
    'R2_volwindow',    'Volatility window sensitivity Wv=[24 48 72 168]';
    'R3_alphamap',     'Alpha mapping alternatives   [Linear Quadratic Threshold]';
    'R4_subperiod',    'Yearly subperiod stability   2020-2023';
    'R5_adaptive_lasso','Adaptive LASSO baseline     gamma=1'
};

nChecks = size(checks,1);
workDir = pwd;

fprintf("Robustness checks:\n");
for i = 1:nChecks
    fprintf("  %s: %s\n", checks{i,1}, checks{i,2});
end
fprintf("\n");

%% ----------------------------------------------------------
%  Storage
%% ----------------------------------------------------------

results  = cell(nChecks,1);
timings  = zeros(nChecks,1);
success  = false(nChecks,1);
errors   = cell(nChecks,1);

%% ----------------------------------------------------------
%  Run all checks in parallel
%% ----------------------------------------------------------

fprintf("Launching parallel robustness jobs...\n\n");
t_start_all = tic;

parfor i = 1:nChecks

    checkName = checks{i,1};
    t_start   = tic;

    fprintf("[Worker] Starting: %s  (%s)\n", checkName, datestr(now,'HH:MM:SS'));

    try

        switch checkName

            case 'R1_window'
                results{i} = robustness_window(workDir);

            case 'R2_volwindow'
                results{i} = robustness_volwindow(workDir);

            case 'R3_alphamap'
                results{i} = robustness_alphamap(workDir);

            case 'R4_subperiod'
                results{i} = robustness_subperiod(workDir);

            case 'R5_adaptive_lasso'
                results{i} = robustness_adaptive_lasso(workDir);

        end

        timings(i) = toc(t_start);
        success(i) = true;

        fprintf("[Worker] DONE: %s  (%.1f min)\n", checkName, timings(i)/60);

    catch ME

        timings(i) = toc(t_start);
        success(i) = false;
        errors{i}  = ME;

        fprintf("[Worker] FAILED: %s — %s\n", checkName, ME.message);

    end

end

t_total = toc(t_start_all);

%% ----------------------------------------------------------
%  Save all results
%% ----------------------------------------------------------

fprintf("\nSaving robustness results...\n");

if success(1) && ~isempty(results{1})
    results_W = results{1};
    save('robustness_W.mat', 'results_W');
    fprintf("  Saved: robustness_W.mat\n");
end

if success(2) && ~isempty(results{2})
    results_Wv = results{2};
    save('robustness_Wv.mat', 'results_Wv');
    fprintf("  Saved: robustness_Wv.mat\n");
end

if success(3) && ~isempty(results{3})
    results_map = results{3};
    save('robustness_map.mat', 'results_map');
    fprintf("  Saved: robustness_map.mat\n");
end

if success(4) && ~isempty(results{4})
    T_sub = results{4};
    save('subperiod.mat', 'T_sub');
    writetable(T_sub, 'subperiod_stability.csv');
    fprintf("  Saved: subperiod.mat + subperiod_stability.csv\n");
end

if success(5) && ~isempty(results{5})
    S_al = results{5};
    save('results_adaptive_lasso.mat', '-struct', 'S_al');
    fprintf("  Saved: results_adaptive_lasso.mat\n");
end

%% ----------------------------------------------------------
%  Print summary
%% ----------------------------------------------------------

fprintf("\n============================================================\n");
fprintf("  ROBUSTNESS PARALLEL RUN SUMMARY\n");
fprintf("  End: %s\n", datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf("  Total wall time: %.1f hours\n", t_total/3600);
fprintf("============================================================\n\n");

fprintf("  %-25s  %-10s  %-12s\n", "Check","Status","Time");
fprintf("  %-25s  %-10s  %-12s\n", "-----","------","----");

for i = 1:nChecks
    if success(i)
        fprintf("  %-25s  %-10s  %.1f min\n", checks{i,1}, "OK", timings(i)/60);
    else
        fprintf("  %-25s  %-10s  %.1f min  [%s]\n", ...
                checks{i,1}, "FAILED", timings(i)/60, errors{i}.message);
    end
end

fprintf("\n");

if all(success)
    fprintf("  All 5 robustness checks complete.\n");
    fprintf("  Run master_empirical.m for full analysis.\n\n");
else
    fprintf("  %d check(s) failed. Re-run failed checks individually.\n\n", sum(~success));
end

%% ============================================================
%  ROBUSTNESS FUNCTIONS  (self-contained, run on worker nodes)
%% ============================================================

function results_W = robustness_window(workDir)
%  R1: Re-runs VA-REN and Static EN for W in [336 504 720 1080 1440]

    addpath(workDir);
    dataFile  = fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat");
    W_grid    = [336, 504, 720, 1080, 1440];
    K         = 5;
    numLambda = 60;
    alphaGrid = [0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];

    load(dataFile);
    X_all    = X; y_all_std = y_eval_z_std(:);
    y_true   = y_eval_level(:); y_lag1 = y_lag1_eval_level(:);
    alpha_all= evaluationData.Alpha_VAREN(:); T_tot=length(y_all_std);

    results_W = struct();

    for wi = 1:numel(W_grid)
        W = W_grid(wi);
        fprintf("  R1: W=%d\n", W);

        %% Select alpha_static
        Xi=X_all(1:W,:); yi=y_all_std(1:W);
        vi=all(~isnan(Xi),2)&~isnan(yi); Xi=Xi(vi,:); yi=yi(vi);
        bMSE=Inf; as=0.50;
        for a=1:numel(alphaGrid)
            try
                [~,FI]=lasso(Xi,yi,"Alpha",alphaGrid(a),"CV",K, ...
                             "NumLambda",numLambda,"Standardize",false);
                if FI.MSE(FI.IndexMinMSE)<bMSE
                    bMSE=FI.MSE(FI.IndexMinMSE); as=alphaGrid(a);
                end
            catch, end
        end

        yv=NaN(T_tot,1); ys=NaN(T_tot,1);
        for t=W+1:T_tot
            Xtr=X_all((t-W):(t-1),:); ytr=y_all_std((t-W):(t-1));
            Xte=X_all(t,:);
            vt=all(~isnan(Xtr),2)&~isnan(ytr); Xtr=Xtr(vt,:); ytr=ytr(vt);
            if size(Xtr,1)<K+10||any(isnan(Xte)), continue; end
            at=max(0.001,min(1,alpha_all(t)));
            try
                [B,FI]=lasso(Xtr,ytr,"Alpha",at,"CV",K,"NumLambda",numLambda,"Standardize",false);
                idx=FI.IndexMinMSE; zh=FI.Intercept(idx)+Xte*B(:,idx);
                yv(t)=y_lag1(t)*exp((z_mean+z_std*zh)/100);
            catch, end
            try
                [B,FI]=lasso(Xtr,ytr,"Alpha",as,"CV",K,"NumLambda",numLambda,"Standardize",false);
                idx=FI.IndexMinMSE; zh=FI.Intercept(idx)+Xte*B(:,idx);
                ys(t)=y_lag1(t)*exp((z_mean+z_std*zh)/100);
            catch, end
        end

        v=~isnan(yv)&~isnan(ys)&~isnan(y_true)&y_true>0&yv>0&ys>0;
        ev=y_true(v)-yv(v); es=y_true(v)-ys(v);
        rv=sqrt(mean(ev.^2)); rs=sqrt(mean(es.^2));
        results_W(wi).W=W; results_W(wi).RMSE_varen=rv;
        results_W(wi).RMSE_static=rs; results_W(wi).MAE_varen=mean(abs(ev));
        results_W(wi).MAE_static=mean(abs(es)); results_W(wi).N=sum(v);
        results_W(wi).alpha_static=as;
        results_W(wi).delta_rmse_pct=100*(rs-rv)/rs;
        fprintf("    W=%d: VA-REN=%.2f  Static=%.2f  Delta=%+.2f%%\n", ...
                W,rv,rs,results_W(wi).delta_rmse_pct);
    end
end

function results_Wv = robustness_volwindow(workDir)
%  R2: Re-computes NormalizedVolatility with Wv in [24 48 72 168]

    addpath(workDir);
    dataFile = fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat");
    Wv_grid  = [24, 48, 72, 168];
    W=720; K=5; numLambda=60;

    load(dataFile);
    X_all=X; y_all_std=y_eval_z_std(:);
    y_true=y_eval_level(:); y_lag1=y_lag1_eval_level(:);
    T_tot=length(y_all_std);
    z_full=allData.Delta_Log_Load(:); yr_full=allData.Year(:);

    results_Wv = struct();

    for vi=1:numel(Wv_grid)
        Wv=Wv_grid(vi);
        fprintf("  R2: Wv=%d\n", Wv);

        rv_=movvar(z_full,[Wv-1,0],0,"omitnan");
        c19=rv_(yr_full==2019); c19=c19(~isnan(c19));
        if isempty(c19)||max(c19)==min(c19), continue; end
        Vn=(rv_-min(c19))./(max(c19)-min(c19));
        Vn=max(0,min(1,Vn));
        at_new=1-Vn;
        [fd,ia]=ismember(evaluationData.Timestamp,allData.Timestamp);
        if ~all(fd), continue; end
        alpha_eval=at_new(ia);

        yv=NaN(T_tot,1);
        for t=W+1:T_tot
            Xtr=X_all((t-W):(t-1),:); ytr=y_all_std((t-W):(t-1));
            Xte=X_all(t,:);
            vt=all(~isnan(Xtr),2)&~isnan(ytr); Xtr=Xtr(vt,:); ytr=ytr(vt);
            at=max(0.001,min(1,alpha_eval(t)));
            if size(Xtr,1)<K+10||any(isnan(Xte))||isnan(at), continue; end
            try
                [B,FI]=lasso(Xtr,ytr,"Alpha",at,"CV",K,"NumLambda",numLambda,"Standardize",false);
                idx=FI.IndexMinMSE; zh=FI.Intercept(idx)+Xte*B(:,idx);
                yv(t)=y_lag1(t)*exp((z_mean+z_std*zh)/100);
            catch, end
        end

        v=~isnan(yv)&~isnan(y_true)&y_true>0&yv>0;
        err=y_true(v)-yv(v);
        results_Wv(vi).Wv=Wv; results_Wv(vi).RMSE=sqrt(mean(err.^2));
        results_Wv(vi).MAE=mean(abs(err));
        results_Wv(vi).MAPE=mean(abs(err./y_true(v)))*100;
        results_Wv(vi).mean_alpha=mean(alpha_eval,"omitnan");
        results_Wv(vi).N=sum(v);
        fprintf("    Wv=%d: RMSE=%.2f  MAPE=%.4f%%\n", Wv,results_Wv(vi).RMSE,results_Wv(vi).MAPE);
    end
end

function results_map = robustness_alphamap(workDir)
%  R3: Three alpha-volatility mappings

    addpath(workDir);
    dataFile = fullfile(workDir, "caiso_final_logchange_dataset_2019_2023.mat");
    W=720; K=5; numLambda=60;

    load(dataFile);
    X_all=X; y_all_std=y_eval_z_std(:);
    y_true=y_eval_level(:); y_lag1=y_lag1_eval_level(:);
    vol=evaluationData.NormalizedVolatility(:); T_tot=length(y_all_std);

    maps={"Linear",    @(v)1-v;
          "Quadratic", @(v)(1-v).^2;
          "Threshold", @(v)double(v<1/3)*0.9+double(v>=1/3&v<2/3)*0.5+double(v>=2/3)*0.1};

    results_map=struct();
    for mi=1:size(maps,1)
        nm=maps{mi,1}; fn_=maps{mi,2};
        am=max(0.001,min(1,fn_(vol)));
        fprintf("  R3: mapping=%s\n", nm);

        yv=NaN(T_tot,1);
        for t=W+1:T_tot
            Xtr=X_all((t-W):(t-1),:); ytr=y_all_std((t-W):(t-1));
            Xte=X_all(t,:); at=am(t);
            vt=all(~isnan(Xtr),2)&~isnan(ytr); Xtr=Xtr(vt,:); ytr=ytr(vt);
            if size(Xtr,1)<K+10||any(isnan(Xte)), continue; end
            try
                [B,FI]=lasso(Xtr,ytr,"Alpha",at,"CV",K,"NumLambda",numLambda,"Standardize",false);
                idx=FI.IndexMinMSE; zh=FI.Intercept(idx)+Xte*B(:,idx);
                yv(t)=y_lag1(t)*exp((z_mean+z_std*zh)/100);
            catch, end
        end

        v=~isnan(yv)&~isnan(y_true)&y_true>0&yv>0;
        err=y_true(v)-yv(v);
        results_map(mi).mapping=nm; results_map(mi).RMSE=sqrt(mean(err.^2));
        results_map(mi).MAE=mean(abs(err));
        results_map(mi).MAPE=mean(abs(err./y_true(v)))*100;
        results_map(mi).mean_alpha=mean(am,"omitnan"); results_map(mi).N=sum(v);
        fprintf("    %s: RMSE=%.2f  MAPE=%.4f%%\n", nm,results_map(mi).RMSE,results_map(mi).MAPE);
    end
end

function T_sub = robustness_subperiod(workDir)
%  R4: Year-by-year RMSE — no re-estimation, just subsets errors

    addpath(workDir);
    SV = load(fullfile(workDir,"results_varen_logchange.mat"));
    SS = load(fullfile(workDir,"results_baseline_static_en_logchange.mat"));

    [tC,ia,ib] = intersect(SV.timestamps_valid(:), SS.timestamps_valid(:));
    eV  = SV.error_level(ia);
    eS  = SS.error_level(ib);
    yt  = SV.y_true_level_valid(ia);
    yrs = year(tC);

    years=2020:2023; n=4;
    Year=years(:); N_obs=NaN(n,1);
    RV=NaN(n,1); RS=NaN(n,1); MV=NaN(n,1); MS=NaN(n,1);
    PV=NaN(n,1); PS=NaN(n,1); dR=NaN(n,1); dM=NaN(n,1);

    for i=1:n
        yr=years(i); idx=yrs==yr&~isnan(eV)&~isnan(eS)&yt>0;
        if sum(idx)<10, continue; end
        ev=eV(idx); es=eS(idx); y=yt(idx);
        N_obs(i)=sum(idx);
        RV(i)=sqrt(mean(ev.^2)); RS(i)=sqrt(mean(es.^2));
        MV(i)=mean(abs(ev));     MS(i)=mean(abs(es));
        PV(i)=mean(abs(ev./y))*100; PS(i)=mean(abs(es./y))*100;
        dR(i)=100*(RS(i)-RV(i))/RS(i); dM(i)=100*(MS(i)-MV(i))/MS(i);
        fprintf("  R4: %d  RMSE VA-REN=%.2f  Static=%.2f  Delta=%+.2f%%\n", ...
                yr,RV(i),RS(i),dR(i));
    end

    T_sub = table(Year,N_obs,RV,RS,MV,MS,PV,PS,dR,dM, ...
        "VariableNames",["Year","N","RMSE_VAREN","RMSE_StaticEN", ...
        "MAE_VAREN","MAE_StaticEN","MAPE_VAREN","MAPE_StaticEN", ...
        "Delta_RMSE_pct","Delta_MAE_pct"]);
end

function S = robustness_adaptive_lasso(workDir)
%  R5: Adaptive LASSO baseline (gamma=1)

    addpath(workDir);
    dataFile  = fullfile(workDir,"caiso_final_logchange_dataset_2019_2023.mat");
    gamma_val = 1;
    W=720; K=5; numLambda=60; lp=0.01;

    load(dataFile);
    X_all=X; y_all_std=y_eval_z_std(:);
    y_true=y_eval_level(:); y_lag1=y_lag1_eval_level(:);
    T_tot=length(y_all_std); p_dim=size(X_all,2);

    yv=NaN(T_tot,1); as=NaN(T_tot,1);

    for t=W+1:T_tot
        if mod(t,1000)==0, fprintf("  R5: step %d/%d\n",t,T_tot); end
        Xtr=X_all((t-W):(t-1),:); ytr=y_all_std((t-W):(t-1));
        Xte=X_all(t,:);
        vt=all(~isnan(Xtr),2)&~isnan(ytr); Xtr=Xtr(vt,:); ytr=ytr(vt);
        if size(Xtr,1)<K+10||any(isnan(Xte)), continue; end
        try
            n_tr=size(Xtr,1);
            Xi=[ones(n_tr,1),Xtr]; P=eye(p_dim+1); P(1,1)=0;
            bf=(Xi'*Xi+lp*P)\(Xi'*ytr); bp=bf(2:end);
            wj=1./(abs(bp)+1e-6).^gamma_val; wj=wj/mean(wj);
            Xw=Xtr./wj'; Xtw=Xte./wj';
            [B,FI]=lasso(Xw,ytr,"Alpha",1,"CV",K,"NumLambda",numLambda,"Standardize",false);
            idx=FI.IndexMinMSE; bw=B(:,idx)./wj;
            zh=FI.Intercept(idx)+Xte*bw;
            yv(t)=y_lag1(t)*exp((z_mean+z_std*zh)/100);
            as(t)=sum(bw~=0);
        catch, end
    end

    v=~isnan(yv)&~isnan(y_true)&~isnan(y_all_std)&y_true>0&yv>0;
    err=y_true(v)-yv(v);
    S.RMSE_level=sqrt(mean(err.^2)); S.MAE_level=mean(abs(err));
    S.MAPE_level=mean(abs(err./y_true(v)))*100;
    S.mean_active_set=mean(as,"omitnan");
    S.gamma=gamma_val; S.error_level=err;
    S.timestamps_valid=timestamps_eval(v);
    S.y_true_level_valid=y_true(v); S.N=sum(v);
    fprintf("  R5 Adaptive LASSO: RMSE=%.2f  MAPE=%.4f%%\n", S.RMSE_level,S.MAPE_level);
end