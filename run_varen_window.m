function results_W = run_varen_window(dataFile, W_grid)
%% run_varen_window.m  |  FILE 24 of 30
%
%  R1: Window length sensitivity.
%  Re-runs VA-REN and Static EN for each W in W_grid.
%  alpha_static re-selected on first W obs at each W.
%
%  Inputs:  dataFile, W_grid (e.g. [336 504 720 1080 1440])
%  Output:  results_W struct array

    load(dataFile);
    K=5; numLambda=60;
    alphaGrid=[0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];
    X_all=X; y_all_std=y_eval_z_std(:);
    y_true=y_eval_level(:); y_lag1=y_lag1_eval_level(:);
    alpha_all=evaluationData.Alpha_VAREN(:); T_tot=length(y_all_std);
    results_W=struct();

    for wi=1:numel(W_grid)
        W=W_grid(wi);
        fprintf("\n=== R1 W=%d ===\n",W);

        Xi=X_all(1:W,:); yi=y_all_std(1:W);
        vi=all(~isnan(Xi),2)&~isnan(yi); Xi=Xi(vi,:); yi=yi(vi);
        bMSE=Inf; as=0.50;
        for a=1:numel(alphaGrid)
            try
                [~,FI]=lasso(Xi,yi,"Alpha",alphaGrid(a),"CV",K, ...
                             "NumLambda",numLambda,"Standardize",false);
                if FI.MSE(FI.IndexMinMSE)<bMSE
                    bMSE=FI.MSE(FI.IndexMinMSE); as=alphaGrid(a); end
            catch, end
        end
        fprintf("  alpha_static=%.2f\n",as);

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
                [yv(t),~]=back_transform(zh,y_lag1(t),z_mean,z_std);
            catch, end
            try
                [B,FI]=lasso(Xtr,ytr,"Alpha",as,"CV",K,"NumLambda",numLambda,"Standardize",false);
                idx=FI.IndexMinMSE; zh=FI.Intercept(idx)+Xte*B(:,idx);
                [ys(t),~]=back_transform(zh,y_lag1(t),z_mean,z_std);
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
        fprintf("  RMSE VA-REN=%.2f  Static=%.2f  Delta=%+.2f%%\n",rv,rs,results_W(wi).delta_rmse_pct);
    end

    fprintf("\nR1 Summary:\n");
    for wi=1:numel(results_W)
        fprintf("  W=%-5d  VA-REN=%.2f  Static=%.2f  Delta=%+.2f%%\n", ...
            results_W(wi).W,results_W(wi).RMSE_varen, ...
            results_W(wi).RMSE_static,results_W(wi).delta_rmse_pct);
    end
end
