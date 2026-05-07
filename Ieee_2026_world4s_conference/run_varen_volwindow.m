function results_Wv = run_varen_volwindow(dataFile, Wv_grid)
%% run_varen_volwindow.m  |  FILE 25 of 30
%
%  R2: Volatility window sensitivity.
%  Recomputes NormalizedVolatility with each Wv, re-runs VA-REN.
%
%  Inputs:  dataFile, Wv_grid (e.g. [24 48 72 168])
%  Output:  results_Wv struct array

    load(dataFile);
    W=720; K=5; numLambda=60;
    X_all=X; y_all_std=y_eval_z_std(:);
    y_true=y_eval_level(:); y_lag1=y_lag1_eval_level(:);
    T_tot=length(y_all_std);
    z_full=allData.Delta_Log_Load(:); yr_full=allData.Year(:);
    results_Wv=struct();

    for vi=1:numel(Wv_grid)
        Wv=Wv_grid(vi);
        fprintf("\n=== R2 Wv=%d ===\n",Wv);

        rv=movvar(z_full,[Wv-1,0],0,"omitnan");
        c19=rv(yr_full==2019); c19=c19(~isnan(c19));
        if isempty(c19)||max(c19)==min(c19), warning("Skip Wv=%d",Wv); continue; end
        Vn=(rv-min(c19))./(max(c19)-min(c19));
        Vn=max(0,min(1,Vn));
        at_new=1-Vn;

        [fd,ia]=ismember(evaluationData.Timestamp,allData.Timestamp);
        if ~all(fd), warning("Some eval timestamps not found."); end
        alpha_eval=at_new(ia);

        yv=NaN(T_tot,1); au=NaN(T_tot,1);
        for t=W+1:T_tot
            Xtr=X_all((t-W):(t-1),:); ytr=y_all_std((t-W):(t-1));
            Xte=X_all(t,:);
            vt=all(~isnan(Xtr),2)&~isnan(ytr); Xtr=Xtr(vt,:); ytr=ytr(vt);
            at=max(0.001,min(1,alpha_eval(t)));
            if size(Xtr,1)<K+10||any(isnan(Xte))||isnan(at), continue; end
            try
                [B,FI]=lasso(Xtr,ytr,"Alpha",at,"CV",K,"NumLambda",numLambda,"Standardize",false);
                idx=FI.IndexMinMSE; zh=FI.Intercept(idx)+Xte*B(:,idx);
                [yv(t),~]=back_transform(zh,y_lag1(t),z_mean,z_std);
                au(t)=at;
            catch, end
        end

        v=~isnan(yv)&~isnan(y_true)&y_true>0&yv>0;
        err=y_true(v)-yv(v);
        results_Wv(vi).Wv=Wv; results_Wv(vi).RMSE=sqrt(mean(err.^2));
        results_Wv(vi).MAE=mean(abs(err));
        results_Wv(vi).MAPE=mean(abs(err./y_true(v)))*100;
        results_Wv(vi).mean_alpha=mean(au,"omitnan");
        results_Wv(vi).std_alpha=std(au,"omitnan");
        results_Wv(vi).N=sum(v);
        fprintf("  RMSE=%.2f  MAE=%.2f  MAPE=%.4f%%  mean_alpha=%.4f\n", ...
            results_Wv(vi).RMSE,results_Wv(vi).MAE, ...
            results_Wv(vi).MAPE,results_Wv(vi).mean_alpha);
    end
    fprintf("\nR2 complete.\n\n");
end
