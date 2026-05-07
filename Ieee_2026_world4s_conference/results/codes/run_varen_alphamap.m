function results_map = run_varen_alphamap(dataFile)
%% run_varen_alphamap.m  |  FILE 26 of 30
%
%  R3: Alpha mapping sensitivity.
%  Tests: (1) Linear alpha=1-V  (2) Quadratic=(1-V)^2  (3) Threshold step
%
%  Input:  dataFile
%  Output: results_map struct array

    load(dataFile);
    W=720; K=5; numLambda=60;
    X_all=X; y_all_std=y_eval_z_std(:);
    y_true=y_eval_level(:); y_lag1=y_lag1_eval_level(:);
    vol=evaluationData.NormalizedVolatility(:); T_tot=length(y_all_std);

    maps = {"Linear",    @(v) 1-v;
            "Quadratic", @(v) (1-v).^2;
            "Threshold", @(v) double(v<1/3)*0.9 + double(v>=1/3&v<2/3)*0.5 + double(v>=2/3)*0.1};

    results_map=struct();
    for mi=1:size(maps,1)
        nm=maps{mi,1}; fn_=maps{mi,2};
        am=max(0.001,min(1,fn_(vol)));
        fprintf("\n=== R3 mapping=%s (mean alpha=%.4f) ===\n",nm,mean(am,"omitnan"));

        yv=NaN(T_tot,1); au=NaN(T_tot,1);
        for t=W+1:T_tot
            Xtr=X_all((t-W):(t-1),:); ytr=y_all_std((t-W):(t-1));
            Xte=X_all(t,:); at=am(t);
            vt=all(~isnan(Xtr),2)&~isnan(ytr); Xtr=Xtr(vt,:); ytr=ytr(vt);
            if size(Xtr,1)<K+10||any(isnan(Xte)), continue; end
            try
                [B,FI]=lasso(Xtr,ytr,"Alpha",at,"CV",K,"NumLambda",numLambda,"Standardize",false);
                idx=FI.IndexMinMSE; zh=FI.Intercept(idx)+Xte*B(:,idx);
                [yv(t),~]=back_transform(zh,y_lag1(t),z_mean,z_std); au(t)=at;
            catch, end
        end

        v=~isnan(yv)&~isnan(y_true)&y_true>0&yv>0;
        err=y_true(v)-yv(v);
        results_map(mi).mapping=nm; results_map(mi).RMSE=sqrt(mean(err.^2));
        results_map(mi).MAE=mean(abs(err));
        results_map(mi).MAPE=mean(abs(err./y_true(v)))*100;
        results_map(mi).mean_alpha=mean(au,"omitnan"); results_map(mi).N=sum(v);
        fprintf("  RMSE=%.2f  MAE=%.2f  MAPE=%.4f%%\n", ...
            results_map(mi).RMSE,results_map(mi).MAE,results_map(mi).MAPE);
    end

    fprintf("\nR3 Summary:\n");
    for mi=1:numel(results_map)
        fprintf("  %-12s  RMSE=%.2f  MAE=%.2f  MAPE=%.4f%%\n", ...
            results_map(mi).mapping,results_map(mi).RMSE, ...
            results_map(mi).MAE,results_map(mi).MAPE);
    end
    fprintf("\nR3 complete.\n\n");
end
