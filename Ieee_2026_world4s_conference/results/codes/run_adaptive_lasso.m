function results_al = run_adaptive_lasso(dataFile, gamma_val)
%% run_adaptive_lasso.m  |  FILE 28 of 30
%
%  R5: Adaptive LASSO extended baseline.
%  Coefficient-specific weights: w_j = 1/|beta_pilot_j|^gamma
%  Pilot: Ridge with small fixed lambda.
%  Implementation via predictor rescaling trick.
%
%  Inputs:  dataFile, gamma_val (default 1)
%  Output:  results_al struct

    if nargin < 2, gamma_val = 1; end
    load(dataFile);
    W=720; K=5; numLambda=60; lp=0.01;
    X_all=X; y_all_std=y_eval_z_std(:);
    y_true=y_eval_level(:); y_lag1=y_lag1_eval_level(:);
    T_tot=length(y_all_std); p_dim=size(X_all,2);

    fprintf("Running Adaptive LASSO (gamma=%.1f)...\n",gamma_val);

    yv=NaN(T_tot,1); as=NaN(T_tot,1);
    for t=W+1:T_tot
        if mod(t,500)==0, fprintf("  step %d/%d\n",t,T_tot); end
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
            [yv(t),~]=back_transform(zh,y_lag1(t),z_mean,z_std);
            as(t)=sum(bw~=0);
        catch ME
            warning("AdapLASSO t=%d: %s",t,ME.message);
        end
    end

    v=~isnan(yv)&~isnan(y_true)&~isnan(y_all_std)&y_true>0&yv>0;
    err=y_true(v)-yv(v);
    results_al.RMSE_level=sqrt(mean(err.^2)); results_al.MAE_level=mean(abs(err));
    results_al.MAPE_level=mean(abs(err./y_true(v)))*100;
    results_al.mean_active_set=mean(as,"omitnan");
    results_al.gamma=gamma_val; results_al.error_level=err;
    results_al.timestamps_valid=timestamps_eval(v);
    results_al.y_true_level_valid=y_true(v); results_al.N=sum(v);

    fprintf("\nAdaptive LASSO (gamma=%.1f): RMSE=%.2f  MAE=%.2f  MAPE=%.4f%%  N=%d\n", ...
        gamma_val,results_al.RMSE_level,results_al.MAE_level, ...
        results_al.MAPE_level,results_al.N);

    save("results_adaptive_lasso.mat","-struct","results_al");
    fprintf("Saved: results_adaptive_lasso.mat\n\n");
end
