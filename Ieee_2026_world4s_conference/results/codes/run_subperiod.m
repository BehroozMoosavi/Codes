function T_sub = run_subperiod(varen_file, staticen_file)
%% run_subperiod.m  |  FILE 27 of 30
%
%  R4: Subperiod stability. RMSE/MAE/MAPE by year 2020-2023.
%  No re-estimation — subsets saved errors by year.
%
%  Inputs:  varen_file, staticen_file paths
%  Output:  T_sub MATLAB table

    SV=load(varen_file); SS=load(staticen_file);
    [eV,eS,tC]=align_errors(SV,SS);
    yrs=year(tC);
    [~,ia]=ismember(tC,SV.timestamps_valid);
    yt=SV.y_true_level_valid(ia);

    years=2020:2023; n=4;
    Year=years(:); N_obs=NaN(n,1);
    RV=NaN(n,1); RS=NaN(n,1); MV=NaN(n,1); MS=NaN(n,1);
    PV=NaN(n,1); PS=NaN(n,1); dR=NaN(n,1); dM=NaN(n,1);

    fprintf("R4: Subperiod stability\n");
    fprintf("  %-6s  %-6s  %-12s  %-12s  %-12s\n","Year","N","RMSE VA-REN","RMSE Static","Delta%%");

    for i=1:n
        yr=years(i); idx=yrs==yr&~isnan(eV)&~isnan(eS)&yt>0;
        if sum(idx)<10, fprintf("  %d insufficient\n",yr); continue; end
        ev=eV(idx); es=eS(idx); y=yt(idx);
        N_obs(i)=sum(idx);
        RV(i)=sqrt(mean(ev.^2)); RS(i)=sqrt(mean(es.^2));
        MV(i)=mean(abs(ev)); MS(i)=mean(abs(es));
        PV(i)=mean(abs(ev./y))*100; PS(i)=mean(abs(es./y))*100;
        dR(i)=100*(RS(i)-RV(i))/RS(i); dM(i)=100*(MS(i)-MV(i))/MS(i);
        fprintf("  %-6d  %-6d  %-12.2f  %-12.2f  %+.2f%%\n",yr,N_obs(i),RV(i),RS(i),dR(i));
    end

    fprintf("\n  VA-REN positive in %d/4 years.\n\n", sum(dR>0,"omitnan"));

    T_sub=table(Year,N_obs,RV,RS,MV,MS,PV,PS,dR,dM, ...
        "VariableNames",["Year","N","RMSE_VAREN","RMSE_StaticEN", ...
        "MAE_VAREN","MAE_StaticEN","MAPE_VAREN","MAPE_StaticEN", ...
        "Delta_RMSE_pct","Delta_MAE_pct"]);
end
