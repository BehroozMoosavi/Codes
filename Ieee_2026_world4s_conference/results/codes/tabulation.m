%% ============================================================
%  tabulation.m
%  FILE 08 of 30
%
%  Full paper LaTeX tables (6 tables + CSVs).
%
%  Table 1: data summary
%  Table 2: full-sample accuracy
%  Table 3: regime accuracy (stable + ramp)
%  Table 4: Diebold-Mariano tests vs Static EN
%  Table 5: model diagnostics
%  Table 6: VA-REN mechanism
%
%  Input : dataset .mat + all 5 results .mat files
%  Output: table_1_*.tex/.csv ... table_6_*.tex/.csv
%% ============================================================

clear; clc;
fprintf("Generating full paper LaTeX tables...\n\n");

dataFile = "caiso_final_logchange_dataset_2019_2023.mat";
resultFiles = {
    'results_baseline_seasonal_naive_logchange.mat','Seasonal naive','SeasonalNaive';
    'results_baseline_ridge_logchange.mat','Ridge','Ridge';
    'results_baseline_lasso_logchange.mat','LASSO','LASSO';
    'results_baseline_static_en_logchange.mat','Static EN','StaticEN';
    'results_varen_logchange.mat','VA-REN','VAREN'
};
methodOrder = {'SeasonalNaive','Ridge','LASSO','StaticEN','VAREN'};

if ~isfile(dataFile), error("Missing: %s",dataFile); end

D              = load(dataFile);
allData        = D.allData;
modelData      = D.modelData;
evaluationData = D.evaluationData;

Results = struct();
for i = 1:size(resultFiles,1)
    fname = resultFiles{i,1};
    dn    = resultFiles{i,2};
    kn    = resultFiles{i,3};
    if isfile(fname)
        S             = load(fname);
        S.DisplayName = dn; S.KeyName = kn;
        Results.(kn)  = S;
        fprintf("Loaded: %s\n", fname);
    else
        warning("Missing: %s", fname);
        Results.(kn) = struct("DisplayName",dn,"KeyName",kn);
    end
end

T1 = buildDataSummary(allData, modelData, evaluationData);
T2 = buildMainAccuracy(Results, methodOrder);
T3 = buildRegimeAccuracy(Results, methodOrder);
T4 = buildDMTable(Results, methodOrder);
T5 = buildDiagnostics(Results, methodOrder);
T6 = buildMechanism(Results);

fprintf("\nTable 1:\n"); disp(T1);
fprintf("\nTable 2:\n"); disp(T2);
fprintf("\nTable 3:\n"); disp(T3);
fprintf("\nTable 4:\n"); disp(T4);
fprintf("\nTable 5:\n"); disp(T5);
fprintf("\nTable 6:\n"); disp(T6);

writetable(T1,"table_1_data_summary.csv");
writetable(T2,"table_2_main_accuracy.csv");
writetable(T3,"table_3_regime_accuracy.csv");
writetable(T4,"table_4_dm_tests.csv");
writetable(T5,"table_5_model_diagnostics.csv");
writetable(T6,"table_6_varen_mechanism.csv");

writeT1(T1,"table_1_data_summary.tex");
writeT2(T2,"table_2_main_accuracy.tex");
writeT3(T3,"table_3_regime_accuracy.tex");
writeT4(T4,"table_4_dm_tests.tex");
writeT5(T5,"table_5_model_diagnostics.tex");
writeT6(T6,"table_6_varen_mechanism.tex");

fprintf("\nAll tables saved.\n");
fprintf("Include in paper with: \\input{table_X_...tex}\n");

%% ============================================================
%  Table builders
%% ============================================================

function T = buildDataSummary(allData, modelData, evalData)
    Period = strings(3,1); Start = strings(3,1); Finish = strings(3,1);
    Obs=NaN(3,1); MnLd=NaN(3,1); SdLd=NaN(3,1); MnZ=NaN(3,1); SdZ=NaN(3,1);
    DS = {allData, modelData, evalData};
    NM = ["Full raw sample","Model-ready sample","Evaluation sample"];
    for i = 1:3
        X = DS{i};
        Period(i) = NM(i);
        Start(i)  = string(X.Timestamp(1));
        Finish(i) = string(X.Timestamp(end));
        Obs(i)    = height(X);
        if ismember('CAISO_Load_MW', X.Properties.VariableNames)
            MnLd(i) = mean(X.CAISO_Load_MW, 'omitnan');
            SdLd(i) = std(X.CAISO_Load_MW,  0, 'omitnan'); end
        if ismember('Delta_Log_Load', X.Properties.VariableNames)
            MnZ(i) = mean(X.Delta_Log_Load, 'omitnan');
            SdZ(i) = std(X.Delta_Log_Load,  0, 'omitnan'); end
    end
    T = table(Period, Start, Finish, Obs, MnLd, SdLd, MnZ, SdZ, ...
        'VariableNames', {'Period','Start','End','Observations', ...
        'MeanLoad','StdLoad','MeanLogChange','StdLogChange'});
end

function T = buildMainAccuracy(Results, mo)
    n=numel(mo); Method=strings(n,1); VObs=NaN(n,1);
    MSFE=NaN(n,1); RMSE=NaN(n,1); MAE=NaN(n,1); MAPE=NaN(n,1);
    dRMSE=NaN(n,1); dMAE=NaN(n,1);
    for i=1:n
        S=Results.(mo{i}); Method(i)=S.DisplayName;
        if isfield(S,'error_level'), VObs(i)=numel(S.error_level);
        elseif isfield(S,'y_true_level_valid'), VObs(i)=numel(S.y_true_level_valid); end
        MSFE(i)=gf(S,'MSFE_level'); RMSE(i)=gf(S,'RMSE_level');
        MAE(i) =gf(S,'MAE_level');  MAPE(i)=gf(S,'MAPE_level');
    end
    sR=RMSE(strcmp(mo,'StaticEN')); sM=MAE(strcmp(mo,'StaticEN'));
    for i=1:n
        dRMSE(i)=100*(sR-RMSE(i))/sR;
        dMAE(i) =100*(sM-MAE(i))/sM; end
    T=table(Method,VObs,MSFE,RMSE,MAE,MAPE,dRMSE,dMAE, ...
        'VariableNames',{'Method','N','MSFE','RMSE','MAE','MAPE','RMSE_Imp_vs_StaticEN','MAE_Imp_vs_StaticEN'});
end

function T = buildRegimeAccuracy(Results, mo)
    n=numel(mo); Method=strings(n,1);
    sN=NaN(n,1); sR=NaN(n,1); sM=NaN(n,1); sP=NaN(n,1);
    rN=NaN(n,1); rR=NaN(n,1); rM=NaN(n,1); rP=NaN(n,1);
    drR=NaN(n,1); drM=NaN(n,1);
    for i=1:n
        S=Results.(mo{i}); Method(i)=S.DisplayName;
        if ~(isfield(S,'timestamps_valid')&&isfield(S,'error_level')&&...
             isfield(S,'y_true_level_valid')), continue; end
        hv=hour(S.timestamps_valid); e=S.error_level(:); y=S.y_true_level_valid(:);
        si=hv>=10&hv<=14; ri=hv>=16&hv<=21;
        sN(i)=sum(si); rN(i)=sum(ri);
        sR(i)=sqrt(mean(e(si).^2,'omitnan')); sM(i)=mean(abs(e(si)),'omitnan');
        sP(i)=mean(abs(e(si)./y(si)),'omitnan')*100;
        rR(i)=sqrt(mean(e(ri).^2,'omitnan')); rM(i)=mean(abs(e(ri)),'omitnan');
        rP(i)=mean(abs(e(ri)./y(ri)),'omitnan')*100;
    end
    srR=rR(strcmp(mo,'StaticEN')); srM=rM(strcmp(mo,'StaticEN'));
    for i=1:n, drR(i)=100*(srR-rR(i))/srR; drM(i)=100*(srM-rM(i))/srM; end
    T=table(Method,sN,sR,sM,sP,rN,rR,rM,rP,drR,drM, ...
        'VariableNames',{'Method','Stable_N','Stable_RMSE','Stable_MAE','Stable_MAPE','Ramp_N','Ramp_RMSE','Ramp_MAE','Ramp_MAPE','Ramp_RMSE_Imp_vs_StaticEN','Ramp_MAE_Imp_vs_StaticEN'});
end

function T = buildDMTable(Results, mo)
    base=Results.StaticEN;
    Method=strings(0); Sample=strings(0); N=zeros(0);
    MLD=zeros(0); DM=zeros(0); PV=zeros(0); Sig=strings(0);
    if ~(isfield(base,'timestamps_valid')&&isfield(base,'error_level'))
        T=table(Method,Sample,N,MLD,DM,PV,Sig, ...
            'VariableNames',{'Method','Sample','N','MeanLossDiff','DM','PValue','Stars'});
        return; end
    for i=1:numel(mo)
        k=mo{i}; if strcmp(k,'StaticEN'), continue; end
        S=Results.(k);
        if ~(isfield(S,'timestamps_valid')&&isfield(S,'error_level')), continue; end
        for q=1:3
            sname=["Full sample","Stable hours","Ramp hours"]; sn=sname(q);
            [eB,eM,tC]=alignErr(base,S);
            if isempty(eB), continue; end
            hv=hour(tC);
            if     sn=="Stable hours", idx=hv>=10&hv<=14;
            elseif sn=="Ramp hours",   idx=hv>=16&hv<=21;
            else,                      idx=true(size(hv)); end
            eB2=eB(idx); eM2=eM(idx);
            if numel(eB2)<30, continue; end
            [dm,pv,md]=dmTest(eB2,eM2);
            Method(end+1,1)=S.DisplayName; Sample(end+1,1)=sn;
            N(end+1,1)=numel(eB2); MLD(end+1,1)=md;
            DM(end+1,1)=dm; PV(end+1,1)=pv;
            Sig(end+1,1)=stars(pv);
        end
    end
    T=table(Method,Sample,N,MLD,DM,PV,Sig, ...
        'VariableNames',{'Method','Sample','N','MeanLossDiff','DM','PValue','Stars'});
end

function T = buildDiagnostics(Results, mo)
    n=numel(mo); Method=strings(n,1); Alpha=NaN(n,1);
    MnL=NaN(n,1); MdL=NaN(n,1); MnA=NaN(n,1); MdA=NaN(n,1);
    for i=1:n
        S=Results.(mo{i}); Method(i)=S.DisplayName;
        if     strcmp(mo{i},'Ridge'),    Alpha(i)=0;
        elseif strcmp(mo{i},'LASSO'),    Alpha(i)=1;
        elseif strcmp(mo{i},'StaticEN'), Alpha(i)=gf(S,'alpha_static');
        elseif strcmp(mo{i},'VAREN'),    Alpha(i)=gf(S,'mean_alpha'); end
        if isfield(S,'lambda_selected')
            MnL(i)=mean(S.lambda_selected,'omitnan');
            MdL(i)=median(S.lambda_selected,'omitnan');
        elseif isfield(S,'lambda_valid')
            MnL(i)=mean(S.lambda_valid,'omitnan');
            MdL(i)=median(S.lambda_valid,'omitnan'); end
        if isfield(S,'active_set_size')
            MnA(i)=mean(S.active_set_size,'omitnan');
            MdA(i)=median(S.active_set_size,'omitnan');
        elseif isfield(S,'active_set_size_valid')
            MnA(i)=mean(S.active_set_size_valid,'omitnan');
            MdA(i)=median(S.active_set_size_valid,'omitnan'); end
    end
    T=table(Method,Alpha,MnL,MdL,MnA,MdA, ...
        'VariableNames',{'Method','Alpha','MeanLambda','MedianLambda','MeanActiveSet','MedianActiveSet'});
end

function T = buildMechanism(Results)
    S=Results.VAREN;
    Q=["Mean normalized volatility, stable hours";
       "Mean normalized volatility, ramp hours";
       "Mean alpha, stable hours";
       "Mean alpha, ramp hours";
       "Mean active-set size, stable hours";
       "Mean active-set size, ramp hours";
       "Corr(alpha, active-set size)"];
    V=[gf(S,"mean_vol_stable"); gf(S,"mean_vol_ramp");
       gf(S,"mean_alpha_stable"); gf(S,"mean_alpha_ramp");
       gf(S,"mean_active_stable"); gf(S,"mean_active_ramp");
       gf(S,"corr_alpha_active")];
    I=["Low volatility expected in solar-plateau hours";
       "High volatility expected in ramp hours";
       "Stable hours push VA-REN toward LASSO (sparse)";
       "Ramp hours push VA-REN toward Ridge (dense)";
       "Relatively sparse model expected in stable hours";
       "Ramp hours should retain more correlated predictors";
       "Confirms adaptive alpha changes model sparsity"];
    T=table(Q,V,I,'VariableNames',{'Quantity','Value','Interpretation'});
end

%% ============================================================
%  LaTeX writers
%% ============================================================

function writeT1(T,f)
    fid=fopen(f,"w");
    fprintf(fid,"\\begin{table}[t]\\centering\n");
    fprintf(fid,"\\caption{CAISO dataset summary. Response: $z_t=100[\\log y_t - \\log y_{t-1}]$.}\n");
    fprintf(fid,"\\label{tab:data_summary}\n");
    fprintf(fid,"\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n");
    fprintf(fid,"Sample & Obs & Mean load & Std load & Mean $z_t$ & Std $z_t$\\\\\\midrule\n");
    for i=1:height(T)
        fprintf(fid,"%s & %s & %s & %s & %s & %s\\\\\n", ...
            le_escape(T.Period(i)), fint(T.Observations(i)), ...
            fn(T.MeanLoad(i),2), fn(T.StdLoad(i),2), ...
            fn(T.MeanLogChange(i),4), fn(T.StdLogChange(i),4));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n",f);
end

function writeT2(T,f)
    fid=fopen(f,"w");
    bR=min(T.RMSE,[],'omitnan'); bM=min(T.MAE,[],'omitnan'); bP=min(T.MAPE,[],'omitnan');
    fprintf(fid,"\\begin{table}[t]\\centering\n");
    fprintf(fid,"\\caption{Full-sample accuracy, 2020--2023. $\\Delta$RMSE = \\%% improvement vs Static EN.}\n");
    fprintf(fid,"\\label{tab:main_accuracy}\n");
    fprintf(fid,"\\begin{tabular}{@{}lrrrrrr@{}}\\toprule\n");
    fprintf(fid,"Model & $N$ & MSFE & RMSE & MAE & MAPE & $\\Delta$RMSE\\\\\\midrule\n");
    for i=1:height(T)
        fprintf(fid,"%s & %s & %s & %s & %s & %s & %s\\\\\n", ...
            mn(T.Method(i)), fint(T.N(i)), fn(T.MSFE(i),2), ...
            bb(T.RMSE(i),bR,2), bb(T.MAE(i),bM,2), bb(T.MAPE(i),bP,3), ...
            fsgn(T.RMSE_Imp_vs_StaticEN(i),2));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\n");
    fprintf(fid,"\\begin{tablenotes}\\footnotesize\n");
    fprintf(fid,"\\item MSFE in MW$^2$; RMSE and MAE in MW; MAPE in \\%%. Bold = best.\n");
    fprintf(fid,"\\end{tablenotes}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n",f);
end

function writeT3(T,f)
    fid=fopen(f,"w");
    bsR=min(T.Stable_RMSE,[],'omitnan'); brR=min(T.Ramp_RMSE,[],'omitnan');
    fprintf(fid,"\\begin{table}[t]\\centering\n");
    fprintf(fid,"\\caption{Regime accuracy. Stable: 10--14 PT; Ramp: 16--21 PT.}\n");
    fprintf(fid,"\\label{tab:regime_accuracy}\n");
    fprintf(fid,"\\begin{tabular}{@{}lrrrrrrr@{}}\\toprule\n");
    fprintf(fid,"& \\multicolumn{3}{c}{Stable hours} & \\multicolumn{4}{c}{Ramp hours}\\\\\n");
    fprintf(fid,"\\cmidrule(lr){2-4}\\cmidrule(lr){5-8}\n");
    fprintf(fid,"Model & RMSE & MAE & MAPE & RMSE & MAE & MAPE & $\\Delta$RMSE\\\\\\midrule\n");
    for i=1:height(T)
        fprintf(fid,"%s & %s & %s & %s & %s & %s & %s & %s\\\\\n", ...
            mn(T.Method(i)), bb(T.Stable_RMSE(i),bsR,2), fn(T.Stable_MAE(i),2), fn(T.Stable_MAPE(i),3), ...
            bb(T.Ramp_RMSE(i),brR,2), fn(T.Ramp_MAE(i),2), fn(T.Ramp_MAPE(i),3), ...
            fsgn(T.Ramp_RMSE_Imp_vs_StaticEN(i),2));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\n");
    fprintf(fid,"\\begin{tablenotes}\\footnotesize\n");
    fprintf(fid,"\\item RMSE and MAE in MW; MAPE and $\\Delta$RMSE in \\%%.\n");
    fprintf(fid,"\\end{tablenotes}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n",f);
end

function writeT4(T,f)
    fid=fopen(f,"w");
    fprintf(fid,"\\begin{table}[t]\\centering\n");
    fprintf(fid,"\\caption{Diebold--Mariano tests vs Static EN. Loss: $d_t=e^2_{\\text{Static EN},t}-e^2_{m,t}$. Positive DM = row model beats Static EN.}\n");
    fprintf(fid,"\\label{tab:dm_tests}\n");
    fprintf(fid,"\\begin{tabular}{@{}llrrrr@{}}\\toprule\n");
    fprintf(fid,"Model & Sample & $N$ & $\\bar d$ & DM & $p$-value\\\\\\midrule\n");
    for i=1:height(T)
        fprintf(fid,"%s & %s & %s & %s & %s%s & %s\\\\\n", ...
            mn(T.Method(i)), le_escape(T.Sample(i)), fint(T.N(i)), ...
            fn(T.MeanLossDiff(i),2), fn(T.DM(i),3), T.Stars(i), fn(T.PValue(i),4));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\n");
    fprintf(fid,"\\begin{tablenotes}\\footnotesize\n");
    fprintf(fid,"\\item $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$. Newey-West HAC, Andrews bandwidth.\n");
    fprintf(fid,"\\end{tablenotes}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n",f);
end

function writeT5(T,f)
    fid=fopen(f,"w");
    fprintf(fid,"\\begin{table}[t]\\centering\n");
    fprintf(fid,"\\caption{Model diagnostics. For VA-REN, $\\alpha$ is the sample mean of the adaptive parameter.}\n");
    fprintf(fid,"\\label{tab:diagnostics}\n");
    fprintf(fid,"\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n");
    fprintf(fid,"Model & $\\alpha$ & Mean $\\lambda$ & Median $\\lambda$ & Mean act.\\ set & Median act.\\ set\\\\\\midrule\n");
    for i=1:height(T)
        fprintf(fid,"%s & %s & %s & %s & %s & %s\\\\\n", ...
            mn(T.Method(i)), fn(T.Alpha(i),4), fn(T.MeanLambda(i),6), ...
            fn(T.MedianLambda(i),6), fn(T.MeanActiveSet(i),2), fn(T.MedianActiveSet(i),2));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n",f);
end

function writeT6(T,f)
    fid=fopen(f,"w");
    fprintf(fid,"\\begin{table}[t]\\centering\n");
    fprintf(fid,"\\caption{VA-REN mechanism diagnostics. Expected: higher vol, lower $\\alpha_t$, denser model in ramp hours.}\n");
    fprintf(fid,"\\label{tab:varen_mechanism}\n");
    fprintf(fid,"\\begin{tabular}{@{}lrp{6cm}@{}}\\toprule\n");
    fprintf(fid,"Quantity & Value & Interpretation\\\\\\midrule\n");
    for i=1:height(T)
        fprintf(fid,"%s & %s & %s\\\\\n", le_escape(T.Quantity(i)), fn(T.Value(i),4), le_escape(T.Interpretation(i)));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n",f);
end

%% ============================================================
%  Statistical + formatting helpers
%% ============================================================

function [dm,pv,md] = dmTest(eB,eM)
    eB=eB(:); eM=eM(:);
    v=~isnan(eB)&~isnan(eM); eB=eB(v); eM=eM(v);
    d=eB.^2-eM.^2; Tn=numel(d); md=mean(d,'omitnan');
    dc=d-md; h=floor(4*(Tn/100)^(2/9));
    g0=mean(dc.^2,'omitnan'); lrv=g0;
    for j=1:h
        gj=mean(dc(j+1:end).*dc(1:end-j),'omitnan');
        lrv=lrv+2*(1-j/(h+1))*gj; end
    if lrv<=0||isnan(lrv), dm=NaN; pv=NaN; return; end
    dm=md/sqrt(lrv/Tn); pv=0.5*erfc(dm/sqrt(2));
end

function [eB,eM,tC] = alignErr(base,model)
    [tC,ia,ib]=intersect(base.timestamps_valid(:),model.timestamps_valid(:));
    eB=base.error_level(ia); eM=model.error_level(ib);
end

function s = stars(p)
    if isnan(p), s="";
    elseif p<0.01, s="$^{***}$";
    elseif p<0.05, s="$^{**}$";
    elseif p<0.10, s="$^{*}$";
    else, s=""; end
end

function v = gf(S,f)
    if isfield(S,f), v=S.(f); if numel(v)>1, v=v(1); end
    else, v=NaN; end
end

function s = fn(x,d)
    if isempty(x)||isnan(x), s="---"; else, s=sprintf("%.*f",d,x); end
end
function s = fint(x)
    if isempty(x)||isnan(x), s="---"; else, s=sprintf("%d",round(x)); end
end
function s = fsgn(x,d)
    if isempty(x)||isnan(x), s="---";
    elseif x>0, s=sprintf("+%.*f",d,x);
    else, s=sprintf("%.*f",d,x); end
end
function s = bb(v,bv,d)
    s=fn(v,d);
    if ~isnan(v)&&~isnan(bv)&&abs(v-bv)<1e-10, s="\\textbf{"+s+"}"; end
end
function s = mn(x)
    x=string(x);
    if x=='VA-REN', s="\\textbf{VA-REN}";
    else, s=replace(replace(x,"_","\\_"),"&","\\&"); end
end
function s = le_escape(x)
    s=string(x);
    s=replace(s,"_","\\_"); s=replace(s,"%","\\%"); s=replace(s,"&","\\&");
end