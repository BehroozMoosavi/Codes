%% ============================================================
%  comparison_tables.m
%  FILE 07 of 30
%
%  Quick Overleaf-ready LaTeX comparison document.
%  Produces 4 tables: compact, full accuracy, regime, diagnostics.
%
%  Input : all 5 results_*.mat files
%  Output: clean_method_comparison_tables_overleaf.tex
%% ============================================================

clear; clc;
fprintf("Creating comparison tables...\n");

methods = {
    "Seasonal naive", "results_baseline_seasonal_naive_logchange.mat";
    "Ridge",          "results_baseline_ridge_logchange.mat";
    "LASSO",          "results_baseline_lasso_logchange.mat";
    "Static EN",      "results_baseline_static_en_logchange.mat";
    "VA-REN",         "results_varen_logchange.mat"
};

M                 = size(methods,1);
Method            = strings(M,1);
N_obs             = NaN(M,1);
MSFE              = NaN(M,1); RMSE   = NaN(M,1);
MAE               = NaN(M,1); MAPE   = NaN(M,1);
Stable_RMSE       = NaN(M,1); Stable_MAE = NaN(M,1); Stable_MAPE = NaN(M,1);
Ramp_RMSE         = NaN(M,1); Ramp_MAE   = NaN(M,1); Ramp_MAPE   = NaN(M,1);
Mean_Alpha        = NaN(M,1);
Mean_Lambda       = NaN(M,1); Median_Lambda     = NaN(M,1);
Mean_Active_Set   = NaN(M,1); Median_Active_Set = NaN(M,1);
Imp_vs_StaticEN   = NaN(M,1);

for i = 1:M
    Method(i) = methods{i,1};
    fname     = methods{i,2};
    if ~isfile(fname), warning("Missing: %s",fname); continue; end
    S = load(fname);
    fprintf("Loaded: %s\n", fname);

    if isfield(S,"error_level"), N_obs(i)=numel(S.error_level);
    elseif isfield(S,"y_true_level_valid"), N_obs(i)=numel(S.y_true_level_valid); end

    MSFE(i) = gf(S,"MSFE_level"); RMSE(i) = gf(S,"RMSE_level");
    MAE(i)  = gf(S,"MAE_level");  MAPE(i) = gf(S,"MAPE_level");

    if isfield(S,"timestamps_valid") && isfield(S,"error_level") && ...
       isfield(S,"y_true_level_valid")
        hv = hour(S.timestamps_valid);
        e  = S.error_level(:); y = S.y_true_level_valid(:);
        si = hv>=10&hv<=14; ri = hv>=16&hv<=21;
        Stable_RMSE(i) = sqrt(mean(e(si).^2,"omitnan"));
        Stable_MAE(i)  = mean(abs(e(si)),"omitnan");
        Stable_MAPE(i) = mean(abs(e(si)./y(si)),"omitnan")*100;
        Ramp_RMSE(i)   = sqrt(mean(e(ri).^2,"omitnan"));
        Ramp_MAE(i)    = mean(abs(e(ri)),"omitnan");
        Ramp_MAPE(i)   = mean(abs(e(ri)./y(ri)),"omitnan")*100;
    end

    if     Method(i)=="Ridge",    Mean_Alpha(i)=0;
    elseif Method(i)=="LASSO",    Mean_Alpha(i)=1;
    elseif Method(i)=="Static EN",Mean_Alpha(i)=gf(S,"alpha_static");
    elseif Method(i)=="VA-REN",   Mean_Alpha(i)=gf(S,"mean_alpha"); end

    if isfield(S,"lambda_selected")
        Mean_Lambda(i)=mean(S.lambda_selected,"omitnan");
        Median_Lambda(i)=median(S.lambda_selected,"omitnan");
    elseif isfield(S,"lambda_valid")
        Mean_Lambda(i)=mean(S.lambda_valid,"omitnan");
        Median_Lambda(i)=median(S.lambda_valid,"omitnan"); end

    if isfield(S,"active_set_size")
        Mean_Active_Set(i)=mean(S.active_set_size,"omitnan");
        Median_Active_Set(i)=median(S.active_set_size,"omitnan");
    elseif isfield(S,"active_set_size_valid")
        Mean_Active_Set(i)=mean(S.active_set_size_valid,"omitnan");
        Median_Active_Set(i)=median(S.active_set_size_valid,"omitnan"); end
end

idx_s = find(Method=="Static EN",1);
if ~isempty(idx_s) && ~isnan(RMSE(idx_s))
    Imp_vs_StaticEN = 100*(RMSE(idx_s)-RMSE)./RMSE(idx_s);
end

fullT   = table(Method,N_obs,MSFE,RMSE,MAE,MAPE,Imp_vs_StaticEN);
regimeT = table(Method,Stable_RMSE,Stable_MAE,Stable_MAPE,Ramp_RMSE,Ramp_MAE,Ramp_MAPE);
diagT   = table(Method,Mean_Alpha,Mean_Lambda,Median_Lambda,Mean_Active_Set,Median_Active_Set);

fullT_s   = sortrows(fullT,"RMSE","ascend");
regimeT_s = sortrows(regimeT,"Ramp_RMSE","ascend");

disp(fullT_s); disp(regimeT_s); disp(diagT);

%% Write LaTeX
outFile = "clean_method_comparison_tables_overleaf.tex";
fid     = fopen(outFile,"w");

fprintf(fid,"\\documentclass[11pt]{article}\n");
fprintf(fid,"\\usepackage[margin=1in]{geometry}\n");
fprintf(fid,"\\usepackage{booktabs,threeparttable,float,amsmath}\n");
fprintf(fid,"\\renewcommand{\\arraystretch}{1.15}\n");
fprintf(fid,"\\begin{document}\n\n");

%% Table 1 compact
bRMSE=min(fullT_s.RMSE,[],"omitnan");
bMAE =min(fullT_s.MAE,[],"omitnan");
bsR  =min(regimeT.Stable_RMSE,[],"omitnan");
brR  =min(regimeT.Ramp_RMSE,[],"omitnan");

fprintf(fid,"\\begin{table}[H]\\centering\n");
fprintf(fid,"\\caption{Compact method comparison}\n\\label{tab:compact}\n");
fprintf(fid,"\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n");
fprintf(fid,"Model & RMSE & MAE & Stable RMSE & Ramp RMSE & $\\Delta$RMSE\\\\\\midrule\n");
for i=1:height(fullT_s)
    m=fullT_s.Method(i);
    jj=find(regimeT.Method==m,1);
    sR=NaN; rR=NaN;
    if ~isempty(jj), sR=regimeT.Stable_RMSE(jj); rR=regimeT.Ramp_RMSE(jj); end
    fprintf(fid,"%s & %s & %s & %s & %s & %s\\\\\n", ...
        mname(m), bb(fullT_s.RMSE(i),bRMSE,2), bb(fullT_s.MAE(i),bMAE,2), ...
        bb(sR,bsR,2), bb(rR,brR,2), fsgn(fullT_s.Imp_vs_StaticEN(i),2));
end
fprintf(fid,"\\bottomrule\\end{tabular}\n");
fprintf(fid,"\\begin{tablenotes}\\footnotesize\n");
fprintf(fid,"\\item RMSE and MAE in MW. $\\Delta$RMSE = \\%% improvement vs Static EN. Bold = best.\n");
fprintf(fid,"\\end{tablenotes}\\end{table}\n\n");

%% Table 2 full accuracy
bR2=min(fullT_s.RMSE,[],"omitnan");
bM2=min(fullT_s.MAE,[],"omitnan");
bP2=min(fullT_s.MAPE,[],"omitnan");
fprintf(fid,"\\begin{table}[H]\\centering\n");
fprintf(fid,"\\caption{Full-sample accuracy}\\label{tab:full}\n");
fprintf(fid,"\\begin{tabular}{@{}lrrrrrr@{}}\\toprule\n");
fprintf(fid,"Model & $N$ & MSFE & RMSE & MAE & MAPE & $\\Delta$RMSE\\\\\\midrule\n");
for i=1:height(fullT_s)
    fprintf(fid,"%s & %s & %s & %s & %s & %s & %s\\\\\n", ...
        mname(fullT_s.Method(i)), fint(fullT_s.N_obs(i)), ...
        fn(fullT_s.MSFE(i),2), bb(fullT_s.RMSE(i),bR2,2), ...
        bb(fullT_s.MAE(i),bM2,2), bb(fullT_s.MAPE(i),bP2,3), ...
        fsgn(fullT_s.Imp_vs_StaticEN(i),2));
end
fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n\n");

%% Table 3 regime
bsR2=min(regimeT_s.Stable_RMSE,[],"omitnan");
brR2=min(regimeT_s.Ramp_RMSE,[],"omitnan");
fprintf(fid,"\\begin{table}[H]\\centering\n");
fprintf(fid,"\\caption{Regime accuracy}\\label{tab:regime}\n");
fprintf(fid,"\\begin{tabular}{@{}lrrrrrrr@{}}\\toprule\n");
fprintf(fid,"& \\multicolumn{3}{c}{Stable (10--14 PT)} & \\multicolumn{3}{c}{Ramp (16--21 PT)}\\\\\n");
fprintf(fid,"\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n");
fprintf(fid,"Model & RMSE & MAE & MAPE & RMSE & MAE & MAPE\\\\\\midrule\n");
for i=1:height(regimeT_s)
    fprintf(fid,"%s & %s & %s & %s & %s & %s & %s\\\\\n", ...
        mname(regimeT_s.Method(i)), ...
        bb(regimeT_s.Stable_RMSE(i),bsR2,2), fn(regimeT_s.Stable_MAE(i),2), fn(regimeT_s.Stable_MAPE(i),3), ...
        bb(regimeT_s.Ramp_RMSE(i),brR2,2),   fn(regimeT_s.Ramp_MAE(i),2),   fn(regimeT_s.Ramp_MAPE(i),3));
end
fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n\n");

%% Table 4 diagnostics
fprintf(fid,"\\begin{table}[H]\\centering\n");
fprintf(fid,"\\caption{Regularization diagnostics}\\label{tab:diag}\n");
fprintf(fid,"\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n");
fprintf(fid,"Model & $\\bar{\\alpha}$ & Mean $\\lambda$ & Median $\\lambda$ & Mean act.\\ set & Median act.\\ set\\\\\\midrule\n");
for i=1:height(diagT)
    fprintf(fid,"%s & %s & %s & %s & %s & %s\\\\\n", ...
        mname(diagT.Method(i)), fn(diagT.Mean_Alpha(i),3), ...
        fn(diagT.Mean_Lambda(i),6), fn(diagT.Median_Lambda(i),6), ...
        fn(diagT.Mean_Active_Set(i),2), fn(diagT.Median_Active_Set(i),2));
end
fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n\n\\end{document}\n");
fclose(fid);

fprintf("Saved: %s\n", outFile);

%% Helpers
function v = gf(S,f)
    if isfield(S,f), v=S.(f); if numel(v)>1, v=v(1); end
    else, v=NaN; end
end
function s = fn(x,d)
    if isnan(x)||isempty(x), s="---"; else, s=sprintf("%.*f",d,x); end
end
function s = fint(x)
    if isnan(x)||isempty(x), s="---"; else, s=sprintf("%d",round(x)); end
end
function s = fsgn(x,d)
    if isnan(x)||isempty(x), s="---";
    elseif x>0, s=sprintf("+%.*f",d,x);
    else, s=sprintf("%.*f",d,x); end
end
function s = bb(v,bv,d)
    s=fn(v,d);
    if ~isnan(v)&&~isnan(bv)&&abs(v-bv)<1e-10, s="\\textbf{"+s+"}"; end
end
function s = mname(x)
    x=string(x);
    if x=="VA-REN", s="\\textbf{VA-REN}";
    else, s=replace(replace(x,"_","\\_"),"&","\\&"); end
end
