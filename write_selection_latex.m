function write_selection_latex(result, outDir)
%% write_selection_latex.m  |  FILE 23 of 30
%
%  Writes two LaTeX tables from regime_selection output.
%  table_7a_selection_top10.tex  — top-10 ramp vs stable side by side
%  table_7b_selection_groups.tex — group-level summary

    if nargin < 2 || isempty(outDir), outDir = "output/tables"; end
    if ~exist(outDir,"dir"), mkdir(outDir); end

    %% Table 7a: top-10
    f7a  = fullfile(outDir,"table_7a_selection_top10.tex");
    fid  = fopen(f7a,"w");
    R10  = result.top10_ramp;
    S10  = result.top10_stable;
    n    = min(height(R10),height(S10));

    fprintf(fid,"\\begin{table}[t]\\centering\\small\n");
    fprintf(fid,"\\caption{Top-10 predictors by regime-conditional selection differential ");
    fprintf(fid,"$\\Delta_j=f_j^{\\mathrm{ramp}}-f_j^{\\mathrm{stable}}$. ");
    fprintf(fid,"$^{*}$ = 95\\%% bootstrap CI excludes zero.}\n");
    fprintf(fid,"\\label{tab:selection_top10}\n");
    fprintf(fid,"\\begin{tabular}{@{}p{2.8cm}lrrr@{\\quad}p{2.8cm}lrrr@{}}\\toprule\n");
    fprintf(fid,"\\multicolumn{4}{c}{\\textit{Ramp-preferred}} & & \\multicolumn{4}{c}{\\textit{Stable-preferred}}\\\\\n");
    fprintf(fid,"\\cmidrule(r){1-4}\\cmidrule(l){6-9}\n");
    fprintf(fid,"Predictor & Group & $f^{\\mathrm{ramp}}$ & $\\Delta_j$ & & ");
    fprintf(fid,"Predictor & Group & $f^{\\mathrm{stable}}$ & $\\Delta_j$\\\\\n\\midrule\n");
    for i=1:n
        r=R10(i,:); s=S10(i,:);
        rs=""; if r.SigRamp,   rs="$^{*}$"; end
        ss=""; if s.SigStable, ss="$^{*}$"; end
        fprintf(fid,"%s%s & %s & %.3f & %+.3f & & %s%s & %s & %.3f & %+.3f\\\\\n", ...
            le(r.Predictor),rs,le(r.Group),r.FreqRamp,r.Delta, ...
            le(s.Predictor),ss,le(s.Group),s.FreqStable,s.Delta);
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n", f7a);

    %% Table 7b: group summary
    f7b = fullfile(outDir,"table_7b_selection_groups.tex");
    fid = fopen(f7b,"w");
    Tg  = result.T_group;

    fprintf(fid,"\\begin{table}[t]\\centering\n");
    fprintf(fid,"\\caption{Group-level regime-conditional selection frequencies. ");
    fprintf(fid,"$\\overline{\\Delta}$= group mean of $\\Delta_j$.}\n");
    fprintf(fid,"\\label{tab:selection_groups}\n");
    fprintf(fid,"\\begin{tabular}{@{}lrrrrrrrr@{}}\\toprule\n");
    fprintf(fid,"Group & $p$ & $\\bar f^{\\text{full}}$ & $\\bar f^{\\text{stable}}$ & ");
    fprintf(fid,"$\\bar f^{\\text{ramp}}$ & $\\overline{\\Delta}$ & ");
    fprintf(fid,"$N^*_{\\text{ramp}}$ & $N^*_{\\text{stable}}$\\\\\n\\midrule\n");
    for i=1:height(Tg)
        fprintf(fid,"%s & %d & %.3f & %.3f & %.3f & %+.3f & %d & %d\\\\\n", ...
            le(Tg.Group(i)), Tg.Size(i), Tg.FreqFull(i), Tg.FreqStable(i), ...
            Tg.FreqRamp(i), Tg.MeanDelta(i), Tg.NSig_Ramp(i), Tg.NSig_Stable(i));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\\end{table}\n");
    fclose(fid); fprintf("Saved: %s\n", f7b);

end

function s = le(x)
    s=string(x); s=replace(s,"_","\\_"); s=replace(s,"%","\\%"); s=replace(s,"&","\\&");
end
