function plot_all_results(dataFile, varen_file, staticen_file, naive_file, ...
                           results_W, outDir)
%% plot_all_results.m  |  FILE 21 of 30
%
%  Generates 5 publication figures.
%
%  Fig 1: NormalizedVolatility + alpha_t dual-axis time series
%  Fig 2: Monthly rolling RMSE — VA-REN vs Static EN vs Naive
%  Fig 3: Regime RMSE grouped bar chart (stable vs ramp)
%  Fig 4: RMSE vs window W  (R1 robustness)
%  Fig 5: alpha_t distribution histogram
%
%  Inputs:
%    dataFile     — dataset .mat path
%    varen_file   — results_varen_logchange.mat path
%    staticen_file— results_baseline_static_en_logchange.mat path
%    naive_file   — results_baseline_seasonal_naive_logchange.mat path
%    results_W    — struct array from run_varen_window (or [])
%    outDir       — output folder (default "output/figures")

    if nargin < 6 || isempty(outDir), outDir = "output/figures"; end
    if ~exist(outDir,"dir"), mkdir(outDir); end

    set(groot,"DefaultAxesFontSize",11);
    set(groot,"DefaultLineLineWidth",1.5);

    D  = load(dataFile);
    SV = load(varen_file);
    SS = load(staticen_file);
    SN = load(naive_file);

    vol    = D.evaluationData.NormalizedVolatility;
    alpha_v= D.evaluationData.Alpha_VAREN;
    t_eval = D.evaluationData.Timestamp;

    %% Fig 1: Volatility + alpha dual-axis
    fig1 = figure("Position",[80 80 960 280]);
    yyaxis left;
    plot(t_eval, vol, "Color",[0.20 0.45 0.75],"LineWidth",0.7);
    ylabel("V_t^{norm}"); ylim([-0.04 1.12]);
    ax = gca; ax.YColor = [0.20 0.45 0.75];
    yyaxis right;
    plot(t_eval, alpha_v, "Color",[0.85 0.33 0.10],"LineWidth",0.7);
    ylabel("\alpha_t"); ylim([-0.04 1.12]);
    ax = gca; ax.YColor = [0.85 0.33 0.10];
    xlabel("Date");
    title("VA-REN adaptive parameter — evaluation period 2020--2023");
    legend(["V_t^{norm}","\alpha_t"],"Location","northeast","FontSize",9);
    grid on;
    save_fig(fig1, fullfile(outDir,"fig1_volatility"));

    %% Fig 2: Monthly rolling RMSE
    fig2 = figure("Position",[80 80 940 320]);
    mods  = {SV,SS,SN};
    mlabs = ["VA-REN","Static EN","Seasonal Naive"];
    mclrs = {[0.85 0.33 0.10],[0.20 0.45 0.75],[0.55 0.55 0.55]};
    mstyl = {"-","--",":"};
    hold on;
    for mi=1:3
        S=mods{mi}; ts=S.timestamps_valid; e=S.error_level;
        ms=unique(dateshift(ts,"start","month")); rm=NaN(numel(ms),1);
        for mo=1:numel(ms)
            idx=dateshift(ts,"start","month")==ms(mo);
            rm(mo)=sqrt(mean(e(idx).^2,"omitnan")); end
        plot(ms,rm,"Color",mclrs{mi},"LineStyle",mstyl{mi},"DisplayName",mlabs(mi));
    end
    hold off;
    xlabel("Month"); ylabel("RMSE (MW)");
    title("Monthly RMSE — VA-REN vs baselines");
    legend("Location","northwest","FontSize",10); grid on;
    save_fig(fig2, fullfile(outDir,"fig2_monthly_rmse"));

    %% Fig 3: Regime bar chart
    fig3 = figure("Position",[80 80 680 360]);
    mods3  = {SV,SS,SN};
    labs3  = ["VA-REN","Static EN","Naive"];
    rs = NaN(3,1); rr = NaN(3,1);
    for mi=1:3
        S=mods3{mi}; hv=hour(S.timestamps_valid); e=S.error_level;
        rs(mi)=sqrt(mean(e(hv>=10&hv<=14).^2,"omitnan"));
        rr(mi)=sqrt(mean(e(hv>=16&hv<=21).^2,"omitnan")); end
    b=bar([rs,rr],"grouped");
    b(1).FaceColor=[0.20 0.45 0.75]; b(1).EdgeColor="none";
    b(2).FaceColor=[0.85 0.33 0.10]; b(2).EdgeColor="none";
    xticks(1:3); xticklabels(labs3);
    ylabel("RMSE (MW)");
    legend(["Stable (10-14 PT)","Ramp (16-21 PT)"],"Location","northwest","FontSize",10);
    title("Regime-specific RMSE by method"); ylim([0,max([rs;rr])*1.18]); grid on; box on;
    save_fig(fig3, fullfile(outDir,"fig3_regime_bar"));

    %% Fig 4: Window sensitivity (optional)
    if ~isempty(results_W) && numel(results_W) > 1
        fig4 = figure("Position",[80 80 680 320]);
        Wv=[results_W.W]; rv=[results_W.RMSE_varen]; rs_=[results_W.RMSE_static];
        plot(Wv,rv,"-o","Color",[0.85 0.33 0.10],"DisplayName","VA-REN"); hold on;
        plot(Wv,rs_,"--s","Color",[0.20 0.45 0.75],"DisplayName","Static EN");
        xline(720,"k:","LineWidth",1.4,"Label","W=720","LabelHorizontalAlignment","left","FontSize",9);
        hold off;
        xlabel("W (hours)"); ylabel("RMSE (MW)");
        title("RMSE sensitivity to rolling window W");
        legend("Location","best","FontSize",10); grid on;
        save_fig(fig4, fullfile(outDir,"fig4_window_robust"));
    end

    %% Fig 5: Alpha distribution
    fig5 = figure("Position",[80 80 560 300]);
    histogram(alpha_v, 30, "FaceColor",[0.85 0.33 0.10], ...
              "EdgeColor","none","Normalization","probability");
    mn_a = mean(alpha_v,"omitnan");
    xline(mn_a,"k--","LineWidth",1.4, ...
          "Label",sprintf("Mean=%.3f",mn_a),"LabelHorizontalAlignment","right","FontSize",9);
    xlabel("\alpha_t"); ylabel("Relative frequency");
    title("Distribution of VA-REN adaptive parameter \alpha_t");
    grid on;
    save_fig(fig5, fullfile(outDir,"fig5_alpha_distribution"));

    fprintf("\nAll figures saved to: %s\n\n", outDir);
end