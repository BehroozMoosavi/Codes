function plot_regime_selection(result, outDir)
%% ============================================================
%  plot_regime_selection.m  (fixed)
%  Fixed: containers.Map keys must be char, not string.
%% ============================================================

    if nargin < 2 || isempty(outDir), outDir = 'output/figures'; end
    if ~exist(outDir,'dir'), mkdir(outDir); end

    set(groot,'DefaultAxesFontSize',11);

    %% Color map — keys must be char for containers.Map
    cmap = containers.Map( ...
        {'Lag log-change','Log-level anchor','Hour dummy', ...
         'Weekday dummy','Month dummy','Other'}, ...
        {[0.20 0.45 0.75],[0.85 0.33 0.10],[0.47 0.67 0.19], ...
         [0.63 0.08 0.18],[0.93 0.69 0.13],[0.50 0.50 0.50]});

    %% --------------------------------------------------------
    %  Figure A: diverging bar chart
    %% --------------------------------------------------------

    [ds, si] = sort(result.delta, 'descend');
    ns  = result.predictor_names(si);
    gs  = result.group(si);
    lo  = result.ci_lo(si);
    hi  = result.ci_hi(si);
    srp = result.sig_ramp_flag(si);
    sst = result.sig_stable_flag(si);
    p   = numel(ds);

    figA = figure('Position',[50 50 480 880]);
    hold on;

    for j = 1:p
        clr = getColor(cmap, char(gs(j)));
        barh(j, ds(j), 'FaceColor',clr, 'EdgeColor','none', 'BarWidth',0.75);
        plot([lo(j),hi(j)],[j,j],'-','Color',[0.25 0.25 0.25],'LineWidth',0.8);
        plot(lo(j),j,'|','Color',[0.25 0.25 0.25],'MarkerSize',3.5);
        plot(hi(j),j,'|','Color',[0.25 0.25 0.25],'MarkerSize',3.5);
    end

    xline(0,'k-','LineWidth',1.2);

    li = unique([1:5, p-4:p],'stable');
    for j = li
        xo = 0.008*sign(ds(j)); if xo==0, xo=0.008; end
        ha = 'left'; if ds(j)<0, ha='right'; end
        st = ''; if srp(j)||sst(j), st='*'; end
        text(ds(j)+xo, j, [char(ns(j)), st], 'FontSize',7.5, ...
             'HorizontalAlignment',ha, 'VerticalAlignment','middle');
    end

    gu = unique(gs,'stable');
    lh = gobjects(numel(gu),1);
    for g = 1:numel(gu)
        lh(g) = patch(NaN,NaN,getColor(cmap,char(gu(g))),'EdgeColor','none');
    end
    legend(lh, cellstr(gu), 'Location','southeast','FontSize',8.5);

    yticks([]); grid on; box on;
    set(gca,'YDir','normal');
    xlabel('\Delta_j  =  f_j^{ramp}  -  f_j^{stable}','FontSize',11);
    title({'Regime-conditional selection differential', ...
           'Positive \rightarrow selected more in ramp hours (16-21 PT)'}, ...
          'FontSize',11);
    xlim([min(ds)-0.02, max(ds)+0.02]);
    hold off;

    save_fig(figA, fullfile(outDir,'figA_selection_diverging'));

    %% --------------------------------------------------------
    %  Figure B: group bar chart
    %% --------------------------------------------------------

    Tg = result.T_group;
    ng = height(Tg);

    figB = figure('Position',[100 100 720 360]);
    b = bar([Tg.FreqStable, Tg.FreqRamp],'grouped');
    b(1).FaceColor = [0.20 0.45 0.75]; b(1).EdgeColor = 'none';
    b(2).FaceColor = [0.85 0.33 0.10]; b(2).EdgeColor = 'none';

    gw  = min(0.8, 2/(2+1.5));
    xrp = (1:ng) + gw/4;
    for i = 1:ng
        dg = Tg.FreqRamp(i) - Tg.FreqStable(i);
        text(xrp(i), Tg.FreqRamp(i)+0.01, sprintf('%+.4f',dg), ...
             'HorizontalAlignment','center','FontSize',9, ...
             'FontWeight','bold','Color',[0.3 0.3 0.3]);
    end

    xticks(1:ng);
    xticklabels(cellstr(Tg.Group));
    xtickangle(18);
    ylabel('Mean selection frequency');
    legend({'Stable (10-14 PT)','Ramp (16-21 PT)'}, ...
           'Location','northwest','FontSize',10);
    title('Group-level selection frequency by regime','FontSize',11);
    ylim([0, max([Tg.FreqStable; Tg.FreqRamp])*1.22]);
    grid on; box on;

    save_fig(figB, fullfile(outDir,'figB_selection_group'));

    %% --------------------------------------------------------
    %  Figure C: heatmap
    %% --------------------------------------------------------

    if isfield(result,'sel_by_hour_group') && ~isempty(result.sel_by_hour_group)

        H   = result.sel_by_hour_group;
        go  = result.group_names_ordered;
        ngo = numel(go);

        figC = figure('Position',[100 100 820 300]);
        imagesc(0:23, 1:ngo, H');
        colormap(flipud(hot));
        cb = colorbar;
        cb.Label.String = 'Selection frequency';
        cb.FontSize = 10;
        clim([0 1]);

        yticks(1:ngo);
        yticklabels(cellstr(go));
        xticks(0:2:23);
        xlabel('Hour of day (PT)','FontSize',11);
        title('VA-REN: group \times hour selection frequency','FontSize',11);

        hold on;
        xline(9.5,'w--','LineWidth',2);
        xline(14.5,'w--','LineWidth',2);
        xline(15.5,'w:','LineWidth',2);
        xline(21.5,'w:','LineWidth',2);
        text(12,  ngo+0.65,'Stable','Color','w','HorizontalAlignment','center','FontSize',9);
        text(18.5,ngo+0.65,'Ramp',  'Color','w','HorizontalAlignment','center','FontSize',9);
        hold off;
        set(gca,'YDir','reverse');

        save_fig(figC, fullfile(outDir,'figC_selection_heatmap'));

    end

end

function clr = getColor(cmap, g)
    if isKey(cmap, g)
        clr = cmap(g);
    else
        clr = [0.5 0.5 0.5];
    end
end