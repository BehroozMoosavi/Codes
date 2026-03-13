% =========================================================================
% Section 4 Numerical Illustrations — Polished Final
% 540x400 figures, 400 DPI, white background, LaTeX labels
% =========================================================================
clear; clc; close all;

set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize',14);
set(groot,'defaultAxesLineWidth',1.2);
set(groot,'defaultLineLineWidth',2.2);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

rng(2025);
outdir = 'figures_privacy_mechanism';
if ~exist(outdir,'dir'), mkdir(outdir); end

% Distinct, colorblind-friendly palette
c1 = [0.00 0.45 0.74];   % blue
c2 = [0.85 0.33 0.10];   % red-orange
c3 = [0.93 0.69 0.13];   % gold
c4 = [0.49 0.18 0.56];   % purple
c5 = [0.47 0.67 0.19];   % green
c6 = [0.30 0.75 0.93];   % light blue
c7 = [0.64 0.08 0.18];   % dark red
pal = {c1,c2,c3,c4,c5,c6,c7};

fprintf('=== Generating polished Section 4 figures ===\n\n');

%% =====================================================================
%  FIG 1: Uniform-Gaussian Revenue
%  n=100, rho in {0.5,0.7,0.85,0.95}
%  CLT lines + MC dots, positive-revenue region only
% ======================================================================
fprintf('Fig 1: Uniform-Gaussian revenue\n');
n1 = 100; Nmc1 = 4e5;
rho_list = [0.5 0.7 0.85 0.95];
kgrid = linspace(1, n1, 400);

f1 = figure('Position',[80 80 560 420]); hold on;
legs1 = {}; amax = 0;
for idx = 1:length(rho_list)
    rho = rho_list(idx); sig = sqrt(1-rho^2);
    vh = var_h_exact(rho);
    Sig_e = sqrt(n1*vh);
    Rclt = 2*Sig_e*normpdf(kgrid/Sig_e) - n1*(1-normcdf(kgrid/Sig_e));
    amax = max(amax, max(Rclt));

    % Monte Carlo benchmark
    xm = 2*rand(Nmc1,n1)-1;
    ym = rho*xm + sig*randn(Nmc1,n1);
    hm = trunc_mean_mat(ym,rho);
    Sem = sum(hm,2); Spm = sum(2*xm-1,2);
    ksub = linspace(5, n1-5, 18);
    Rmc = arrayfun(@(k) mean(Spm.*(Sem>=k)), ksub);

    plot(kgrid, Rclt, '-', 'Color', pal{idx}, 'LineWidth', 2.5);
    plot(ksub, Rmc, 'o', 'Color', pal{idx}, 'MarkerSize', 5, ...
        'MarkerFaceColor', pal{idx}, 'HandleVisibility', 'off');
    legs1{end+1} = ['$\rho=' num2str(rho) '$'];
end
ytop = amax*1.18;
plot([n1/2 n1/2],[0 ytop],'k--','LineWidth',1.2,'HandleVisibility','off');
text(n1/2+2, 0.93*ytop, '$k_{\max}=n/2$', 'FontSize', 14);
plot([0 n1],[0 0],'k-','LineWidth',0.5,'HandleVisibility','off');
xlabel('Aggregate threshold $k^*$');
ylabel('Expected revenue $R(k^*)$');
legend(legs1,'Location','northeast','FontSize',12,'Box','on');
grid on; box on;
xlim([0 n1]); ylim([-0.015*ytop ytop]);
set(gca,'GridAlpha',0.15);
exportgraphics(f1,fullfile(outdir,'fig1_uniform_revenue.png'),'Resolution',400);
close(f1); fprintf('   Done.\n');

%% =====================================================================
%  FIG 2a: Beta-Gaussian Dynamic Boundary
% ======================================================================
fprintf('Fig 2a: Beta dynamic boundary\n');
ab = 2.5; bb = 3.5; rhob = 0.75; sigb = sqrt(1-rhob^2);
yg2 = linspace(-2.5, 3.5, 250);
xh2 = zeros(size(yg2)); Hh2 = zeros(size(yg2));
for j=1:length(yg2)
    [xh2(j),Hh2(j)] = beta_post(yg2(j),rhob,sigb,ab,bb);
end

f2a = figure('Position',[80 80 560 420]); hold on;
plot(Hh2, xh2, '-', 'Color', c1, 'LineWidth', 3);
Hl = linspace(0, max(Hh2)*1.18, 200);
lams2 = [0.3 1.0 3.0]; lstyles = {'--','-.','--'};
bcols = {c2, c3, c4};
for i=1:3
    plot(Hl, lams2(i)/(1+lams2(i))*Hl, lstyles{i}, ...
        'Color', bcols{i}, 'LineWidth', 2);
end
xlabel('$\hat{H}(y_i)$'); ylabel('$\hat{x}(y_i)$');
legend({'Posterior locus','$\lambda^*\!=\!0.3$','$\lambda^*\!=\!1.0$',...
    '$\lambda^*\!=\!3.0$'},'Location','southeast','FontSize',11,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([0 max(Hl)]); ylim([0 max(xh2)*1.18]);
exportgraphics(f2a,fullfile(outdir,'fig2a_beta_boundary.png'),'Resolution',400);
close(f2a); fprintf('   Done.\n');

%% =====================================================================
%  FIG 2b: Jensen Gap
% ======================================================================
fprintf('Fig 2b: Jensen gap\n');
Hofx = zeros(size(yg2));
for j=1:length(yg2)
    xc = min(max(xh2(j),1e-6),1-1e-6);
    Hofx(j) = beta_ih(xc,ab,bb);
end
gap2 = max(Hh2 - Hofx, 0);

f2b = figure('Position',[80 80 560 420]); hold on;
fill([yg2 fliplr(yg2)], [gap2 zeros(size(gap2))], c5, ...
    'FaceAlpha', 0.35, 'EdgeColor', 'none');
plot(yg2, Hh2, '-', 'Color', c1, 'LineWidth', 2.5);
plot(yg2, Hofx, '--', 'Color', c2, 'LineWidth', 2.5);
xlabel('Signal $y_i$'); ylabel('Inverse hazard rate');
legend({'Jensen gap','$E[H(x_i) \mid y_i]$','$H(E[x_i \mid y_i])$'},...
    'Location','northeast','FontSize',11,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([yg2(1) yg2(end)]); ylim([0 max(Hh2)*1.1]);
exportgraphics(f2b,fullfile(outdir,'fig2b_jensen_gap.png'),'Resolution',400);
close(f2b); fprintf('   Done.\n');

%% =====================================================================
%  FIG 3a: Clamped Posterior Plateauing
% ======================================================================
fprintf('Fig 3a: Clamped posterior\n');
Ccl = 1.5; scl = 1.0; psig = 1.0;
ygc = linspace(-8,8,500);
Flo = normcdf(-Ccl,0,psig); Fhi = 1-normcdf(Ccl,0,psig);
[hxc, hpc] = clamp_post(ygc,Ccl,scl,psig,Flo,Fhi);

f3a = figure('Position',[80 80 560 420]); hold on;
% Pooling shading
fill([ 5.5  8  8  5.5],[-2.1 -2.1 2.1 2.1],[.91 .91 .91],...
    'EdgeColor','none','HandleVisibility','off');
fill([-8 -5.5 -5.5 -8],[-2.1 -2.1 2.1 2.1],[.91 .91 .91],...
    'EdgeColor','none','HandleVisibility','off');
plot(ygc, hxc, '-', 'Color', c1, 'LineWidth', 2.8);
plot([-8 8],[ Ccl  Ccl],'--','Color',c2,'LineWidth',1.5);
plot([-8 8],[-Ccl -Ccl],'--','Color',c2,'LineWidth',1.5);
text( 6.75, 0, 'Pooling','HorizontalAlignment','center','FontSize',13,...
    'FontWeight','bold','Color',[.4 .4 .4]);
text(-6.75, 0, 'Pooling','HorizontalAlignment','center','FontSize',13,...
    'FontWeight','bold','Color',[.4 .4 .4]);
text(-7.4, Ccl+0.14,'$+C$','Color',c2,'FontSize',13);
text(-7.4,-Ccl-0.26,'$-C$','Color',c2,'FontSize',13);
xlabel('Signal $y_i$'); ylabel('Posterior mean $h_x(y_i)$');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([-8 8]); ylim([-2.1 2.1]);
exportgraphics(f3a,fullfile(outdir,'fig3a_clamped_posterior.png'),'Resolution',400);
close(f3a); fprintf('   Done.\n');

%% =====================================================================
%  FIG 3b: Clamped Score Boundedness
% ======================================================================
fprintf('Fig 3b: Clamped scores\n');
f3b = figure('Position',[80 80 560 420]); hold on;
lc = [0.0 0.3 0.8]; lg3 = {};
scols = {c1, c2, c3};
for i=1:length(lc)
    Sc = hxc + lc(i)*hpc;
    plot(ygc,Sc,'-','Color',scols{i},'LineWidth',2.5);
    % Dotted asymptotic limits
    plot([-8 8],[max(Sc) max(Sc)],':','Color',scols{i},...
        'LineWidth',1.2,'HandleVisibility','off');
    plot([-8 8],[min(Sc) min(Sc)],':','Color',scols{i},...
        'LineWidth',1.2,'HandleVisibility','off');
    lg3{end+1} = ['$\lambda=' num2str(lc(i)) '$'];
end
xlabel('Signal $y_i$'); ylabel('Score $S_i(y_i,\lambda)$');
legend(lg3,'Location','southeast','FontSize',12,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([-8 8]);
% Tight y-limits based on data
allSc = hxc + 0.8*hpc;
ylim([min(allSc)*1.15 max(allSc)*1.15]);
exportgraphics(f3b,fullfile(outdir,'fig3b_clamped_scores.png'),'Resolution',400);
close(f3b); fprintf('   Done.\n');

%% =====================================================================
%  FIG 4a: Nonlinear Gaussian Scores
% ======================================================================
fprintf('Fig 4a: Nonlinear scores\n');
snl = 0.4; anl = 2.5;
ynl = linspace(-2.5,2.5,300);
hxnl = zeros(size(ynl)); hpnl = zeros(size(ynl));
for j=1:length(ynl)
    yv = ynl(j);
    kf = @(x) 0.5*normpdf(yv,tanh(anl*x),snl);
    d = integral(kf,-1,1,'RelTol',1e-10);
    hxnl(j) = integral(@(x) x.*kf(x),-1,1,'RelTol',1e-10)/d;
    hpnl(j) = integral(@(x) (2*x-1).*kf(x),-1,1,'RelTol',1e-10)/d;
end

f4a = figure('Position',[80 80 560 420]); hold on;
lnl = [0.0 0.5 1.0 2.0]; lg4 = {};
ncols = {c1,c2,c3,c4};
for i=1:length(lnl)
    plot(ynl, hxnl+lnl(i)*hpnl, '-', 'Color', ncols{i}, 'LineWidth', 2.5);
    lg4{end+1} = ['$\lambda^*=' num2str(lnl(i)) '$'];
end
plot([-2.5 2.5],[0 0],'k-','LineWidth',0.6,'HandleVisibility','off');
xlabel('Signal $y_i$'); ylabel('Score $S_i(y_i,\lambda^*)$');
legend(lg4,'Location','northwest','FontSize',11,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
exportgraphics(f4a,fullfile(outdir,'fig4a_nonlinear_scores.png'),'Resolution',400);
close(f4a); fprintf('   Done.\n');

%% =====================================================================
%  FIG 4b: Pivotal Threshold Map
% ======================================================================
fprintf('Fig 4b: Pivotal threshold\n');
ls = 1.0; Ss = hxnl + ls*hpnl;
[Su,iu] = unique(Ss); yu = ynl(iu);
Kr = linspace(Su(2),Su(end-1),300);
ysm = interp1(Su,yu,Kr,'pchip');

f4b = figure('Position',[80 80 560 420]);
plot(Kr,ysm,'-','Color',c4,'LineWidth',2.8);
xlabel('Hurdle $K_{-i}(\lambda)$'); ylabel('Pivotal threshold $y_i^*$');
grid on; box on; set(gca,'GridAlpha',0.15);
exportgraphics(f4b,fullfile(outdir,'fig4b_pivotal_threshold.png'),'Resolution',400);
close(f4b); fprintf('   Done.\n');

%% =====================================================================
%  FIG 5: Critical Noise Scaling (bisect in rho)
% ======================================================================
fprintf('Fig 5: Critical noise scaling\n');
nlist5 = [10 20 50 100 200];
Rtgt = 0.05; Nmc5 = 2e5;
rho_crit = zeros(size(nlist5));

for in=1:length(nlist5)
    nn = nlist5(in); fprintf('   n=%d',nn);
    rlo = 0.001; rhi = 0.999;
    for it=1:50
        rm = (rlo+rhi)/2; sm = sqrt(1-rm^2);
        xm = 2*rand(Nmc5,nn)-1;
        ym = rm*xm + sm*randn(Nmc5,nn);
        hm = trunc_mean_mat(ym,rm);
        Rm = mean(sum(2*xm-1,2).*(sum(hm,2)>=nn/2));
        if Rm > Rtgt, rhi = rm; else, rlo = rm; end
    end
    rho_crit(in) = (rlo+rhi)/2;
    fprintf(' -> sig=%.4f\n', sqrt(1-rho_crit(in)^2));
end
sc = sqrt(1 - rho_crit.^2);
Cf = sc(end)*sqrt(nlist5(end))/sqrt(log(1/Rtgt));
sthy = Cf*sqrt(log(1/Rtgt)./nlist5);
Crf = sc(1)*sqrt(nlist5(1));

f5 = figure('Position',[80 80 560 420]); hold on;
plot(nlist5, sc, 's-', 'Color', c1, 'MarkerSize', 10, ...
    'MarkerFaceColor', c1, 'LineWidth', 2.5);
plot(nlist5, sthy, '--', 'Color', c2, 'LineWidth', 2.5);
plot(nlist5, Crf./sqrt(nlist5), ':', 'Color', [.5 .5 .5], 'LineWidth', 2);
set(gca,'XScale','log','YScale','log');
xlabel('Population size $n$'); ylabel('Critical noise $\sigma_{\mathrm{crit}}$');
legend({'MC bisection','$C\sqrt{\ln(1/R)/n}$','$O(n^{-1/2})$'},...
    'Location','southwest','FontSize',11,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
exportgraphics(f5,fullfile(outdir,'fig5_noise_scaling.png'),'Resolution',400);
close(f5); fprintf('   Done.\n');

%% =====================================================================
%  FIG 6: Gaussian vs Laplace MLRP
% ======================================================================
fprintf('Fig 6: MLRP comparison\n');
x1r=-0.8; x2r=0.8; bl=1.0; sg=1.2;
ylr = linspace(-5,5,500);
Ll = (abs(ylr-x1r)-abs(ylr-x2r))/bl;
Lg = (2*ylr*(x2r-x1r)+x1r^2-x2r^2)/(2*sg^2);

f6 = figure('Position',[80 80 560 420]); hold on;
% Shading for constant-LR regions
fill([-5 x1r x1r -5],[-3 -3 3 3],[.94 .91 .97],...
    'EdgeColor','none','HandleVisibility','off');
fill([x2r 5 5 x2r],[-3 -3 3 3],[.94 .91 .97],...
    'EdgeColor','none','HandleVisibility','off');
plot(ylr,Lg,'--','Color',c1,'LineWidth',2.5);
plot(ylr,Ll,'-','Color',c4,'LineWidth',2.8);
% Vertical markers at x1, x2
plot([x1r x1r],[-2.5 2.5],':','Color',[.45 .45 .45],...
    'LineWidth',1.2,'HandleVisibility','off');
plot([x2r x2r],[-2.5 2.5],':','Color',[.45 .45 .45],...
    'LineWidth',1.2,'HandleVisibility','off');
text(x1r-0.05,-2.25,'$x_1$','HorizontalAlignment','center','FontSize',14);
text(x2r+0.05,-2.25,'$x_2$','HorizontalAlignment','center','FontSize',14);
text(-3.2, 1.3,'Constant LR','HorizontalAlignment','center',...
    'FontSize',12,'Color',c4,'FontAngle','italic');
text( 3.2,-1.3,'Constant LR','HorizontalAlignment','center',...
    'FontSize',12,'Color',c4,'FontAngle','italic');
xlabel('Signal $y_i$'); ylabel('Log-likelihood ratio');
legend({'Gaussian (strict MLRP)','Laplace (weak MLRP)'},...
    'Location','northwest','FontSize',12,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([-5 5]); ylim([-2.5 2.5]);
exportgraphics(f6,fullfile(outdir,'fig6_laplace_gaussian.png'),'Resolution',400);
close(f6); fprintf('   Done.\n');

%% =====================================================================
%  FIG 7a: Conditional Allocation
% ======================================================================
fprintf('Fig 7a: Conditional allocation\n');
rt = 0.7; st = sqrt(1-rt^2);
ysl = [-0.5 0.0 0.5 1.0];
xt = linspace(0.001,0.999,300);
tcols = {c1,c2,c3,c4};

f7a = figure('Position',[80 80 560 420]); hold on;
lg7 = {};
for i=1:length(ysl)
    plot(xt, normcdf((rt*xt-ysl(i))/st), '-', 'Color', tcols{i}, 'LineWidth', 2.5);
    lg7{end+1} = ['$y_i^*=' num2str(ysl(i)) '$'];
end
xlabel('Type $x_i$'); ylabel('$Q_i(x_i\mid y_{-i})$');
legend(lg7,'Location','southeast','FontSize',11,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([0 1]); ylim([0 1.05]);
exportgraphics(f7a,fullfile(outdir,'fig7a_cond_allocation.png'),'Resolution',400);
close(f7a); fprintf('   Done.\n');

%% =====================================================================
%  FIG 7b: Conditional Transfer
% ======================================================================
fprintf('Fig 7b: Conditional transfer\n');
f7b = figure('Position',[80 80 560 420]); hold on;
for i=1:length(ysl)
    Tc = cond_transfer(xt, ysl(i), rt, st);
    plot(xt, Tc, '-', 'Color', tcols{i}, 'LineWidth', 2.5);
end
xlabel('Type $x_i$'); ylabel('$T_i(x_i\mid y_{-i})$');
legend(lg7,'Location','northwest','FontSize',11,'Box','on');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([0 1]);
exportgraphics(f7b,fullfile(outdir,'fig7b_cond_transfer.png'),'Resolution',400);
close(f7b); fprintf('   Done.\n');

%% =====================================================================
%  FIG 8: Shadow Price R(lambda)
% ======================================================================
fprintf('Fig 8: Shadow price\n');
n8 = 50; rho8 = 0.6;
vh8 = var_h_exact(rho8); Sig8 = sqrt(n8*vh8);
lam_grid = linspace(0.01, 8, 500);
R_lam = zeros(size(lam_grid));
for il = 1:length(lam_grid)
    lam = lam_grid(il);
    k = n8*lam/(1+2*lam);
    R_lam(il) = 2*Sig8*normpdf(k/Sig8) - n8*(1-normcdf(k/Sig8));
end

f8 = figure('Position',[80 80 560 420]); hold on;
plot(lam_grid, R_lam, '-', 'Color', c1, 'LineWidth', 3);

Rts = [0.01 0.05 0.10 0.20];
rcols = {c2,c3,c4,c5};
for i=1:length(Rts)
    Rt = Rts(i);
    plot([0 8.5],[Rt Rt],'--','Color',rcols{i},'LineWidth',1.2,...
        'HandleVisibility','off');
    % Find crossing
    idx_above = find(R_lam >= Rt);
    if ~isempty(idx_above)
        ic = idx_above(end);
        if ic < length(lam_grid)
            ls = interp1(R_lam(ic:ic+1), lam_grid(ic:ic+1), Rt, 'linear');
            plot(ls, Rt, 'o', 'Color', rcols{i}, 'MarkerSize', 8, ...
                'MarkerFaceColor', rcols{i}, 'HandleVisibility', 'off');
            % Place label above the dot, offset to avoid overlap
            text(ls+0.2, Rt+0.012, ...
                sprintf('$\\lambda^* \\approx %.2f$', ls), ...
                'FontSize', 11, 'Color', rcols{i});
        end
    end
    % Revenue label on the right margin
    text(7.5, Rt+0.006, sprintf('$R=%.2f$', Rt), ...
        'FontSize', 11, 'Color', rcols{i});
end

xlabel('Multiplier $\lambda$'); ylabel('Revenue $R(\lambda)$');
grid on; box on; set(gca,'GridAlpha',0.15);
xlim([0 8]); ylim([-0.015 max(R_lam)*1.18]);
exportgraphics(f8,fullfile(outdir,'fig8_shadow_price.png'),'Resolution',400);
close(f8); fprintf('   Done.\n');

fprintf('\n=== All figures saved to ./%s/ ===\n', outdir);

%% ======================== HELPER FUNCTIONS ==============================

function H = trunc_mean_mat(Y, rho)
    sig = sqrt(1-rho^2); s = sig/rho;
    mu = Y/rho; a = (-1-mu)/s; b = (1-mu)/s;
    Z = normcdf(b)-normcdf(a); Z(Z<1e-15) = 1e-15;
    H = mu + s*(normpdf(a)-normpdf(b))./Z;
    H = max(-1,min(1,H));
end

function h = trunc_mean_sc(y, rho)
    sig = sqrt(1-rho^2); s = sig/rho;
    mu = y/rho; a = (-1-mu)/s; b = (1-mu)/s;
    Z = normcdf(b)-normcdf(a);
    if Z<1e-15
        if mu>0, h=1; elseif mu<0, h=-1; else, h=0; end
        return;
    end
    h = mu + s*(normpdf(a)-normpdf(b))/Z;
    h = max(-1,min(1,h));
end

function v = var_h_exact(rho)
    sig = sqrt(1-rho^2);
    my = @(y) 0.5*(normcdf((y+rho)/sig)-normcdf((y-rho)/sig))/rho;
    hf = @(y) trunc_mean_sc(y,rho);
    E1 = integral(@(y) hf(y).*my(y),-10,10,'ArrayValued',true,'RelTol',1e-10);
    E2 = integral(@(y) hf(y).^2.*my(y),-10,10,'ArrayValued',true,'RelTol',1e-10);
    v = E2 - E1^2;
end

function [xh,Hh] = beta_post(y,rho,sig,al,be)
    ep = 1e-7;
    fp = @(x) betapdf(x,al,be).*normpdf(y,rho*x,sig);
    d = integral(fp,ep,1-ep,'RelTol',1e-10,'AbsTol',1e-12);
    xh = integral(@(x) x.*fp(x),ep,1-ep,'RelTol',1e-10,'AbsTol',1e-12)/d;
    Hh = integral(@(x) beta_ih(x,al,be).*fp(x),ep,1-ep,...
        'RelTol',1e-10,'AbsTol',1e-12)/d;
end

function H = beta_ih(x,al,be)
    H = (1-betacdf(x,al,be))./max(betapdf(x,al,be),1e-15);
end

function [hx,hp] = clamp_post(yg,C,sig,psig,Flo,Fhi)
    alo = Flo*normpdf(yg,-C,sig);
    ahi = Fhi*normpdf(yg,C,sig);
    st = sqrt(sig^2+psig^2); sp = sig*psig/st;
    mp = yg*psig^2/st^2;
    atn = (-C-mp)/sp; btn = (C-mp)/sp;
    Ztn = max(normcdf(btn)-normcdf(atn),1e-15);
    intm = normpdf(yg,0,st).*Ztn;
    my = alo+ahi+intm;
    Exi = mp + sp*(normpdf(atn)-normpdf(btn))./Ztn;
    hx = ((-C)*alo + C*ahi + Exi.*intm)./my;

    phif = @(x) x - (1-normcdf(x,0,psig))./max(normpdf(x,0,psig),1e-12);
    pmC = phif(-C); ppC = phif(C);
    hp = zeros(size(yg));
    for j=1:length(yg)
        wlo = alo(j)/my(j); whi = ahi(j)/my(j); win = intm(j)/my(j);
        if win > 1e-12
            mup = mp(j);
            Zloc = max(normcdf((C-mup)/sp)-normcdf((-C-mup)/sp),1e-15);
            ppdf = @(x) normpdf(x,mup,sp)/Zloc;
            Epi = integral(@(x) phif(x).*ppdf(x),-C+1e-6,C-1e-6,...
                'RelTol',1e-8,'AbsTol',1e-10);
        else
            Epi = 0;
        end
        hp(j) = wlo*pmC + whi*ppC + win*Epi;
    end
end

function T = cond_transfer(xg, ys, rho, sigma)
    u = (rho*xg - ys)/sigma;
    Cys = -ys*normcdf(-ys/sigma) + sigma*normpdf(-ys/sigma);
    T = (ys/rho)*normcdf(u) - (sigma/rho)*normpdf(u) + Cys/rho;
end