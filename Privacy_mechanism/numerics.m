% =========================================================================
% Numerical Illustrations for "Optimal Mechanism Design under Privacy Constraints"
% Complete MATLAB Script (Including Ex-Post Transfers)
% =========================================================================

clear; clc; close all;

%% =========================================================================
% Global plotting setup
% =========================================================================
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultAxesLineWidth', 1.0);
set(groot, 'defaultLineLineWidth', 2.0);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

rng(42); % For reproducibility

outdir = 'figures_privacy_mechanism';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% MATLAB-safe publication color palette
c1 = [0.0000, 0.4470, 0.7410]; % blue
c2 = [0.8500, 0.3250, 0.0980]; % orange
c3 = [0.9290, 0.6940, 0.1250]; % yellow
c4 = [0.4940, 0.1840, 0.5560]; % purple
c5 = [0.4660, 0.6740, 0.1880]; % green
colors = {c1, c2, c3, c4, c5};

disp('Generating corrected numerical figures for the paper...');

%% =========================================================================
% Figure 1. Uniform-Gaussian revenue curve (Exact Deterministic Variance)
% =========================================================================
disp('-> Figure 1: Uniform-Gaussian revenue curve');
n = 20;
k_vals = linspace(0, n, 250);
rho_vals = [0.3, 0.6, 0.9];

fig1 = figure('Position', [100, 100, 560, 420]);
hold on;
maxRevAll = -inf;

for i = 1:length(rho_vals)
    rho = rho_vals(i);
    s = sqrt(1 - rho^2);
    
    % Exact marginal density m(y) for Uniform[-1,1] prior
    m_y = @(y) (normcdf((y + rho)/s) - normcdf((y - rho)/s)) / (2*rho);
    
    % Exact Variance via Numerical Integration
    integrand = @(y) (uniform_exact_posterior_mean(y, rho).^2) .* m_y(y);
    var_h = integral(integrand, -10, 10, 'ArrayValued', true, 'RelTol', 1e-6);
    
    Sigma_e = sqrt(n * var_h); 
    
    % Revenue Curve
    z = k_vals / Sigma_e;
    Rev = 2 * Sigma_e .* normpdf(z) - n .* (1 - normcdf(z));
    maxRevAll = max(maxRevAll, max(Rev));
    
    plot(k_vals, Rev, 'Color', colors{i}, 'DisplayName', ['$\rho = ', num2str(rho), '$']);
end

plot([n/2, n/2], [0, maxRevAll*1.1], 'k--', 'LineWidth', 1.3, 'HandleVisibility', 'off');
text(n/2 + 0.35, 0.82 * maxRevAll, '$k^* = n/2$', 'Interpreter', 'latex');

xlabel('Aggregate threshold $k^*$');
ylabel('Expected revenue $\mathcal{R}(k^*)$');
title('Revenue--Threshold Relation under Uniform-Gaussian Signals');
legend('Location', 'northeast');
grid on; box on;
xlim([0, n]); ylim([0, 1.08 * maxRevAll]);
exportgraphics(fig1, fullfile(outdir, 'fig1_uniform_gaussian_revenue_exact.png'), 'Resolution', 400);

%% =========================================================================
% Figure 2. Beta-Gaussian posterior objects and dynamic boundary
% =========================================================================
disp('-> Figure 2: Beta-Gaussian dynamic boundary');
alpha_beta = 2.5;
beta_beta  = 3.5;
rho_beta   = 0.75;
sigma_beta = sqrt(1 - rho_beta^2);

y_grid_beta = linspace(-2.5, 3.0, 100);
xhat_beta = zeros(size(y_grid_beta));
Hhat_beta = zeros(size(y_grid_beta));

for j = 1:length(y_grid_beta)
    y = y_grid_beta(j);
    [xhat_beta(j), Hhat_beta(j)] = beta_posterior_moments(y, rho_beta, sigma_beta, alpha_beta, beta_beta);
end

fig2 = figure('Position', [130, 130, 560, 420]);
hold on;

plot(Hhat_beta, xhat_beta, 'Color', c1, 'LineWidth', 2.6, 'DisplayName', 'Posterior locus $(\hat{H}(y),\hat{x}(y))$');

lambda_vals = [0.5, 1.0, 3.0];
Hline = linspace(0, max(Hhat_beta)*1.1, 200);

for i = 1:length(lambda_vals)
    lam = lambda_vals(i);
    slope = lam / (1 + lam);
    plot(Hline, slope * Hline, '--', 'Color', colors{i+1}, ...
        'DisplayName', ['$\sum \hat{x} = \frac{\lambda^*}{1+\lambda^*}\sum \hat{H},\ \lambda^*=', num2str(lam), '$']);
end

xlabel('Posterior inverse-hazard statistic $\hat{H}(y)$');
ylabel('Posterior mean $\hat{x}(y)$');
title('Beta-Gaussian Posterior Objects and Dynamic Boundary');
legend('Location', 'southeast');
grid on; box on;
xlim([0, max(Hline)]); ylim([0, max(xhat_beta)*1.15]);
exportgraphics(fig2, fullfile(outdir, 'fig2_beta_dynamic_boundary_numerical.png'), 'Resolution', 400);

%% =========================================================================
% Figure 3. Clamped Gaussian posterior plateauing / pooling
% =========================================================================
disp('-> Figure 3: Clamped Gaussian pooling');
C = 1.5;
sig = 1.0;
prior_sig = 1.0; 
y_grid = linspace(-6, 6, 320);

F_mC = normcdf(-C, 0, prior_sig);
F_pC = 1 - normcdf(C, 0, prior_sig);

sig_p = (sig * prior_sig) / sqrt(sig^2 + prior_sig^2);
mu_p = y_grid .* (prior_sig^2) / (sig^2 + prior_sig^2);

marginal_term = normpdf(y_grid, 0, sqrt(sig^2 + prior_sig^2));
a = (-C - mu_p) ./ sig_p;
b = ( C - mu_p) ./ sig_p;
Z_trunc = max(normcdf(b) - normcdf(a), 1e-14);

m_y = F_mC .* normpdf(y_grid, -C, sig) + F_pC .* normpdf(y_grid,  C, sig) + marginal_term .* Z_trunc;
E_x_trunc = mu_p + sig_p .* (normpdf(a) - normpdf(b)) ./ Z_trunc;
num = -C .* F_mC .* normpdf(y_grid, -C, sig) + C .* F_pC .* normpdf(y_grid,  C, sig) + E_x_trunc .* marginal_term .* Z_trunc;

hx = num ./ m_y;

fig3 = figure('Position', [160, 160, 560, 420]);
hold on;
plot(y_grid, hx, 'Color', c1, 'LineWidth', 2.4);

plot([-6, 6], [C, C], '--', 'Color', c2, 'LineWidth', 1.3);
plot([-6, 6], [-C, -C], '--', 'Color', c2, 'LineWidth', 1.3);

fill([3.5, 6, 6, 3.5], [-3, -3, 3, 3], 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([-6, -3.5, -3.5, -6], [-3, -3, 3, 3], 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');

text(-5.25, C + 0.18, 'Upper clamp $+C$', 'Color', c2, 'Interpreter', 'latex');
text(-5.25, -C - 0.28, 'Lower clamp $-C$', 'Color', c2, 'Interpreter', 'latex');
text(4.75, 0, 'Pooling', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
text(-4.75, 0, 'Pooling', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

xlabel('Observed privatized signal $y_i$');
ylabel('Posterior estimator $\hat{x}_i(y_i)$');
title('Posterior Plateauing under Clamped Types');
grid on; box on;
xlim([-6, 6]); ylim([-2.5, 2.5]);
exportgraphics(fig3, fullfile(outdir, 'fig3_clamped_pooling.png'), 'Resolution', 400);

%% =========================================================================
% Figure 4A & 4B. Nonlinear Gaussian Transformation
% =========================================================================
disp('-> Figure 4A/B: Nonlinear Gaussian scores & threshold');
sigma_nl = 0.35;
a_nl = 2.5; 
y_grid_nl = linspace(-2.2, 2.2, 260);
lambda_list_nl = [0.0, 0.5, 1.0, 2.0];

hx_nl = zeros(size(y_grid_nl));
hphi_nl = zeros(size(y_grid_nl));

for j = 1:length(y_grid_nl)
    y = y_grid_nl(j);
    denom = integral(@(x) nonlinear_kernel_uniform(x, y, a_nl, sigma_nl), -1, 1, 'ArrayValued', true);
    num_x = integral(@(x) x .* nonlinear_kernel_uniform(x, y, a_nl, sigma_nl), -1, 1, 'ArrayValued', true);
    num_phi = integral(@(x) (2*x - 1) .* nonlinear_kernel_uniform(x, y, a_nl, sigma_nl), -1, 1, 'ArrayValued', true);
    hx_nl(j) = num_x / denom;
    hphi_nl(j) = num_phi / denom;
end

fig4a = figure('Position', [190, 190, 560, 420]);
hold on;
for i = 1:length(lambda_list_nl)
    lam = lambda_list_nl(i);
    S_nl = hx_nl + lam * hphi_nl;
    plot(y_grid_nl, S_nl, 'Color', colors{i}, 'DisplayName', ['$\lambda^* = ', num2str(lam), '$']);
end
xlabel('Observed signal $y_i$');
ylabel('$S_i(y_i,\lambda^*)=\hat{x}_i(y_i)+\lambda^*\hat{\phi}_i(y_i)$');
title(['Posterior Scores under $T(x)=\tanh(', num2str(a_nl), 'x)$']);
legend('Location', 'best');
grid on; box on;
exportgraphics(fig4a, fullfile(outdir, 'fig4a_nonlinear_scores.png'), 'Resolution', 400);

lambda_star_nl = 1.0;
S_star = hx_nl + lambda_star_nl * hphi_nl;
[S_sort, idx_sort] = unique(S_star);
y_sort = y_grid_nl(idx_sort);
K_vals = linspace(min(S_sort) + 1e-4, max(S_sort) - 1e-4, 200);
y_star_vals = interp1(S_sort, y_sort, K_vals, 'linear');

fig4b = figure('Position', [220, 220, 560, 420]);
plot(K_vals, y_star_vals, 'Color', c4);
xlabel('Opponent-induced hurdle $K_{-i}$');
ylabel('Unique pivotal threshold $y_i^*$');
title(['Pivotal Threshold Map ($\lambda^* = ', num2str(lambda_star_nl), '$)']);
grid on; box on;
exportgraphics(fig4b, fullfile(outdir, 'fig4b_nonlinear_pivotal_threshold.png'), 'Resolution', 400);

%% =========================================================================
% Figure 5. Critical noise scaling
% =========================================================================
disp('-> Figure 5: Critical noise scaling');
n_list = round(logspace(2, 5, 60));
R_vals = [0.1, 0.01, 0.001];

fig5 = figure('Position', [250, 250, 560, 420]);
hold on;
for i = 1:length(R_vals)
    R = R_vals(i);
    sigma_crit = sqrt(log(1/R) ./ n_list);
    plot(n_list, sigma_crit, 'Color', colors{i}, 'DisplayName', ['$R = ', num2str(R), '$']);
end

plot(n_list, 2 ./ sqrt(n_list), 'k--', 'LineWidth', 1.3, 'DisplayName', '$\mathcal{O}(n^{-1/2})$ reference');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Population size $n$');
ylabel('Critical noise limit $\sigma_{\mathrm{crit}}$');
title('Asymptotic Scaling of Maximal Permissible Noise');
legend('Location', 'northeast');
grid on; box on;
exportgraphics(fig5, fullfile(outdir, 'fig5_critical_noise_scaling.png'), 'Resolution', 400);

%% =========================================================================
% Figure 6. Breakdown of strict MLRP under Laplace noise
% =========================================================================
disp('-> Figure 6: Laplace breakdown of strict MLRP');
x1 = -1.0; x2 = 1.0; b = 1.0; sig_cmp = 1.5;
y_grid_dp = linspace(-5, 5, 350);

LLR_Laplace = (abs(y_grid_dp - x1) - abs(y_grid_dp - x2)) / b;
LLR_Gaussian = (-(y_grid_dp - x2).^2 + (y_grid_dp - x1).^2) / (2 * sig_cmp^2);

fig6 = figure('Position', [280, 280, 560, 420]);
hold on;
plot(y_grid_dp, LLR_Gaussian, '--', 'Color', c1, 'LineWidth', 2.0, 'DisplayName', 'Gaussian channel');
plot(y_grid_dp, LLR_Laplace, '-', 'Color', c4, 'LineWidth', 2.4, 'DisplayName', 'Laplace channel (pure $\epsilon$-DP)');

plot([x1, x1], [-4, 4], 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot([x2, x2], [-4, 4], 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');

fill([-5, x1, x1, -5], [-3, -3, 3, 3], c4, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([x2, 5, 5, x2], [-3, -3, 3, 3], c4, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');

text(x1 - 0.05, -2.2, '$x_1$', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
text(x2 + 0.05, -2.2, '$x_2$', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
text(-3.05, 1.35, {'Zero marginal', 'information'}, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
text(3.05, -1.35, {'Zero marginal', 'information'}, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

xlabel('Observed privatized signal $y_i$');
ylabel('Log-likelihood ratio $\ln \frac{k(y_i\mid x_2)}{k(y_i\mid x_1)}$');
title('Breakdown of Strict MLRP under Laplace Noise');
legend('Location', 'northwest');
grid on; box on;
xlim([-5, 5]); ylim([-2.5, 2.5]);
exportgraphics(fig6, fullfile(outdir, 'fig6_laplace_mlrp.png'), 'Resolution', 400);

%% =========================================================================
% Figure 7. Exact Ex-Post Transfers (Uniform-Gaussian)
% =========================================================================
disp('-> Figure 7: Ex-post transfer mechanisms');
n_agents = 20;
k_star = 6.0; % Fixed feasible threshold for illustration
budget_factor = (1 - 2*k_star/n_agents); % Scaling factor

opp_sums = linspace(0, 10, 300); % Sum of opponents: S_{-i}
h_i_vals = [-0.5, 0.0, 0.5];     % Different posterior signals for agent i

fig7 = figure('Position', [310, 310, 560, 420]);
hold on;

for i = 1:length(h_i_vals)
    h_i = h_i_vals(i);
    
    % Allocation condition: S_{-i} + h_i >= k* =>  S_{-i} >= k* - h_i
    allocation_idx = opp_sums >= (k_star - h_i);
    
    % Payment rule: t_i = scaling * (k* - S_{-i}) if allocated, else 0
    t_i = zeros(size(opp_sums));
    t_i(allocation_idx) = budget_factor * (k_star - opp_sums(allocation_idx));
    
    plot(opp_sums, t_i, 'Color', colors{i}, 'DisplayName', ['$\hat{x}_i(y_i) = ', num2str(h_i), '$']);
end

yline(0, 'k-', 'HandleVisibility', 'off');
plot([k_star, k_star], [-2.5, 2.5], 'k:', 'HandleVisibility', 'off');
text(k_star + 0.1, 1.5, '$S_{-i} = k^*$', 'Interpreter', 'latex');

xlabel('Aggregate Opponent Estimator $S_{-i} = \sum_{j \neq i} \hat{x}_j(y_j)$');
ylabel('Ex-Post Payment $t_i(y)$');
title('Pivotal Ex-Post Transfer Rule (Uniform-Gaussian)');
legend('Location', 'southwest');
grid on; box on;
xlim([0, 10]); ylim([-2.5, 2.5]);
exportgraphics(fig7, fullfile(outdir, 'fig7_ex_post_transfers.png'), 'Resolution', 400);

disp('Done. All figures successfully generated!');

%% =========================================================================
% Local helper functions
% =========================================================================
function h = uniform_exact_posterior_mean(y, rho)
    s = sqrt(1 - rho^2) / rho;
    mu = y / rho;
    alpha = (-1 - mu) ./ s;
    beta  = ( 1 - mu) ./ s;
    numer = normpdf(alpha) - normpdf(beta);
    denom = max(normcdf(beta) - normcdf(alpha), 1e-14);
    h = mu + s .* (numer ./ denom);
end

function val = nonlinear_kernel_uniform(x, y, a, sigma)
    T = tanh(a * x);
    val = 0.5 * normpdf(y, T, sigma);
end

function [xhat, Hhat] = beta_posterior_moments(y, rho, sigma, alpha, betaPar)
    epsx = 1e-6; 
    fpost = @(x) beta_prior_pdf(x, alpha, betaPar) .* normpdf(y, rho*x, sigma);
    denom = integral(fpost, epsx, 1 - epsx, 'ArrayValued', true, 'RelTol', 1e-8, 'AbsTol', 1e-10);
    num_x = integral(@(x) x .* fpost(x), epsx, 1 - epsx, 'ArrayValued', true, 'RelTol', 1e-8, 'AbsTol', 1e-10);
    num_H = integral(@(x) beta_inverse_hazard(x, alpha, betaPar) .* fpost(x), epsx, 1 - epsx, 'ArrayValued', true, 'RelTol', 1e-8, 'AbsTol', 1e-10);
    xhat = num_x / denom;
    Hhat = num_H / denom;
end

function f = beta_prior_pdf(x, alpha, betaPar)
    B = beta(alpha, betaPar);
    f = (x.^(alpha - 1) .* (1 - x).^(betaPar - 1)) ./ B;
end

function H = beta_inverse_hazard(x, alpha, betaPar)
    f = beta_prior_pdf(x, alpha, betaPar);
    F = betainc(x, alpha, betaPar);
    H = (1 - F) ./ f;
end