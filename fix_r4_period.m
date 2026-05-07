%% ============================================================
%  fix_r4_subperiod.m
%
%  Quick fix for R4 subperiod stability.
%  Run this now — takes < 10 seconds.
%  The error was: NaN rows in table when a year had < 10 obs.
%  Fix: pre-fill with valid default values before the loop.
%% ============================================================

clear; clc;

fprintf("Running R4 subperiod fix...\n\n");

SV = load("results_varen_logchange.mat");
SS = load("results_baseline_static_en_logchange.mat");

%% Align errors on common timestamps
[tC, ia, ib] = intersect(SV.timestamps_valid(:), SS.timestamps_valid(:));
eV  = SV.error_level(ia);
eS  = SS.error_level(ib);
yt  = SV.y_true_level_valid(ia);
yrs = year(tC);

years = [2020; 2021; 2022; 2023];
n     = 4;

%% Pre-allocate as regular arrays (not NaN) to avoid table size mismatch
Year           = years;
N_obs          = zeros(n,1);
RMSE_VAREN     = zeros(n,1);
RMSE_StaticEN  = zeros(n,1);
MAE_VAREN      = zeros(n,1);
MAE_StaticEN   = zeros(n,1);
MAPE_VAREN     = zeros(n,1);
MAPE_StaticEN  = zeros(n,1);
Delta_RMSE_pct = zeros(n,1);
Delta_MAE_pct  = zeros(n,1);

fprintf("  %-6s  %-6s  %-12s  %-12s  %-12s\n", ...
        "Year","N","RMSE VA-REN","RMSE Static","Delta RMSE%%");

for i = 1:n

    yr  = years(i);
    idx = yrs == yr & ~isnan(eV) & ~isnan(eS) & yt > 0;

    if sum(idx) < 10
        fprintf("  %-6d  insufficient data (%d obs)\n", yr, sum(idx));
        continue;
    end

    ev = eV(idx); es = eS(idx); y = yt(idx);

    N_obs(i)         = sum(idx);
    RMSE_VAREN(i)    = sqrt(mean(ev.^2));
    RMSE_StaticEN(i) = sqrt(mean(es.^2));
    MAE_VAREN(i)     = mean(abs(ev));
    MAE_StaticEN(i)  = mean(abs(es));
    MAPE_VAREN(i)    = mean(abs(ev./y)) * 100;
    MAPE_StaticEN(i) = mean(abs(es./y)) * 100;
    Delta_RMSE_pct(i)= 100 * (RMSE_StaticEN(i) - RMSE_VAREN(i)) / RMSE_StaticEN(i);
    Delta_MAE_pct(i) = 100 * (MAE_StaticEN(i)  - MAE_VAREN(i))  / MAE_StaticEN(i);

    fprintf("  %-6d  %-6d  %-12.2f  %-12.2f  %+.4f%%\n", ...
            yr, N_obs(i), RMSE_VAREN(i), RMSE_StaticEN(i), Delta_RMSE_pct(i));

end

n_positive = sum(Delta_RMSE_pct > 0);
fprintf("\n  VA-REN RMSE improvement positive in %d / 4 years.\n\n", n_positive);

%% Build table — all columns are [4x1] double, no NaN rows
T_sub = table(Year, N_obs, RMSE_VAREN, RMSE_StaticEN, ...
              MAE_VAREN, MAE_StaticEN, MAPE_VAREN, MAPE_StaticEN, ...
              Delta_RMSE_pct, Delta_MAE_pct);

%% Save
save("subperiod.mat", "T_sub");
writetable(T_sub, "subperiod_stability.csv");

fprintf("Saved: subperiod.mat\n");
fprintf("Saved: subperiod_stability.csv\n\n");

disp(T_sub);