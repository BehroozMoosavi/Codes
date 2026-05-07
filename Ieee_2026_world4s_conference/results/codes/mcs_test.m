function [mcs_set, elim_table] = mcs_test(errors_cell, labels, alpha_level, n_boot)
%% ============================================================
%  mcs_test.m
%
%  Model Confidence Set — Hansen, Lunde & Nason (2011).
%  T_max elimination rule with bootstrap critical values.
%
%  Inputs:
%    errors_cell : {M x 1} cell of [T x 1] aligned level errors
%    labels      : [M x 1] string array of method names
%    alpha_level : significance level e.g. 0.10
%    n_boot      : bootstrap replications e.g. 1000
%
%  Outputs:
%    mcs_set    : string array of surviving methods
%    elim_table : table with Step, Eliminated, PValue
%% ============================================================

    M = numel(errors_cell);

    %% Build [T x M] squared-loss matrix robustly
    %  Do NOT use {: } indexing — fails in some MATLAB versions
    %  Build column by column instead
    T_len = numel(errors_cell{1});
    L_mat = zeros(T_len, M);

    for col = 1:M
        e = errors_cell{col};
        if numel(e) ~= T_len
            error("mcs_test: errors_cell{%d} has length %d, expected %d", ...
                  col, numel(e), T_len);
        end
        L_mat(:, col) = e(:).^2;
    end

    %% Keep only rows where ALL models have valid (non-NaN) loss
    valid = all(~isnan(L_mat), 2);
    L_mat = L_mat(valid, :);
    T     = size(L_mat, 1);

    fprintf("MCS: M=%d  T=%d  alpha=%.2f  n_boot=%d\n\n", ...
            M, T, alpha_level, n_boot);

    rng(42);

    in_set = true(M, 1);
    Steps  = zeros(0,1);
    Elim   = strings(0,1);
    Pvals  = zeros(0,1);

    while sum(in_set) > 1

        idx  = find(in_set);
        m    = numel(idx);
        h_bw = floor(4*(T/100)^(2/9));

        %% Pairwise HAC t-statistics
        d_bar = zeros(m, m);
        t_ij  = zeros(m, m);

        for ii = 1:m
            for jj = 1:m
                if ii == jj, continue; end

                d     = L_mat(:, idx(ii)) - L_mat(:, idx(jj));
                dmean = mean(d);
                dc    = d - dmean;

                g0  = mean(dc.^2);
                lrv = g0;
                for lag = 1:h_bw
                    wt  = 1 - lag/(h_bw+1);
                    gj  = mean(dc(lag+1:end) .* dc(1:end-lag));
                    lrv = lrv + 2*wt*gj;
                end
                lrv = max(lrv, 1e-12);

                d_bar(ii,jj) = dmean;
                t_ij(ii,jj)  = dmean / sqrt(lrv/T);
            end
        end

        %% Observed T_max
        T_obs = max(t_ij(:));

        %% Bootstrap null distribution
        T_boot = zeros(n_boot, 1);

        for b = 1:n_boot
            bi   = randi(T, T, 1);
            Lb   = L_mat(bi, idx);
            tb   = zeros(m, m);

            for ii = 1:m
                for jj = 1:m
                    if ii == jj, continue; end

                    db   = Lb(:,ii) - Lb(:,jj);
                    dbc  = db - d_bar(ii,jj);
                    g0b  = mean(dbc.^2);
                    lrvb = g0b;
                    for lag = 1:h_bw
                        wt   = 1 - lag/(h_bw+1);
                        gjb  = mean(dbc(lag+1:end) .* dbc(1:end-lag));
                        lrvb = lrvb + 2*wt*gjb;
                    end
                    lrvb = max(lrvb, 1e-12);
                    tb(ii,jj) = mean(db) / sqrt(lrvb/T);
                end
            end

            T_boot(b) = max(tb(:));
        end

        pv   = mean(T_boot >= T_obs);
        step = numel(Steps) + 1;

        if pv <= alpha_level
            %% Eliminate model with highest average relative loss
            mr = mean(d_bar, 2);
            [~, wl] = max(mr);
            wg = idx(wl);
            en = labels(wg);
            in_set(wg) = false;

            Steps(end+1) = step;
            Elim(end+1)  = en;
            Pvals(end+1) = pv;

            fprintf("  Step %d: eliminated %-18s  T_max=%.4f  p=%.4f\n", ...
                    step, en, T_obs, pv);
        else
            Steps(end+1) = step;
            Elim(end+1)  = "(stopped)";
            Pvals(end+1) = pv;

            fprintf("  Step %d: H0 not rejected  T_max=%.4f  p=%.4f  => stop\n", ...
                    step, T_obs, pv);
            break;
        end

    end

    mcs_set = labels(in_set);

    fprintf("\nMCS (alpha=%.2f) contains: %s\n\n", ...
            alpha_level, strjoin(mcs_set, ", "));

    elim_table = table(Steps(:), Elim(:), Pvals(:), ...
        "VariableNames", ["Step","Eliminated","PValue"]);

end