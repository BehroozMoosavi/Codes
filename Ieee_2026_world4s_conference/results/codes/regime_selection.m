function result = regime_selection(varen_file, predictor_names, n_boot)
%% ============================================================
%  regime_selection.m  (fixed)
%
%  Regime-conditional predictor selection frequency analysis.
%  Fixed: table constructor uses char cell array for VariableNames,
%  and all column vectors are explicitly forced to [p x 1].
%% ============================================================

    if nargin < 3, n_boot = 500; end

    S = load(varen_file);

    if ~isfield(S, 'selection_matrix_valid')
        error("selection_matrix_valid not found in %s.\nRe-run va_ren.m first.", varen_file);
    end

    sel_mat  = double(S.selection_matrix_valid);
    ts_valid = S.timestamps_valid;
    h_vec    = hour(ts_valid);
    [T, p]   = size(sel_mat);

    fprintf("Regime selection: T=%d  p=%d\n", T, p);

    if numel(predictor_names) ~= p
        error("predictor_names has %d entries but selection matrix has %d columns.", ...
              numel(predictor_names), p);
    end

    stableIdx = h_vec >= 10 & h_vec <= 14;
    rampIdx   = h_vec >= 16 & h_vec <= 21;
    otherIdx  = ~stableIdx & ~rampIdx;
    N_s = sum(stableIdx); N_r = sum(rampIdx);

    fprintf("  Stable N=%d  Ramp N=%d  Other N=%d\n\n", N_s, N_r, sum(otherIdx));

    %% Selection frequencies — force [p x 1] column vectors
    freq_full   = reshape(mean(sel_mat,              1), p, 1);
    freq_stable = reshape(mean(sel_mat(stableIdx,:), 1), p, 1);
    freq_ramp   = reshape(mean(sel_mat(rampIdx,  :), 1), p, 1);
    freq_other  = reshape(mean(sel_mat(otherIdx, :), 1), p, 1);
    delta       = freq_ramp - freq_stable;

    %% Bootstrap CIs
    fprintf("  Bootstrapping CIs (n_boot=%d)...\n", n_boot);
    rng(42);
    S_rows = sel_mat(stableIdx, :);
    R_rows = sel_mat(rampIdx,   :);
    delta_boot = zeros(n_boot, p);

    for b = 1:n_boot
        f_s = mean(S_rows(randi(N_s, N_s, 1), :), 1);
        f_r = mean(R_rows(randi(N_r, N_r, 1), :), 1);
        delta_boot(b, :) = f_r - f_s;
    end

    ci_lo = reshape(prctile(delta_boot,  2.5, 1), p, 1);
    ci_hi = reshape(prctile(delta_boot, 97.5, 1), p, 1);
    sig_r = ci_lo > 0;
    sig_s = ci_hi < 0;

    fprintf("  Sig ramp=%d  Sig stable=%d\n\n", sum(sig_r), sum(sig_s));

    %% Predictor groups
    dn  = reshape(replace(predictor_names, "Z_", ""), p, 1);
    grp = repmat("Other", p, 1);

    for j = 1:p
        nm = predictor_names(j);
        if     contains(nm, "LogDeltaLag"),     grp(j) = "Lag log-change";
        elseif contains(nm, "LogLoadLag"),      grp(j) = "Log-level anchor";
        elseif contains(nm, "HourDummy"),       grp(j) = "Hour dummy";
        elseif contains(nm, "WeekdayDummy"),    grp(j) = "Weekday dummy";
        elseif contains(nm, "MonthDummy"),      grp(j) = "Month dummy";
        end
    end

    %% Build full table using 'VariableNames' char syntax (compatible)
    T_full = table( ...
        dn, grp, freq_full, freq_stable, freq_ramp, freq_other, ...
        delta, ci_lo, ci_hi, sig_r, sig_s, ...
        'VariableNames', { ...
        'Predictor','Group','FreqFull','FreqStable','FreqRamp', ...
        'FreqOther','Delta','CI_lo','CI_hi','SigRamp','SigStable'});

    %% Top-10 tables
    T_sr = sortrows(T_full, 'Delta', 'descend');
    T_ss = sortrows(T_full, 'Delta', 'ascend');
    top10_ramp   = T_sr(1:min(10,height(T_sr)), :);
    top10_stable = T_ss(1:min(10,height(T_ss)), :);

    %% Group summary
    gu = unique(grp, 'stable');
    ng = numel(gu);
    GN = strings(ng,1); GS = zeros(ng,1);
    GF = zeros(ng,1); Gs = zeros(ng,1); Gr = zeros(ng,1);
    Gd = zeros(ng,1); Nr = zeros(ng,1); Ns_ = zeros(ng,1);

    for g = 1:ng
        gi = grp == gu(g);
        GN(g)  = gu(g);       GS(g)  = sum(gi);
        GF(g)  = mean(freq_full(gi));
        Gs(g)  = mean(freq_stable(gi));
        Gr(g)  = mean(freq_ramp(gi));
        Gd(g)  = mean(delta(gi));
        Nr(g)  = sum(sig_r(gi));
        Ns_(g) = sum(sig_s(gi));
    end

    T_group = table(GN, GS, GF, Gs, Gr, Gd, Nr, Ns_, ...
        'VariableNames', {'Group','Size','FreqFull','FreqStable', ...
        'FreqRamp','MeanDelta','NSig_Ramp','NSig_Stable'});
    T_group = sortrows(T_group, 'MeanDelta', 'descend');

    %% Hourly selection per group [24 x ng]
    sel_hour = NaN(24, ng);
    for hr = 0:23
        hi = h_vec == hr;
        if ~any(hi), continue; end
        for g = 1:ng
            gi = grp == gu(g);
            sel_hour(hr+1, g) = mean(mean(sel_mat(hi, gi), 2), 'omitnan');
        end
    end

    %% Print
    fprintf("Top 10 RAMP-preferred:\n");
    for i = 1:height(top10_ramp)
        r = top10_ramp(i,:);
        st = ""; if r.SigRamp, st = " *"; end
        fprintf("  %-28s  %-18s  ramp=%.3f  stable=%.3f  D=%+.3f%s\n", ...
            r.Predictor, r.Group, r.FreqRamp, r.FreqStable, r.Delta, st);
    end

    fprintf("\nTop 10 STABLE-preferred:\n");
    for i = 1:height(top10_stable)
        r = top10_stable(i,:);
        st = ""; if r.SigStable, st = " *"; end
        fprintf("  %-28s  %-18s  stable=%.3f  ramp=%.3f  D=%+.3f%s\n", ...
            r.Predictor, r.Group, r.FreqStable, r.FreqRamp, r.Delta, st);
    end

    fprintf("\nGroup summary:\n");
    disp(T_group);

    %% Pack output
    result.T_full              = T_full;
    result.T_group             = T_group;
    result.top10_ramp          = top10_ramp;
    result.top10_stable        = top10_stable;
    result.freq_full           = freq_full;
    result.freq_stable         = freq_stable;
    result.freq_ramp           = freq_ramp;
    result.delta               = delta;
    result.ci_lo               = ci_lo;
    result.ci_hi               = ci_hi;
    result.sig_ramp_flag       = sig_r;
    result.sig_stable_flag     = sig_s;
    result.group               = grp;
    result.group_names_ordered = T_group.Group;
    result.predictor_names     = dn;
    result.N_stable            = N_s;
    result.N_ramp              = N_r;
    result.stableIdx           = stableIdx;
    result.rampIdx             = rampIdx;
    result.delta_boot          = delta_boot;
    result.sel_by_hour_group   = sel_hour;
    result.hours_axis          = (0:23)';

end