function result = gw_test(SA, SB, volatility, timestamps_vol, label_A, label_B)
%% gw_test.m  |  FILE 19 of 30
%
%  Giacomini-White (2006) conditional predictive ability test.
%  Regresses d_t = eA_t^2 - eB_t^2  on  V_{t-1}^{norm}.
%  H0: gamma_1 = 0    H1: gamma_1 > 0  (advantage grows with vol)
%  HAC standard errors via newey_west_hac.m.
%
%  Inputs:
%    SA, SB         — result structs
%    volatility     — [N x 1] NormalizedVolatility
%    timestamps_vol — [N x 1 datetime]
%    label_A/B      — display strings
%  Output:
%    result — struct (gamma0/1, se, t_stat, p_val, R2, stars)

    [eA, eB, tCommon] = align_errors(SA, SB);
    [found, ia] = ismember(tCommon, timestamps_vol);

    vol_lag = NaN(numel(tCommon),1);
    for k = 1:numel(tCommon)
        if found(k) && ia(k) > 1
            vol_lag(k) = volatility(ia(k)-1);
        end
    end

    d     = eA.^2 - eB.^2;
    valid = ~isnan(d) & ~isnan(vol_lag);
    d_v   = d(valid); vol_v = vol_lag(valid);
    T     = numel(d_v);

    X_gw  = [ones(T,1), vol_v];
    beta  = X_gw \ d_v;
    resid = d_v - X_gw*beta;

    SS_res = sum(resid.^2);
    SS_tot = sum((d_v-mean(d_v)).^2);
    R2     = 1 - SS_res/max(SS_tot,1e-12);

    h_bw = floor(4*(T/100)^(2/9));
    V    = newey_west_hac(X_gw, resid, h_bw);
    se   = sqrt(diag(V));

    t_stat = beta(2)/se(2);
    p_val  = 1 - normcdf(t_stat);

    result.label_base  = label_A;
    result.label_model = label_B;
    result.gamma0      = beta(1);
    result.gamma1      = beta(2);
    result.se_gamma0   = se(1);
    result.se_gamma1   = se(2);
    result.t_stat      = t_stat;
    result.p_val       = p_val;
    result.N           = T;
    result.bandwidth   = h_bw;
    result.R2          = R2;
    result.stars       = significance_stars(p_val);

    fprintf("GW test  |  %s vs %s\n", label_A, label_B);
    fprintf("  N=%d  h=%d\n", T, h_bw);
    fprintf("  gamma_0=%+.4f (se=%.4f)\n", beta(1), se(1));
    fprintf("  gamma_1=%+.4f (se=%.4f)  t=%+.4f  p=%.4f  %s\n", ...
            beta(2), se(2), t_stat, p_val, result.stars);
    fprintf("  R^2=%.4f\n\n", R2);

    if beta(2)>0 && p_val<0.10
        fprintf("  >> gamma_1>0 and significant: mechanism confirmed.\n\n");
    elseif beta(2)>0
        fprintf("  >> gamma_1>0 but not significant at 10%%.\n\n");
    else
        fprintf("  >> gamma_1<=0: mechanism not confirmed.\n\n");
    end
end
