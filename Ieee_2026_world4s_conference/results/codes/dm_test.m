function result = dm_test(SA, SB, label_A, label_B)
%% dm_test.m  |  FILE 16 of 30
%
%  One-sided Diebold-Mariano (1995) test.
%  H0: E[eA^2] = E[eB^2]
%  H1: E[eB^2] < E[eA^2]   (B beats A)
%  d_t = eA_t^2 - eB_t^2  — positive DM => B improves over A.
%  Long-run variance: Newey-West, Andrews bandwidth h=floor(4*(T/100)^(2/9)).
%  p-value (one-sided): P(Z > DM).
%
%  Inputs:  SA, SB  — result structs (.timestamps_valid, .error_level)
%           label_A, label_B — display strings
%  Output:  result  — struct (DM, pVal, meanD, N, bandwidth, stars)

    [eA, eB, ~] = align_errors(SA, SB);
    valid = ~isnan(eA) & ~isnan(eB);
    eA = eA(valid); eB = eB(valid);

    T     = numel(eA);
    d     = eA.^2 - eB.^2;
    meanD = mean(d);
    h     = floor(4*(T/100)^(2/9));

    dc  = d - meanD;
    g0  = mean(dc.^2);
    lrv = g0;
    for j = 1:h
        lrv = lrv + 2*(1-j/(h+1))*mean(dc(j+1:end).*dc(1:end-j));
    end
    lrv = max(lrv, 1e-12);

    dmStat = meanD / sqrt(lrv/T);
    pVal   = 0.5 * erfc(dmStat/sqrt(2));

    result.label_base  = label_A;
    result.label_model = label_B;
    result.DM          = dmStat;
    result.pVal        = pVal;
    result.meanD       = meanD;
    result.N           = T;
    result.bandwidth   = h;
    result.stars       = significance_stars(pVal);

    fprintf("DM  |  base: %-16s  model: %-16s\n", label_A, label_B);
    fprintf("    N=%d  h=%d  meanD=%+.4f  DM=%+.4f  p=%.4f  %s\n\n", ...
            T, h, meanD, dmStat, pVal, result.stars);
end
