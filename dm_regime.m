function [res_stable, res_ramp, res_full] = dm_regime(SA, SB, label_A, label_B)
%% dm_regime.m  |  FILE 17 of 30
%
%  Regime-specific DM tests (full / stable 10-14 / ramp 16-21).
%  Tests H2: VA-REN advantage is concentrated in ramp hours.
%
%  Inputs:  SA, SB, label_A, label_B  (same as dm_test.m)
%  Outputs: res_stable, res_ramp, res_full  — DM result structs

    [eA, eB, tCommon] = align_errors(SA, SB);
    h_vec     = hour(tCommon);
    stableIdx = h_vec >= 10 & h_vec <= 14;
    rampIdx   = h_vec >= 16 & h_vec <= 21;

    fprintf("=== Regime DM: %s vs %s ===\n", label_A, label_B);
    fprintf("  Full N=%d  Stable N=%d  Ramp N=%d\n\n", ...
            numel(eA), sum(stableIdx), sum(rampIdx));

    mk = @(e,t) struct("timestamps_valid",t,"error_level",e);

    fprintf("--- Full ---\n");
    res_full   = dm_test(mk(eA,tCommon), mk(eB,tCommon), label_A, label_B);

    fprintf("--- Stable (10-14) ---\n");
    res_stable = dm_test(mk(eA(stableIdx),tCommon(stableIdx)), ...
                         mk(eB(stableIdx),tCommon(stableIdx)), ...
                         label_A+" [stable]", label_B+" [stable]");

    fprintf("--- Ramp (16-21) ---\n");
    res_ramp   = dm_test(mk(eA(rampIdx),tCommon(rampIdx)), ...
                         mk(eB(rampIdx),tCommon(rampIdx)), ...
                         label_A+" [ramp]", label_B+" [ramp]");

    fprintf("Summary: DM_full=%+.4f  DM_stable=%+.4f  DM_ramp=%+.4f\n", ...
            res_full.DM, res_stable.DM, res_ramp.DM);
    if res_ramp.DM > res_stable.DM
        fprintf(">> DM_ramp > DM_stable: consistent with H2.\n\n");
    else
        fprintf(">> DM_ramp <= DM_stable: H2 not supported.\n\n");
    end
end
