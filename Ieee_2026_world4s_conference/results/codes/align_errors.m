function [eA, eB, tCommon] = align_errors(SA, SB)
%% align_errors.m  |  FILE 11 of 30
%  Aligns level forecast errors from two result structs onto
%  their common timestamp intersection.
%  Required before any pairwise statistical test.
%
%  Inputs:  SA, SB  — structs with .timestamps_valid and .error_level
%  Outputs: eA, eB  — aligned [N x 1] error vectors
%           tCommon — [N x 1 datetime] common timestamps

    [tCommon, ia, ib] = intersect(SA.timestamps_valid(:), SB.timestamps_valid(:));
    eA = SA.error_level(ia);
    eB = SB.error_level(ib);
end
