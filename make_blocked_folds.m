function foldID = make_blocked_folds(n, K)
%% make_blocked_folds.m  |  FILE 10 of 30
%  Splits n obs into K contiguous (temporally ordered) folds.
%  Prevents future leakage in time-series cross-validation.

    foldID = zeros(n,1);
    edges  = round(linspace(1, n+1, K+1));
    for k = 1:K
        foldID(edges(k):edges(k+1)-1) = k;
    end
end
