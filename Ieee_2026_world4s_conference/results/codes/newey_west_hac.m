function V = newey_west_hac(X, resid, h)
%% newey_west_hac.m  |  FILE 09 of 30
%  Newey-West HAC covariance matrix.
%  V = (X'X)^{-1} * Omega * (X'X)^{-1}
%  Bartlett kernel, bandwidth h (Andrews rule: floor(4*(T/100)^(2/9))).
%  Usage: V = newey_west_hac(X, y-X*beta, h);  se = sqrt(diag(V));

    T      = size(X,1);
    scores = X .* resid;
    bread  = inv(X'*X);
    meat   = (scores'*scores)/T;
    for j = 1:h
        w    = 1 - j/(h+1);
        Gj   = (scores(j+1:end,:)'*scores(1:end-j,:))/T;
        meat = meat + w*(Gj+Gj');
    end
    V = bread * meat * bread * T;
end
