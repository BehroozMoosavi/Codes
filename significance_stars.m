function s = significance_stars(p)
%% significance_stars.m  |  FILE 13 of 30
%  Returns *** / ** / * / "" for a given p-value.
    if isnan(p),      s = "";
    elseif p < 0.01,  s = "***";
    elseif p < 0.05,  s = "**";
    elseif p < 0.10,  s = "*";
    else,             s = "";
    end
end
