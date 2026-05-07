function out = latex_escape(s)
%% latex_escape.m  |  FILE 14 of 30
%  Escapes _ % & for safe use in LaTeX table cells.
    out = string(s);
    out = replace(out, "_", "\\_");
    out = replace(out, "%", "\\%");
    out = replace(out, "&", "\\&");
end
