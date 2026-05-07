function save_fig(fig, basepath)
%% ============================================================
%  save_fig.m  (robust version)
%
%  Saves a figure as PNG (300 DPI) and PDF.
%  Uses print() instead of exportgraphics() for compatibility.
%
%  Inputs:
%    fig      : figure handle
%    basepath : full path without extension
%% ============================================================

    %% Make sure figure is valid before saving
    if ~ishandle(fig) || ~isvalid(fig)
        fprintf("  Warning: figure handle invalid, skipping save: %s\n", basepath);
        return;
    end

    %% Make figure the current figure
    figure(fig);

    %% PNG — 300 DPI raster
    png_path = char(basepath + ".png");
    try
        print(fig, png_path, '-dpng', '-r300');
        fprintf("  Saved: %s.png\n", basepath);
    catch ME
        fprintf("  Warning: PNG save failed for %s: %s\n", basepath, ME.message);
    end

    %% PDF — vector
    pdf_path = char(basepath + ".pdf");
    try
        print(fig, pdf_path, '-dpdf', '-bestfit');
        fprintf("  Saved: %s.pdf\n", basepath);
    catch ME
        fprintf("  Warning: PDF save failed for %s: %s\n", basepath, ME.message);
    end

end