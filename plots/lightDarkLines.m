function [] = lightDarkLines(n, hax)
    if ~exist('hax', 'var'), hax = gca; end
    if ~exist('n', 'var'), n = 5; end

    set(hax, 'ColorOrder', ...
             brighten(cbrewer('seq','Reds',n), -0.5));

    % this makes sure that ColorOrder isn't reset when you call
    % plot after using this script.
    hold on;
end