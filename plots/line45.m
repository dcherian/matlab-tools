% Draw 45 degree line. takes axis handle as optional input
%     [hline] = line45(hAxis)

function [hline] = line45(hAxis)
    if ~exist('hAxis', 'var'); hAxis = gca; end

    hold all;

    for ii=1:length(hAxis)
        axes(hAxis(ii));
        limx = xlim;
        limy = ylim;
        lim(1) = min([limx limy]);
        lim(2) = max([limx limy]);
        hline = plot(lim, lim, '--', 'Color', [1 1 1]*0.75);
        axis square;
    end
end