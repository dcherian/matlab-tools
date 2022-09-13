function [] = MakeKozak(hax)

    if ~exist('hax', 'var') | isempty(hax)
        hax = gca;
    end

    hlines = findobj(hax.Children, '-property', 'xdata', '-property', 'ydata');

    xlo = []; xhi = []; ylo = []; yhi = [];
    for ii=1:length(hlines)
        xlo = min([xlo hlines(ii).XData]);
        xhi = max([xhi hlines(ii).XData]);

        ylo = min([ylo hlines(ii).YData]);
        yhi = max([yhi hlines(ii).YData]);
    end

    hax.XTick((hax.XTick < xlo) | (hax.XTick > xhi)) = [];
    hax.YTick((hax.YTick < ylo) | (hax.YTick > yhi)) = [];

    hax.XTick = sort(unique([hax.XTick xlo xhi]));
    hax.YTick = sort(unique([hax.YTick ylo yhi]));

    hax.XAxis.Axle.VertexData(1,:) = single([xlo, xhi]);
    hax.YAxis.Axle.VertexData(2,:) = single([ylo, yhi]);
end