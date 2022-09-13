function [] = ggplot()

    hax = gca;

    grid on;
    hax.Color = [1 1 1]*0.90;
    hax.GridColor = [1 1 1];
    hax.GridAlpha = 0.65;
    hax.TickLength = [0 0];
    hax.LineWidth = 2;

    hax.XAxis.Axle.LineStyle = 'none';
    hax.YAxis.Axle.LineStyle = 'none';
    hax.XAxis.TickLabelGapOffset = 0;
    hax.XAxis.TickLabelGapMultiplier = 0;
    hax.YAxis.TickLabelGapOffset = 0;
    hax.YAxis.TickLabelGapMultiplier = 0;

    set(gcf, 'renderer', 'opengl');
end