function [c, lags] = plotcorr(in, hax)

    [c, lags] = GappyCorrelation(in);

    if ~exist('hax', 'var') | isempty(hax)
        figure;
    else
        axes(hax);
        hold on;
    end
    plot(lags, c);
    linex(0); liney(0);
    xlabel('Lags'); ylabel('Autocorrelation');