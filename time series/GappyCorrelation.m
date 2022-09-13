%      [c, lags, delay] = GappyCorrelation(var1, var2)
% if not provided with var2 returns auto-correlation
% delay uses signal processing toolbox's finddelay to find *lowest*
% lag with maximum absolute auto-correlation.
function [c, lags, delay] = GappyCorrelation(var1, var2)

    nan1 = isnan(var1);

    if ~exist('var2', 'var')
        var2 = var1;
        nan12 = nan1;
    else
        nan2 = isnan(var2);
        nan12 = nan1 | nan2;
    end
    [gapstart, gapend] = FindGaps(fillnan(double(nan12), 1));

    if gapstart == gapend % no gaps
        [c, lags] = xcorr(in, 'coeff');
    else % use largest continuous record length
        [~, imax] = max(gapstart(2:end)-gapend(1:end-1));
        range = gapend(imax)+1:gapstart(imax+1)-1;

        if isempty(range)
            c = NaN; lags = NaN; delay = NaN;
            return;
        end

        assert(all(isnan(var1(range)) ~= 1));
        assert(all(isnan(var2(range)) ~= 1));

        [c, lags] = xcorr(var1(range), var2(range), 'coeff');
        delay = finddelay(var1(range), var2(range));
    end
