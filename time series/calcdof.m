% estimate degrees of freedom using integral timescale
%        [dof, IntegralTimeScale] = calcdof(in, filter, cutoff)
% Set filter = 'narrowband' and cutoff to *band frequencies*
% to use that methodology.
function [dof, IntegralTimeScale] = calcdof(in, filter, cutoff)

    if ~exist('filter', 'var')
        filter = 'broadband';
    end

    nans = isnan(in);
    edges = diff(nans);
    gapstart = find(edges == 1) + 1;
    gapend = find(edges == -1);

    if isnan(in(1))
        gapstart(2:end+1) = gapstart(1:end);
        gapstart(1) = 1;
    end
    if isnan(in(end))
        gapend(end+1) = length(in);
    else
        % last gap is not at the end of time series.
        gapstart(end+1) = length(in);
        gapend(end+1) = length(in);
    end

    if isempty(gapstart) & isempty(gapend) ...
            & isequal(nans, zeros(size(nans)))
        % input series has no gaps
        gapstart = length(in) + 1;
        gapend = gapstart;
    end

    assert(length(gapstart) == length(gapend), ...
           ['FilterSeries: gapstart and gapend are not same ' ...
            'size.']);

    start = 1;
    for ii=1:length(gapstart)
        range = start:gapstart(ii)-1;
        % set start for next iteration to be end of current gap.
        start = gapend(ii)+1;

        if isempty(range)
            % first element is NaN
            continue;
        end

        assert(all(isnan(in(range)) == 0), 'NaNs in selected range!');

        if strcmpi(filter, 'narrowband')
            if ~exist('cutoff', 'var')
                error(['Provide filter cutoffs for narrowband time ' ...
                       'series']);
            end

            % length of segment
            N = length(range);
            % number of resolvable frequencies in segment
            numfreq = length(min(cutoff):1/N:max(cutoff));
            % amplitude and frequency of sinusoid at each frequency
            % is number of degrees of freedom
            dof(ii) = numfreq*2;
            IntegralTimeScale = NaN;
        else
            % Using coef and biased are the same because I normalize by
            % cmax later. Sarah Gille's notes uses biased, but says to
            % normalize so that autocorrelation is 1 at 0 lag.
            [c,lags] = xcorr(in(range) - mean(in(range)), 'coef');
            [cmax,imax] = max(c);

            if isempty(find(c < 0))
                % if no zero-crossing in auto-correlation, then the
                % segment isn't long enough to be meaningful
                dof(ii) = NaN;
                continue;
            end

            % calculate a bunch of timescales and take maximum
            % From Talley et al., Descriptive Physical Oceanography.
            IntegralTimeScale = max(cumtrapz(c(imax:length(c))))./cmax;
            % hold on; plot(cumtrapz(c(imax:length(c))))

            % Techincally, for sin(2π/T t), IntegralTimeScale = T/2π
            % The zero-crossing for that sinusoid is at T/4 points
            % T/4 = IntegralTimeScale * π/2
            % Sarah Gille's notes do not do this though. Kyun? Pata nahin.
            dof(ii) = floor( length(in(range))/(IntegralTimeScale) );
        end
    end

    dof = nansum(dof);
end
