% [zwidth, idiff] = findProfilePeakWidth(vec, zvec, debug)
% Finds peak width based on derivative criterion
% ASSUMES ONLY ONE PEAK EXISTS!

function [zwidth, idifflo, idiffhi] = findProfilePeakWidth(vec, zvec, debug, test)

    if ~exist('debug', 'var'), debug = 0; end
    if ~exist('zvec', 'var'), zvec = 1:length(vec); end
    if ~exist('test', 'var'), test = 0; end

    if test, test_function; return; end

    zmax = NaN; imax = NaN; zwidth = NaN; idifflo = []; idiffhi = []; idiff = [];

    if size(vec, 1) == 1, vec = vec'; end
    if size(zvec, 1) == 1, zvec = zvec'; end

    nsmooth = 3;

    dvec = smooth(diff(vec,1,1)./diff(zvec,1,1), nsmooth);
    dzvec = avg1(zvec);

    % peak
    [vmax,ivmax] = max(abs(vec));
    zmax = zvec(ivmax);

    % peakwidth based on derivative drop to half it's max
    % value — look for first depth after maximum
    if ivmax > length(dvec), ivmax = length(dvec); end
    [dvmax, idvmax] = max(abs(dvec(1:ivmax)));
    idifflo = find(abs(dvec(1:min(ivmax,idvmax))) - dvmax/2 > 0, 1, 'first');

    [dvmax, idvmax] = max(abs(dvec(ivmax:end)));
    idvmax = idvmax + ivmax - 1;
    istart = max(ivmax,idvmax);
    idiffhi = find(abs(dvec(istart:end)) - dvmax/2 < 0, 1, 'first') ...
              + istart - 1;

    if idvmax == 1 | idvmax == length(dvec) % min at one end
        idiffhi = idvmax;
    end

    if (isempty(idifflo) | isempty(idiffhi)) % & ~isempty(findstr(runs.name, 'ew-8'))
        % max slope is near bottom, so search up to location of actual maximum.
        %idiff = find(dvec(idmax+10:ivmax) - dvmax/2 < 0, 1, 'first');
        % just look for (peak value)/2
        idifflo = find(vec - vmax/2 < 0, 1, 'last');
        idiffhi = 1;
    else
        % back to actual vector grid (not diff vector)
        if idifflo ~= 1, idifflo = idifflo - 1; end
        idiffhi = idiffhi + 1;
        zwidth = abs(zvec(idifflo) - zvec(idiffhi));
    end

    if debug
        figure;
        hh(1) = plot(vec, zvec);
        liney(zmax, 'max');
        %liney(zbot, 'bot');
        try
            liney([zvec(idifflo) zvec(idiffhi)], {'lo'; 'hi'});
        catch ME; end
        linex(0);
        hh(2) = plot(diff(vec,1,1)./diff(zvec)*10, avg1(zvec));
        legend(hh, {'vector', 'derivative'});
        beautify;
    end

    if isempty(idifflo) & isempty(idiffhi)
        warning('findProfilePeakWidth: Derivative criterion not satisfied.');
        return;
    end


end

function [] = test_function()

    xivec = linspace(0,3,100);
    vec = xivec .* exp(-xivec.^2);

    findProfilePeakWidth(vec, xivec, 1);
end