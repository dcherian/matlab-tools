% Calculate spectrum for gappy time series.
% Provide subset length to average over
%      [S, freq] = GappySpectrum(in, SubsetLength)
function [S, freq] = GappySpectrum(in, dt, SubsetLength, nAvgBands)
    if size(in, 1) == 1, in = in'; end

    if ~exist('dt', 'var'), dt = 1; end

    if ~exist('SubsetLength', 'var')
        warning('SubsetLength not provided. Using longest valid segment.');
        SubsetLength = [];
    end

    if ~exist('nAvgBands', 'var')
        nAvgBands = 5;
    end

    data = ~isnan(in);
    cc = bwconncomp(data);

    kk = 1;
    maxdataidx = 1;
    for ii=1:cc.NumObjects
        if length(cc.PixelIdxList{ii}) > ...
                length(cc.PixelIdxList{maxdataidx})
            maxdataidx = ii;
        end

        % make sure there are no NaNs in segment
        assert(all(~isnan(in(cc.PixelIdxList{ii})) == 1));

        SegmentLength = length(cc.PixelIdxList{ii});

        if SegmentLength < SubsetLength
            continue;
        end

        for zz=1:SubsetLength:SegmentLength
            if zz+SubsetLength-1 > SegmentLength
                continue;
            end
            Subset = in(cc.PixelIdxList{ii}(zz:zz+SubsetLength-1));
            [S(:,kk), freq] = spectrum_band_avg(Subset, 1, 5, 'hann', 0);
            kk = kk + 1;
        end
    end

    if isempty(SubsetLength)
        [S, freq] = spectrum_band_avg(in(cc.PixelIdxList{maxdataidx}), ...
                                      dt, nAvgBands, 'hann', 0);
        disp(['Using longest segment: ' ...
              num2str(length(cc.PixelIdxList{maxdataidx})) ...
              ' elements']);
    end

    if ~exist('freq', 'var')
        error('None of the segments are as long as SubsetLength.');
    end

    S = nanmean(abs(S), 2);
end