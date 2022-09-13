% Plots labelled periodogram
%
% function [freq,pow,cn,ind] = dcpgram(data,dt,plot,numpeaks, factor)
%       data - time series
%       dt - sampling delta_time (default = 1)
%       numpeaks - top 'x' peaks to be marked
%       factor - value to scale labelled peak TIME PERIOD by
% Returns freq,pow,cn (from dcfft) and indexes of peaks in original data set
% Calls beautify.m at the end
% Plots by default.

function [freq,pow,cn,ind] = dcpgram(data,dt,plot,mark_peaks,factor)

    if ~exist('dt','var'), dt = 1; end;
    if ~exist('plot','var'), plot = 1; end

    data = squeeze(data);
    data = data-mean(data);
    
    [freq,pow,cn] = dcfft(data,dt); 

    pgram = abs(cn).^2;
    ind = NaN;

    if plot 
        figure
        loglog(freq,pgram,'LineWidth',1.5)
        title(['Periodogram']);
        grid on

        if exist('mark_peaks','var') & mark_peaks ~= 0
            [pks,ind] = findpeaks(pgram,'SORTSTR','descend');
            loglinex(freq(ind(1:mark_peaks)),factor);
            ind = ind(1:mark_peaks);
        else
            ind = NaN;
        end

        beautify;
    end