function [hplt] = PlotSpectrum(in, dt, SubsetLength, nbands, hax)

    if ~exist('hax', 'var'), hax = gca; end
    if ~exist('dt', 'var'), dt = 1; end
    if ~exist('nbands', 'var'), nbands = 5; end
    if ~exist('SubsetLength', 'var'), SubsetLength = []; end

    [S, F] = GappySpectrum(in, dt, SubsetLength, nbands);
    axes(hax); hold on; length(F);
    hplt = plot(F, S, '.-');
    hax.XScale = 'log';
    hax.YScale = 'log';

    xlabel('Frequency')
    ylabel('Spectral density (m^2/s^2/freq)')

    %[F,S] = mspec(in, [], 'cyclic');
    %plot(F, S);
