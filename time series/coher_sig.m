% significance level for /coherence squared/ not /coherence/
% nf -> effective number of frequency components in spectral window
%    = (degrees of freedom)/2
% alpha = 0.05 for 95% significance
% Paul R Julian (1975) - Comments on the determination of the significance 
%                        levels of the coherence statistic

function [sig] = coher_sig(alpha, nf)

    sig = 1-alpha^(1/(nf -1));
