% Generate white noise time series [-1 1], 0 mean.
%        [noise] = whitenoise(sz)
function [noise] = whitenoise(sz)

    if ~exist('sz', 'var')
        error('Provide shape of output');
    end

    noise = rand(sz);
    noise = noise - mean(noise(:));
    noise = noise./max(abs(noise(:)));
end