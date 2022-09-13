% [x] = rednoise(sz)
% runningmean over rand(sz)
function [x] = rednoise(sz)

    if ~exist('sz', 'var')
        error('Provide shape of output');
    end

    SZ = ceil(sz/50);
    L = max(SZ);

    x = rand(sz + L*(sz ~= 1));
    x = conv(x, ones(SZ)./L, 'valid');
    x = x(1:end-1);

    assert(all(size(x) == sz))
end