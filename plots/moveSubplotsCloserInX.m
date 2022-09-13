function [] = moveSubplotsCloserInX(m,n,hax,factor)
    if ~exist('factor', 'var'), factor = 1/3; end

    Nax = length(hax);

    if n > 2
        error('moveSubplotsCloserInX: only n=2 supported.')
    end

    for ii=1:n:m*n
        dx = -hax(ii).OuterPosition + hax(ii+1).OuterPosition;
        dx = dx(1);

        hax(ii+1).OuterPosition(1) = hax(ii+1).OuterPosition(1) - dx * factor;
    end

end