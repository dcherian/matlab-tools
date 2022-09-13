function [] = moveSubplotsCloserInY(m,n,hax,factor)
    if ~exist('factor', 'var'), factor=-1/8; end

    Nax = length(hax);

    for ii=1:n
        dy = hax(ii).OuterPosition - hax(ii+n).OuterPosition;
        dy = dy(2);

        hax(ii).OuterPosition(2) = hax(ii).OuterPosition(2) + dy * factor;
        hax(ii+n).OuterPosition(2) = hax(ii+n).OuterPosition(2) - dy * factor;
    end

end