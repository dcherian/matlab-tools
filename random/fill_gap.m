function [out] = fill_gap(in,interp,interp_len)

%   Fills gaps (of length < interp_len) along columns with scheme specified in 'interp'
%       [out] = fill_gap(in,interp,interp_len)
%           interp -> 'linear'
%           interp_len -> max. length over which to interpolate (optional
%           if not interpolating)

    if(~strcmpi(interp,'linear'))
        fprintf('\n ERROR: ONLY LINEAR INTERPOLATION ALLOWED.');
        return
    end

    s = size(in);

    if s(1) == 1
        in = in';
        s = size(in);
        transposeFlag = 1;
    else
        transposeFlag = 0;
    end

    for i=1:size(in, 2)
        [gapstart, gapend] = FindGaps(in(:,i));
        gaplens = gapend - gapstart + 1; %[start(1)-1;start(2:end)-stop(1:end-1);s(1)-stop(end)];
        g = find(gaplens <= interp_len);
        %num = zeros(length(ln),1);
        out(:,i) = in(:,i);

        nn=0;
        for j=1:length(g)
            start = gapstart(g(j));
            stop = gapend(g(j));

            if start == 1 | start >= s(1) | stop >= s(1)
                continue
            end

            out(start:stop,i) = (interp1([start-1, stop+1], ...
                                         [in(start-1,i),in(stop+1,i)], ...
                                         [start:stop]))';
            nn=nn+1;
        end
    end

    if transposeFlag
        out = out';
    end
end