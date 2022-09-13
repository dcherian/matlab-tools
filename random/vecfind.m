% returns indices where abs(a)(vector) >= abs(b) (vector). Biased high
%       [out] = vecfind(a,b)
% Freaks out if there are no results

function [out] = vecfind(a,b)
    if length(a) == 1
        out = find_approx(b,a,1);
        return;
    end

    if length(b) == 1
        out = find_approx(a,b,1);
        return;
    end

    try
        out = arrayfun(@(x) find(abs(a) >= abs(x),1,'first'), b );
    catch ME
        disp(ME);
        disp('vecfind: No result found?');
        out = NaN;
    end