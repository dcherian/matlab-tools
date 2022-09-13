% Returns maximum percentage error (in decimal) between A and B relative to
% the variable specified in rel (1 -> A, 2 -> B)
%       function [error] = calc_error(A,B,rel)

function [error] = calc_error(A,B,rel)
    
    if rel == 1
        err = abs(A-B)./A;
    else
        err = abs(A-B)./B;
    end
    
    error = max(addnan(err,Inf));
        