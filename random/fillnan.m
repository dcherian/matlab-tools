% Replaces val with NaN's in 'in'
%       [in] = fillnan(in,val)

function [in] = fillnan(in,val,recursive)
    if islogical(in)
        % logicals cant have NaN
        return;
    end

    in(in == val) = NaN;

    if nargin == 2
        recursive = 0;
    end
    
    if val == Inf || val == -Inf
        if recursive == 1, return; end
        fprintf('\n Running fillnan recursively for +-Inf \n');
        in = fillnan(in,-1*val,1);
    end
        