% addnan(var,val,opt) replaces all values in (var) > val with NaN
% opt = abs, does abs(var) > abs(val). anything else does (var > val).
function [A] = addnan(var,val,opt)
    if ~exist('opt','var'), opt = []; end
    if strcmpi(opt,'abs')
        aa = abs(var) > abs(val);
    else
        aa = (var) > val;
    end
    A = var;
    A(aa) = NaN;