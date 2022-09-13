% returns 10^in;
% USAGE: txp in or txp(in)

function [out] = txp(in)
    if ischar(in)
        in = str2double(in);
    end
    
    out = 10^(in);