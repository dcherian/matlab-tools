% [out] = repnan(in,val) replaces NaN's in 'in' with specified 'val'

function [out] = repnan(in,val)
    aa = isnan(in);
    out = in;
    out(aa) = val;