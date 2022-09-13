function [out] = orderofmagn(in)
    out = floor(log10(max(abs(in(:)))));