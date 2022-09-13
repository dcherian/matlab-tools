function [out] = roundoff(in,decimals)
    
    out = str2double(sprintf(['%.' num2str(decimals) 'f'],in));