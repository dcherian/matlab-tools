% snips out non-nan sections and concatenates them
%   [out] = cut_nan(in)

function [out] = cut_nan(in)
    
    if isvector(in)
        out = in(~isnan(in));
    else
        index = 1:size(in(:,1));
        indices = index(~isnan(in(:,1)));
        out = in(indices,:);
    end