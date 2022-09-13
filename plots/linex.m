% Plot vertical line at x
%       [handles] = linex(x,label,color)

function [handles, htxt] = linex(x,label,color)
    
    if ~exist('label','var'), label = []; end

    if ~exist('color','var'), color = []; end
    
    [handles, htxt] = dcline('x',x,label,color);
