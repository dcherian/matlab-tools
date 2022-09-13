% Plots horizontal line at a given y (can be a vector)
%       [handles] = liney(y,label,color)

function [handles, htxt] = liney(y,label,color)
    
    if ~exist('label','var'), label = []; end
    if ~exist('color','var'), color = []; end
    [handles, htxt] = dcline('y',y,label,color);
