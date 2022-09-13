% one colorbar to rule them all
%   [handle] = subplot_cbar(caxis_vec,position)

function [handle] = subplot_cbar(caxis_vec,position)
    
    if ~exist('position','var'), position = [0.05 0.05 0.93 0.9]; end
    % Create colorbar
    axes('Position', position, 'Visible', 'off');
    handle = colorbar;
    caxis(caxis_vec);