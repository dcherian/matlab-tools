% Given a cell 'volume' with locations, figures out indices on the axes and
% returns correct axes for the section. 
%           [xax,yax,zax,vol] = roms_extract(fname/grd,varname,volume)
%   volume is cell array of form
%           {'x' 'location1' 'location2';
%            'y' 'location1' 'location2';
%            'z' 'location1' 'location2'};
% locations in strings - in units else provide axis index as number

function [xax,yax,zax,vol] = dc_roms_extract(fname,varname,volume,tindex)

    if ~exist('tindex','var'), tindex = 0; end

    [xax,yax,zax,~,~] = dc_roms_var_grid(fname,varname,tindex);

    % Assuming cartesian grid
    xaxis = xax(:,1,end,end);
    yaxis = yax(1,:,end,end)';
    
    vol = [1 Inf; 1 Inf; 1 Inf]; % default - choose all data
    
    if ~exist('volume','var'), return; end
    if isempty(volume), return; end
    
    % parse volume
    sa = size(volume);
    for i=1:sa(1)
        switch volume{i,1}
            case 'x'
                if ischar(volume{i,2}), volume{i,2} = find_approx(xaxis,str2double(volume{i,2}),1); end
                if ischar(volume{i,3}), volume{i,3} = find_approx(xaxis,str2double(volume{i,3}),1); end
                
                if isinf(volume{i,3}), volume{i,3} = length(xaxis); end
                if isinf(volume{i,2}), volume{i,2} = length(xaxis); end
                
                xax = xax(volume{i,2}:volume{i,3},:,:);
                yax = yax(volume{i,2}:volume{i,3},:,:);
                if ~isempty(zax) % for zeta
                    zax = zax(volume{i,2}:volume{i,3},:,:);
                end
                
                vol(1,1) = volume{i,2};
                vol(1,2) = volume{i,3};
                
            case 'y'
                if ischar(volume{i,2}), volume{i,2} = find_approx(yaxis,str2double(volume{i,2}),1); end
                if ischar(volume{i,3}), volume{i,3} = find_approx(yaxis,str2double(volume{i,3}),1); end
                
                if isinf(volume{i,3}), volume{i,3} = length(yaxis); end
                if isinf(volume{i,2}), volume{i,2} = length(yaxis); end
                
                xax = xax(:,volume{i,2}:volume{i,3},:);
                yax = yax(:,volume{i,2}:volume{i,3},:);
                if ~isempty(zax) % for zeta
                    zax = zax(:,volume{i,2}:volume{i,3},:);
                end
                
                vol(2,1) = volume{i,2};
                vol(2,2) = volume{i,3};
                
            case 'z'
                if ischar(volume{i,2}), volume{i,2} = find_approx(zax,str2double(volume{i,2}),1); end
                if ischar(volume{i,3}), volume{i,3} = find_approx(zax,str2double(volume{i,3}),1); end
                
                if isinf(volume{i,3}), volume{i,3} = size(zax,3); end
                if isinf(volume{i,2}), volume{i,2} = size(zax,3); end
                
                if volume{i,3} < volume{i,2}
                    temp = volume{i,3};
                    volume{i,3} = volume{i,2};
                    volume{i,2} = temp;
                end

                try
                    % if no z-dimension, then this fails.
                    xax = xax(:,:,volume{i,2}:volume{i,3});
                    yax = yax(:,:,volume{i,2}:volume{i,3});
                end
                zax = zax(:,:,volume{i,2}:volume{i,3});
                
                vol(3,1) = volume{i,2};
                vol(3,2) = volume{i,3};
                
            otherwise
                error('Unknown option in volume');
        end
    end