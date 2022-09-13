% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
% The plot is of difference w.r.t field at tindices(1) of the variable
%       roms_diffmovie(fname, varname, tindices,axis, index) 
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start end]
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice

function [] = roms_diffmovie(fname, varname, tindices, axis, index,commands) 

%[vars, atts, dims] = ncdfread(fname);vinfo = ncinfo(fname,varname);
if isempty(tindices)
        tindices = [1 Inf];
end

vinfo = ncinfo(fname,varname);
if length(vinfo.Size) == 3
    var  = ncread(fname,varname,[1 1 tindices(1)],[Inf Inf (tindices(2)-tindices(1))]);
else
    var  = ncread(fname,varname,[1 1 1 tindices(1)],[Inf Inf Inf (tindices(2)-tindices(1))]);
end

grid = roms_get_grid(fname,fname,0,1);
%var  = getfield(vars,varname);
s    = size(var);

if strcmp(varname, 'zeta');
    xax = grid.lon_v';
    yax = grid.lat_u';

else
    % from roms_islice.m
    % determine where on the C-grid these values lie 
    varcoords = nc_attget(fname,varname,'coordinates');
    if ~isempty(findstr(varcoords,'_u'))
      pos = 'u';
    elseif ~isempty(findstr(varcoords,'_v'))
      pos = 'v';
    elseif ~isempty(findstr(varcoords,'_rho'))
      pos = 'rho';
    else
      error('Unable to parse the coordinates variables to know where the data fall on C-grid')
    end

    switch pos
        case 'u'
            xax = grid.lon_u';
            yax = grid.lat_u';
            zax = grid.z_u;

        case 'v'
            xax = grid.lon_v';
            yax = grid.lat_v';
            zax = grid.z_v;

        otherwise
            xax = grid.lon_rho';
            yax = grid.lat_rho';
            zax = grid.z_r;
    end
end

% Plot according to options
switch axis
    case 'x'
        dv = zeros([s(2) s(3) s(4)]);
        for i=1:s(end)
            dv(:,:,i) = squeeze(var(index,:,:,i)-var(index,:,:,1));
        end
        animate(yax(1,: ),zax(:,1,1),dv,commands);
    
    case 'y'
        dv = zeros([s(1) s(3) s(4)]);
        for i=1:s(end)
            dv(:,:,i) = squeeze(var(:,index,:,i)-var(:,index,:,1));
        end
        animate(xax(:,1),zax(:,1,1),dv,commands);

    case 'z'
        if strcmp(varname, 'zeta')
            dv = zeros([s(1) s(2) s(3)]);
            for i=1:s(end)
                dv(:,:,i) = squeeze(var(:,:,i)-var(:,:,1));
            end
            animate(xax(:,1),yax(1,:),dv,commands);
        else            
            dv = zeros([s(1) s(2) s(4)]);
            for i=1:s(end)
                dv(:,:,i) = squeeze(var(:,:,index,i)-var(:,:,index,1));
            end
            animate(xax(:,1),yax(1,:),dv,commands);
        end
    otherwise
        fprintf('\n ERROR: Invalid axis label. \n\n');
end
