% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
%       mod_movie(fname, varname, tindices, volume,axis, index, commands, isDir)
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start (stride) end] OR [start count] OR [start] - single snapshot
%           volume - cell array with location of selection edges
%                       {'x' x1 x2; 'y' y1 y2; ...}
%                    x1,x2 etc. can be numbers (index) or string with
%                    actual dimensional values
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice or
%                 - location in proper units (string; finds nearest point)
%                 - 'mid' will find midpoint of 'axis' for you
%           commands - extra commands to be executed after each plot
%                      (passed to animate.m - optional)

% CHANGELOG:
% Changed to use updated roms_var_grid output                   26 Feb 2012
% Now using roms_ncread_params & roms_tindices                  24 Feb 2012
% Bugfix in read_count - skipped one timestep everytime         23 Feb 2012
% Bugfixes - title index when timesteps(1) ~= 1,                22 Feb 2012
%          - change depth to -ve by default and warn user
% Made the animate call generic - woot! & bugfixes              21 Feb 2012
%   - Imposed same colorbar even when multiple strides are used
%   - Fixed bug with Esc quitting when multiple strides are used
%   - Now uses km for x,y axes if applicable
%   - Passes stride number (i) and dt to animate.m
% Extended tindices options and fixed major bug in              14 Feb 2012
%      displaying subsets
% Added 'mid' as possible input for index                       07 Feb 2012
% Fixed 'z' plane sections & added filename check               06 Feb 2012
% Fixed display of single timestep                              26 Jan 2012
% Reducd slab to 100 and fixed many bugs                        26 Jan 2012
%   - titles work perfectly now. Added labels.tmax
%   - throws error if given a 2d simulation and asking for x-y plot
% Changed to mod_movie.m + identifies model from file           26 Jan 2012
% Added slab feature                                            26 Jan 2012
% Added extra commands argument to be passed to animate.m
% Uses NaN's to remove zeros at walls to prevent this from mucking contour plots
% Added labels structure
% Original roms_movie version

function [h_plot] = mod_movie(fname, varname, tindices, volume, axis, index, commands, hax)

if ~exist('hax', 'var')
    figure;
    hax = gca;
else
    axes(hax);
end

% fname = find_file(fname);
% if isempty(fname)
%     error('No file found.');
% else
% 	fprintf('\n Using file: %s\n', fname)
% end

% check inputs
if ~exist('tindices','var'), tindices = []; end
if ~exist('commands','var'), commands = ''; end
if ~exist('isDir','var'), isDir = 0; end

% prevent error by removing skiplevels if file and not directory is
% specified
[flag_skiplevels,commands] = parse_commands({'skiplevels'},commands);

% if folder, loop through all .nc files
if isobject(fname)
    obj = fname;
    fname = obj.dir;
    grd = obj.rgrid;
    outfile = obj.out_file;
else
    grd = [];
end
if isdir(fname)
    if ~strcmpi(varname, 'pv') && ~strcmpi(varname, 'rv')
        files = dir([fname '/*his*.nc*']);
    else
        files = dir([fname '/*vor*.nc*']);
    end
    outfile = [fname '/' files(1).name];
end

insertAnnotation(['mod_movie(' fname ',' varname ')']);

labels.isDir = isDir;
%% model specific setup
gcm = 1;
try
    ncreadatt(fname,'/','MITgcm_version');
catch ME
    gcm = 0;
end

if gcm
    % fix lowercase - uppercase problem
    if length(varname) == 1
        varname = upper(varname);
    else
        varname(1) = upper(varname(1));
    end

    % Set up grid
    [xax,yax,zax,xunits,yunits] = gcm_var_grid(fname,varname);
    vol = [1 Inf; 1 Inf; 1 Inf;];
    time = double(ncread(fname,'T'))./3600/24;

else
    if length(varname) == 1
        varname = lower(varname);
    else
        varname(1) = lower(varname(1));
    end
    % set up grid for first time instant

    if isempty(grd)
        try
            grd = roms_get_grid(outfile,outfile,0,1);
        catch ME
            grd = fname;
        end
    end
    if ~strcmpi(varname, 'pv') && ~strcmpi(varname, 'rv')
        [xax,yax,zax,~] = dc_roms_extract(grd,varname,volume,1);
        [~,~,~,~,xunits,yunits] = dc_roms_var_grid(grd,varname);
    else
        [xax,yax,zax,~] = dc_roms_extract(outfile,varname,volume,1);
        [~,~,~,~,xunits,yunits] = dc_roms_var_grid(outfile,varname);
    end
    time = dc_roms_read_data(fname, 'ocean_time')/86400;
end

%% fix input and get needed info

vinfo  = ncinfo(outfile,varname);
dim    = length(vinfo.Size);
slab   = 100; % slab for ncread. read 'n' records in at a time -
              % faster response + save some memory?

[iend,tindices,dt,~,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

% set only possible axis and index for Eta / zeta / ubar / vbar
if dim == 3
    axis = 'z';
    stride = [1 1 dt];
    index  = 1;
    sliceax(index) = 0;
    plotx = xax;
    ploty = yax;
    axind = 3;
    labels.xax = ['X (' xunits ')'];
    labels.yax = ['Y (' yunits ')'];
    commands = [commands '; image'];
else
    switch axis
        case 'x'
            sliceax = xax(:,1);
            plotx = yax;
            ploty = zax;
            axind = 1;
            labels.xax = ['Y (' yunits ')'];
            labels.yax = 'Z (meter)';

        case 'y'
            sliceax = yax(1,:)';
            plotx = xax;
            ploty = zax;
            axind = 2;
            labels.xax = ['X (' xunits ')'];
            labels.yax = 'Z (meter)';

        case 'z'
            sliceax = squeeze(zax(ceil(vinfo.Size(1)/2),ceil(vinfo.Size(2)/2),:));
            plotx = xax;
            ploty = yax;
            axind = 3;
            labels.xax = ['X (' xunits ')'];
            labels.yax = ['Y (' yunits ')'];
            commands = [commands '; image'];

        case 's'
            sliceax = 1:size(zax,3);
            plotx = xax;
            ploty = yax;
            axind = 3;
            labels.xax = ['X (' xunits ')'];
            labels.yax = ['Y (' yunits ')'];
            commands = [commands '; image'];

%             if ischar(index) && str2double(index) > 0
%                 warning('Changed input depth %s m to -%s m', index, index);
%                 index = num2str(-1 * str2double(index));
%             end

        otherwise
            error('Invalid axis label.');
    end
end

labels.hax = hax;
% check if there is topography
if isempty (zax)
    flags.notopo = 1;
else
    if bsxfun(@minus,zax,zax(1,1,:)) == zeros(size(zax))
        warning('no topography detected.');
        flags.notopo = 1;
    else
        flags.notopo = 0;
    end
end

% given location instead of index
if axis == 'z' && flags.notopo == 0
    if ischar(index),index = abs(str2double(index)); end
else
    if strcmpi(index,'end') || strcmpi(index,'inf'),  index = vinfo.Size(axind); end
    if strcmpi(index,'mid'), index = num2str((sliceax(1)+sliceax(end))/2); end
    if ischar(index), index = find_approx(sliceax,str2double(index),1); end
end

% fix axes
if dim == 3, index = 1; end
switch axind
    case 1
        plotx = squeeze(plotx(index,:,:,:));
        ploty = squeeze(ploty(index,:,:,:));

    case 2
        plotx = squeeze(plotx(:,index,:,:));
        ploty = squeeze(ploty(:,index,:,:));
    case 3
        plotx = squeeze(plotx(:,:,1,:));
        ploty = squeeze(ploty(:,:,1,:));
end

%% Plot according to options

labels.tmax = time(tindices(1) + dt*floor((tindices(2)-tindices(1))/dt));
labels.revz  = 0;
labels.tunits = 'days';
labels.dt = dt;
labels.t0 = tindices(1)-1;
try % shouldn't work only for salt / other stuff i create
    vartitle = [varname ' (' ncreadatt(outfile,varname,'units') ') | '];
catch ME
    vartitle = varname;
end

if exist('obj', 'var')
    vartitle = [obj.name ' | ' vartitle];
end

if ~isempty(strfind(labels.yax,'degree')) ...
        || ~isempty(strfind(labels.xax,'degree'))
    labels.dar = 1;
else
    labels.dar = 0;
end

% fix title string
if axis == 'x' || axis == 'y'
    if sliceax(index) > 1000
        labels.title = [vartitle axis ' = ' sprintf('%5.2f', sliceax(index)/1000) ' km | '];
    else
        labels.title = [vartitle axis ' = ' sprintf('%5.2f', sliceax(index)) ' m | '];
    end
else if axis == 'z'
        labels.title = [vartitle 'z = ' num2str(index) ' m | '];
    else
        labels.title = [vartitle 's = ' num2str(index) ' | '];
    end
end
labels.mm_instance = [];

% add axis info to volume to pass to dc_roms_read_data
if axis == 's', axis = 'z'; end
volume{size(volume,1)+1,1} = axis;
volume{size(volume,1),2} = index;
volume{size(volume,1),3} = index;

labels.time = time(tindices(1):dt:tindices(2));
labels.stride = i;

% read data
if axis ~='z' || (axis == 'z' && (~ischar(index) || flags.notopo)) || dim == 3
    dv = dc_roms_read_data(fname, varname, tindices, volume, [], grd);
else
    % we need to interpolate
    error('upgrade mod_movie interpolation');
    grids.xax = xax;
    grids.yax = yax;
    grids.zax = zax;

    warning off
    read_start(3) = 1;
    read_count(3) = Inf;
    for mmm = 1:read_count(4)
        disp(['reading & interpolating timestep ' ...
              num2str(mmm) '/' num2str(read_count(4))]);
        data = squeeze(double(ncread(fname,varname, ...
                                     [read_start(1:3) read_start(4)+mmm-1], ...
                                     [read_count(1:3) 1],stride)));

        dv(:,:,mmm) = dc_roms_zslice_var(data,index,grids);
    end
    clear data
    warning on
    % pause since this takes a long time and I might go do something else
    commands = [commands; 'pause'];
    [dv(:,:,mmm),~,~] = roms_zslice_var(permute(data,[3 2 1]),NaN,index,grd);
    dv = permute(dv,[2 1 3]);
end

% take care of walls for mitgcm - NEEDS TO BE CHECKED
if gcm
    s = size(dv);
    s(3) = size(dv,3); % correct for single timestep snapshots - in that case s is a 2 element row
    if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:) = NaN; end;
    if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;

    if axind == 3
        if repnan(dv(:,1,:),0)   == zeros([s(1) 1 s(3)]), dv(:,1,:)   = NaN; end;
        if repnan(dv(:,end,:),0) == zeros([s(1) 1 s(3)]), dv(:,end,:) = NaN; end;
    end
end

s = size(dv);
if s(1) == 1 || s(2) == 1
    close(gcf);
    error('2D simulation?');
end

if max(plotx(:)) > 1000
    plotx = plotx/1000;
    labels.xax = [labels.xax ' x 10^3'];
end
if max(ploty(:)) > 1000
    ploty = ploty/1000;
    labels.yax = [labels.yax ' x 10^3'];
end

% send to animate
[labels.mm_instance,h_plot] = animate(plotx,ploty,dv,labels,commands,3);

% for movie
if ~isempty(labels.mm_instance), mm_render(labels.mm_instance); end

end