% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
%       gcm_movie(fname, varname, tindices, axis, index, commands) 
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start end]
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice / location in proper units (finds nearest point)
%           commands - extra commands to be executed after each plot

function [] = gcm_movie4(fname, tindices, axis, index, commands) 

% Set parameters
% if length(varname) == 1
%     varname = upper(varname); 
% else
%     varname(1) = upper(varname(1));
% end

% if strcmp(varname,'Eta')
%     axis = 'z';
%     index = 1;
% end

fname = find_file(fname);

if ~exist('commands','var'), commands = ''; end

varnames = {'V','U','Temp','W'};

for i=1:4
    varname = varnames{i};
    
    vinfo  = ncinfo(fname,varname);
    dim    = length(vinfo.Size);
    stride = 400; % stride for ncread. read 'n' records in at a time

    % Set up grid
    [xax{i},yax{i},zax{i},xunits{i},yunits{i}] = gcm_var_grid(fname,varname);
    
    labels{i}.title = [varname ' (' ncreadatt(fname,varname,'units') ') | '];
    labels{i}.revz  = 0;
    labels{i}.time  = double(ncread(fname,'T'))./3600/24;
    labels{i}.tunits = 'days';
end

% all vars have same time dimension
if isempty(tindices), tindices = [1 vinfo.Size(end)]; end

if tindices(2) == 0 | tindices(2) == 1
    iend   = 1;
else
    iend   = ceil((tindices(2)-tindices(1))/stride);
end

%% Plot according to options

for i=0:iend-1
    read_start = ones(1,dim);
    read_end   = Inf(1,dim);
    read_start(end) = tindices(1)+stride*i;
    if i == (iend-1)
        read_end(end) = tindices(2)-stride*(i);
    else
        read_end(end) = stride*(i+1);
    end
    
    labels{i+1}.time = labels{i+1}.time(read_start(end):(read_start(end)+read_end(end)-1));
            
    switch axis
        case 'x'
            for ii=1:4
                varname = varnames{ii};
                % given location instead of index
                if ischar(index), index = find_approx(xax{ii},str2num(index),1); end
                % fix title string
                labels{ii}.title = [labels{ii}.title axis ' = ' sprintf('%5.2f', xax{ii}(index)) ' m | '];

                read_start(1) = index;
                read_end(1)   = 1;
                var = ncread(fname,varname,read_start,read_end);            
                dv{ii}  = double(squeeze(var));

                % take care of walls
                s   = size(dv{ii});
                if repnan(dv{ii}(1,:,:),0)   == zeros([1 s(2) s(3)]), dv{ii}(1,:,:) = NaN; end;
                if repnan(dv{ii}(end,:,:),0) == zeros([1 s(2) s(3)]), dv{ii}(end,:,:) = NaN; end;

                labels{ii}.xax = ['Y (' yunits ')'];
                labels{ii}.yax = 'Z (m)';
            end
            animate4(yax,zax,dv,labels,commands);

        case 'y'
            % given location instead of index
            if ischar(index), index = find_approx(yax,str2num(index),1); end
            % fix title string
            labels{ii}.title = [labels{ii}.title axis ' = ' sprintf('%5.2f', yax(index)) ' m | '];
            
            read_start(2) = index;
            read_end(2)   = 1;
            var = ncread(fname,varname,read_start,read_end);
            dv  = double(squeeze(var));
            
            % take care of walls
            s   = size(dv);
            if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:) = NaN; end;
            if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
            
            labels.xax = ['X (' xunits ')'];
            labels.yax = 'Z (m)';
            animate(xax,zax,dv,labels,commands);

        case 'z'
            % given location instead of index
            if ischar(index), index = find_approx(zax,str2num(index),1); end
            % fix title string
            labels.title = [labels.title axis ' = ' sprintf('%5.2f', zax(index)) ' m | '];
            
            if dim ~= 3 % catch zeta
                read_start(3) = index;
                read_end(3)   = 1;
            end
            var = ncread(fname,varname,read_start,read_end);
            dv  = double(squeeze(var));
            
            % take care of walls
            s   = size(dv);
            if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:)   = NaN; end;
            if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
            if repnan(dv(:,1,:),0)   == zeros([s(1) 1 s(3)]), dv(:,1,:)   = NaN; end;
            if repnan(dv(:,end,:),0) == zeros([s(1) 1 s(3)]), dv(:,end,:) = NaN; end;
            
            labels.xax = ['X (' xunits ')'];
            labels.yax = ['Y (' yunits ')'];
            animate(xax,yax,dv,labels,commands);

        otherwise
            fprintf('\n ERROR: Invalid axis label. \n\n');
    end
end
