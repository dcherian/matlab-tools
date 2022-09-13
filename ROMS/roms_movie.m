% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
%       roms_movie(fname, varname, tindices, axis, index, commands) 
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start (stride) end]
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice / location in proper units (finds nearest point)
%           commands - extra commands to be exxecuted after each plot or one of below
%               > nocaxis - non-constant colorbar

function [] = roms_movie(fname, varname, tindices, axis, index, commands) 

    mod_movie(fname, varname, tindices, axis, index, commands);

% fname = find_file(fname);
% fprintf('\n Using file: %s\n', fname)
% 
% % Set parameters
% vinfo  = ncinfo(fname,varname);
% dim    = length(vinfo.Size);
% jump   = 400; % jump for ncread. read 'n' records in at a time
% 
% if strcmp(varname,'zeta')
%     axis = 'z';
%     index = 1;
% end
% 
% if ~exist('commands','var'), commands = ''; end
% 
% if ~exist('tindices','var') || isempty(tindices)
%         tindices = [1 Inf];
%         dt = 1;
% else    
%     if length(tindices) == 3
%         dt = tindices(2);
%         tindices(2) = tindices(3);
%         tindices(3) = NaN;
%     else
%         if length(tindices) == 2
%             dt = 1;
%         end    
%     end
% end
% 
% if isinf(tindices(2)), tindices(2) = vinfo.Size(end); end
% 
% stride = [1 1 1 dt];    
% 
% if tindices(2) == 0 | tindices(2) == 1
%     iend   = 1;
% else
%     iend   = ceil((tindices(2)-tindices(1))/jump);
% end
% 
% % Set up grid
% [xax,yax,zax,xunits,yunits] = roms_var_grid(fname,varname);
% 
% %% Plot according to options
% 
% labels.title = [varname ' (' ncreadatt(fname,varname,'units') ') | '];
% labels.revz  = 0;
% labels.time  = double(ncread(fname,'ocean_time'))./3600/24;
% labels.tunits = 'days';
% 
% for i=0:iend-1
%     read_start = ones(1,dim);
%     read_end   = Inf(1,dim);
%     read_start(end) = tindices(1)+jump*i;
%     if i == (iend-1)
%         read_end(end) = ceil(tindices(2)/dt-jump*(i));
%     else
%         read_end(end) = jump*(i+1);
%     end
%     
%     labels.time = labels.time(read_start(end):dt:(read_start(end)+read_end(end)-1)*dt);
%             
%     switch axis
%         case 'x'
%             % given location instead of index
%             if ischar(index), index = find_approx(xax(:,1),str2num(index),1); end
%             % fix title string
%             labels.title = [labels.title axis ' = ' sprintf('%5.2f', xax(index,1)) ' m | '];
%             
%             read_start(1) = index;
%             read_end(1)   = 1;
%             var = ncread(fname,varname,read_start,read_end,stride);
%             dv  = double(squeeze(var));
%             
%             labels.xax = ['Y (' yunits ')'];
%             labels.yax = 'Z (m)';
%             animate(yax(1,: ),zax(:,1,1),dv,labels,commands);
% 
%         case 'y'
%             % given location instead of index
%             if ischar(index), index = find_approx(yax(1,:),str2num(index),1); end
%             % fix title string
%             labels.title = [labels.title axis ' = ' sprintf('%5.2f', yax(1,index)) ' m | '];
%             
%             read_start(2) = index;
%             read_end(2)   = 1;
%             var = ncread(fname,varname,read_start,read_end,stride);
%             dv  = double(squeeze(var));
%             
%             labels.xax = ['X (' xunits ')'];
%             labels.yax = 'Z (m)';
%             animate(xax(:,1),zax(:,1,1),dv,labels,commands);
% 
%         case 'z'
%             % given location instead of index
%             if ischar(index), index = find_approx(zax(:,1,1),str2num(index),1); end
%             % fix title string
%             labels.title = [labels.title axis ' = ' sprintf('%5.2f', zax(index,1,1)) ' m | '];
%             
%             if dim ~= 3 % catch zeta
%                 read_start(3) = index;
%                 read_end(3)   = 1;
%                 stride = [1 1 dt];
%             end
%             var = ncread(fname,varname,read_start,read_end,stride);
%             dv  = double(squeeze(var));
%             
%             labels.xax = ['X (' xunits ')'];
%             labels.yax = ['Y (' yunits ')'];
%             animate(xax(:,1),yax(1,:),dv,labels,commands);
% 
%         otherwise
%             fprintf('\n ERROR: Invalid axis label. \n\n');
%     end
% end
