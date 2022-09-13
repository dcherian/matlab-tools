% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
%       gcm_movie(fname, varname, tindices, axis, index, commands) 
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start (stride) end]
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice / location in proper units (finds nearest point)
%           commands - extra commands to be executed after each plot

function [] = gcm_movie(fname, varname, tindices, axis, index, commands)

    mod_movie(fname, varname, tindices, axis, index, commands);

% fname = find_file(fname);
% fprintf('\n Using file: %s\n', fname)
% 
% % Set parameters
% if length(varname) == 1
%     varname = upper(varname); 
% else
%     varname(1) = upper(varname(1));
% end
% 
% if strcmp(varname,'Eta')
%     axis = 'z';
%     index = 1;
% end
% 
% if ~exist('commands','var'), commands = ''; end
% 
% vinfo  = ncinfo(fname,varname);
% dim    = length(vinfo.Size);
% jump = 400; % jump for ncread. read 'n' records in at a time
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
% [xax,yax,zax,xunits,yunits] = gcm_var_grid(fname,varname);
% 
% %% Plot according to options
% 
% labels.title = [varname ' (' ncreadatt(fname,varname,'units') ') | '];
% labels.revz  = 0;
% labels.time  = double(ncread(fname,'T'))./3600/24;
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
%             if ischar(index), index = find_approx(xax,str2num(index),1); end
%             % fix title string
%             labels.title = [labels.title axis ' = ' sprintf('%5.2f', xax(index)) ' m | '];
%             
%             read_start(1) = index;
%             read_end(1)   = 1;
%             var = ncread(fname,varname,read_start,read_end,stride);            
%             dv  = double(squeeze(var));
%             
%             % take care of walls
%             s   = size(dv);
%             if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:) = NaN; end;
%             if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
%             
%             labels.xax = ['Y (' yunits ')'];
%             labels.yax = 'Z (m)';
%             animate(yax,zax,dv,labels,commands);
% 
%         case 'y'
%             % given location instead of index
%             if ischar(index), index = find_approx(yax,str2num(index),1); end
%             % fix title string
%             labels.title = [labels.title axis ' = ' sprintf('%5.2f', yax(index)) ' m | '];
%             
%             read_start(2) = index;
%             read_end(2)   = 1;
%             var = ncread(fname,varname,read_start,read_end,stride);
%             dv  = double(squeeze(var));
%             
%             % take care of walls
%             s   = size(dv);
%             if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:) = NaN; end;
%             if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
%             
%             labels.xax = ['X (' xunits ')'];
%             labels.yax = 'Z (m)';
%             animate(xax,zax,dv,labels,commands);
% 
%         case 'z'
%             % given location instead of index
%             if ischar(index), index = find_approx(zax,str2num(index),1); end
%                         
%             if dim == 3 % catch zeta
%                 read_start(3) = index;
%                 %read_end(3)   = 1;
%                 stride = [1 1 dt];
%                 % fix title string
%                 labels.title = [labels.title axis ' = 0 m | '];
%             else
%                 % fix title string
%                 labels.title = [labels.title axis ' = ' sprintf('%5.2f', zax(index)) ' m | '];
%             end
%             var = ncread(fname,varname,read_start,read_end,stride);
%             dv  = double(squeeze(var));
%             
%             % take care of walls
%             s   = size(dv);
%             if length(s) == 2
%                 error('2D simulation?');
%             end
%             if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:)   = NaN; end;
%             if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
%             if repnan(dv(:,1,:),0)   == zeros([s(1) 1 s(3)]), dv(:,1,:)   = NaN; end;
%             if repnan(dv(:,end,:),0) == zeros([s(1) 1 s(3)]), dv(:,end,:) = NaN; end;
%             
%             labels.xax = ['X (' xunits ')'];
%             labels.yax = ['Y (' yunits ')'];
%             animate(xax,yax,dv,labels,commands);
% 
%         otherwise
%             fprintf('\n ERROR: Invalid axis label. \n\n');
%     end
% end
