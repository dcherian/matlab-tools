% [floats] = roms_read_floats(type,file,rgrid)
% type = roms / ltrans / tracmass
function [floats] = read_floats(type,file,rgrid)

    if strcmpi(type,'roms')
        floats.x = ncread(file,'x')';
        floats.y = ncread(file,'y')';
        floats.z = ncread(file,'depth')';
        floats.time = ncread(file,'ocean_time');
        floats.temp = ncread(file,'temp')';
        floats.salt = ncread(file,'salt')';

        for ii = 1:size(floats.x,2)
            clear ind
            ind = find(floats.x(:,ii) > 0);
            if isempty(ind)
                ind = 1;
            else
                ind = ind(1);
            end

            floats.init(ii,:) = [floats.x(ind,ii) floats.y(ind,ii) floats.z(ind,ii) floats.time(ind)];
        end
        floats.type = 'roms';
    end
    
    if strcmpi(type,'ltrans')
        floats.y = ncread(file,'lat')';
        floats.x = ncread(file,'lon')';
        floats.z = ncread(file,'depth')';
        floats.age = ncread(file,'age')'; 
        floats.time = ncread(file,'model_time');
        da = -1*diff(floats.age == 0,1,1);
        da(end+1,:) = 0;
        tinit = cut_nan(fillnan(da .* repmat(floats.time,[1 size(da,2)]),0));
        
        floats.type = 'ltrans';
        floats.init = [floats.x(1,:)' floats.y(1,:)' floats.z(1,:)' tinit];
        dfx = (diff(floats.x,1,1) == 0);
        floats.x(dfx) = NaN;
        floats.y(dfx) = NaN;
        floats.z(dfx) = NaN;
        floats.age(dfx) = NaN;
    end
    
    if strcmpi(type,'tracmass')
       floats = tracmass_read(file,rgrid);
    end
    floats.comment = ['init = (x,y,z,t) = initial location, release time in meters, seconds | ' ...
        'fac = number float timesteps per ROMS output timestep'];
    % calculate fac
    dtroms = rgrid.ocean_time(2)-rgrid.ocean_time(1);
    dtltr  = floats.time(2,1)-floats.time(1,1);
    floats.fac = dtroms/dtltr;
    
% ASSUMES REGULAR GRID
function [floats] = tracmass_read(fname,rgrid)

% from the manual
% The trajectories are stored in <outDataDir>/<outDataFile>_run.asc, which has to be
% specified in <project>_run.in
% The trajectories are stored in colums of
%
%              ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens
%
% where
% ntrac is the trajectory number
% niter is the TRACMASS code iteration (only important for TRACMASS modellers)
% x1 is the zoonal position of the trajectory particle - INDEX
% y1 is the meridional position of the trajectory particle - INDEX
% z1 is the vertical position of the trajectory particle - INDEX
% tt the time of the trajectory particle (in days)
% t0 the initial time of the trajectory particle
% subvol is the the "volume" transport in m3/s of the trajectory
% temp is the temperature of the trajectory particle
% salt is the salinity/specific humidity of the trajectory particle
% dens is the density of the trajectory particle
   tic 
   [ntrac,~,ix,iy,iz,tt,t0,subvol,temp,salt,~] = ...
                textread(fname,'%d%d%f%f%f%f%f%f%f%f%s');
   
%    [data(:,1),~,data(:,2),data(:,3),data(:,4),data(:,5),data(:,6), ...
%             data(:,7),data(:,8),data(:,9),~] = ...
%                 textread(fname,'%d%d%f%f%f%f%f%f%f%f%s');
   toc
   disp('Finished opening file. Now processing for unique records');
   tic         
%   data = sortrows(data,1); % sort according to trac number
   floats.time = unique(tt);
   
   xr = rgrid.x_rho(1,:)';
   yr = rgrid.y_rho(:,1);
   
   dx = xr(2)-xr(1); dy = yr(2)-yr(1);
   cpb = progressbar();
   for i = 1:length(unique(ntrac)) % ith drifter
      %k = 1
       j = find(ntrac == i);
      %for j=1:length(ntrac)    
        %if ntrac(j) == i
            floats.t0(i) = t0(j(1));
            fx = floor(ix(j)); fy = floor(iy(j)); fz = floor(iz(j));
            cz = ceil(iz(j));
            fy(fy == 0) = 1; fx(fx == 0) = 1;
            % not all floats start at t=0
            dt = (floats.t0(i)-floats.time(1))./ (floats.time(2)-floats.time(1));
            k = 1:length(j);
            floats.x(k+dt,i) = xr(fx) + (ix(j)-fx) * dx;
            floats.y(k+dt,i) = yr(fy) + (iy(j)-fy) * dy;
            for kk = 1:length(fz)
                try
                    if fz(kk) == 0, z1 = NaN; end
                    z1 = rgrid.z_r(fz(kk),fy(kk),fx(kk));
                catch ME
                    disp('!!!')
                end
                z2 = rgrid.z_r(cz(kk),fy(kk),fx(kk));
            	dz = z2 - z1;
                floats.z(kk+dt,i) = z1 + (iz(j(kk))-fz(kk)) * dz;
            end
            floats.iz(k+dt,i) = iz(j);
            floats.temp(k+dt,i) = temp(j);
            floats.salt(k+dt,i) = salt(j);
            floats.t(k+dt,i)    = tt(j);
            floats.subvol(k+dt,i) = subvol(j);
            progressbarupdate(cpb,i/length(unique(ntrac)) * 100);
        %    k=k+1;
        %end
      %end
   end
   cpb.stop();
   
   for ii = 1:size(floats.x,2)
        clear ind
        ind = find(floats.x(:,ii) > 0);
        ind = ind(1);
        
        floats.init(ii,:) = [floats.x(ind,ii) floats.y(ind,ii) floats.z(ind,ii) floats.time(ind)];
   end
   
   % fill with nans
   names = fieldnames(floats);
   for ii=1:length(names)
      floats.(names{ii}) = fillnan(floats.(names{ii}),0); 
   end
   
   floats.fac = 1; % outputs at model output
   floats.comment = 'init = (x,y,z,t) = initial location, release time in meters, seconds';
   floats.type = 'tracmass';
   toc