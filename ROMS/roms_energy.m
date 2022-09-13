% Returns energy per unit mass (accounts for different cell volumes and divides by R0). Calculates growth rate in s^(-1)
% Use mean_index to say which dirn. you want to take the mean for defining mean, eddy contributions. By default, mat files 
%are created with integrated energy diagnostics.
%                   [EKE,MKE,PE] = roms_energy(fname,tindices,volume,ntavg,mean_index,commands)
%
% Accepts the following commands
%               - ncdf_out : controls whether netcdf file with energy files is written. 
%               - growthrate_u : calculated growth rate for 0/5u^2 instead of 0.5(u^2 + v^2)
% saves time of *first* peak in growth rate


function [EKE,MKE,PE,APE] = roms_energy(fname,tindices,volume,ntavg,mean_index,commands)

if ~exist('fname','var'), fname = 'ocean_his.nc'; end
if ~exist('tindices','var'), tindices = [1 Inf]; end
if ~exist('mean_index','var'), mean_index = 2; end
if ~exist('ntavg','var'), ntavg = 4; end % average over ntavg timesteps
if ~exist('volume','var'), volume = {}; end
if ~exist('commands','var'), commands = ''; end

command_list = {'ncdf_out','growthrate_u'};
[flag,~] = parse_commands(command_list,commands);
write_out = flag(1);
growthrate_u = flag(2);

ax = 'xyzt';

% parameters
vinfo = ncinfo(fname,'u');
dim   = length(vinfo.Size); 
s = vinfo.Size;
slab  = roms_slab(fname,1,ntavg);

% parse input
[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));
[xr,yr,zr,volr] = dc_roms_extract(fname,'rho',volume);
  
% mess around with volr to replicate boundaries as in https://www.myroms.org/wiki/index.php/Grid_Generation
volu = volr; volv = volr; volw = volr;
if ~isinf(volu(1,2)), volu(1,2) = volr(1,2)-1; end
if ~isinf(volv(2,2)), volv(2,2) = volr(2,2)-1; end
volw = volr;

 [xu,~,~,vu  ] = dc_roms_extract(fname,'u'  ,volume);
[~,yv,~,vv   ] = dc_roms_extract(fname,'v'  ,volume);
[~,~,zw,vw   ] = dc_roms_extract(fname,'w'  ,volume);

if vu(1,2) ~= volu(1,2), xu = xu(1:end-1); end
if vv(2,2) ~= volv(2,2), yv = yv(1:end-1); end

% trim to interior rho points
xr = xr(2:end-1);
yr = yr(2:end-1);

% calculate cell volumes use u/v/w points to mark edges of the cell
totalvol = (xr(end)-xr(1))*(yr(end)-yr(1))*max(abs(zr));
cellvol = abs(bsxfun(@times,diff(xu)*diff(yv),permute(diff(zw),[3 2 1])));
if totalvol./sum(cellvol(:)) > 1.5, error('cell volume calculation is wrong!'); end

%% read data

% caps indicates domain integrated values
EKE = nan(ceil(nt/ntavg)-1,1);
MKE = EKE;
PE  = EKE;
if growthrate_u, EKEu = EKE; end

R0  = ncread(fname,'R0');
h   = ncread(fname,'h');
h   = h(2:end-1,2:end-1);
time = ncread(fname,'ocean_time');
time = time([tindices(1):tindices(2)]);

%% create output file

if write_out
    outname = ['ocean_energy-' ax(mean_index) '.nc'];
    if exist(outname,'file')
        %in = input('File exists. Do you want to overwrite (1/0)? ');
        in =1;
        if in == 1, delete(outname); end
    end
    xname = 'x_en'; yname = 'y_en'; zname = 'z_en'; tname = 't_en';
    try
        nccreate(outname,'eke','Dimensions', {xname s(1)-1 yname s(2)-2 zname s(3) tname Inf});
        nccreate(outname,xname,'Dimensions',{xname s(1)-1});
        nccreate(outname,yname,'Dimensions',{yname s(2)-2});
        nccreate(outname,zname,'Dimensions',{zname s(3)});

        nccreate(outname,'mke','Dimensions', {xname s(1)-1 yname s(2)-2 zname s(3) tname Inf});   
        nccreate(outname,'pe','Dimensions', {xname s(1)-1 yname s(2)-2 zname s(3) tname Inf});

        ncwriteatt(outname,'eke','Description','EKE/horizontal area field');
        ncwriteatt(outname,'eke','coordinates','x_en y_en z_en t_en');
        ncwriteatt(outname,'eke','units','J/m^2');

        ncwriteatt(outname,'mke','Description','MKE/horizontal area field');
        ncwriteatt(outname,'mke','coordinates','x_en y_en z_en t_en');
        ncwriteatt(outname,'mke','units','J/m^2');

        ncwriteatt(outname,'pe','Description','PE/horizontal area field');
        ncwriteatt(outname,'pe','coordinates','x_en y_en z_en t_en');
        ncwriteatt(outname,'pe','units','J/m^2');

        ncwriteatt(outname,xname,'units',ncreadatt(fname,'x_u','units'));
        ncwriteatt(outname,yname,'units',ncreadatt(fname,'y_u','units'));
        ncwriteatt(outname,zname,'units','m');
        fprintf('\n Created file : %s\n', outname);
    catch ME
        fprintf('\n Appending to existing file.\n');
    end
end

%% Calculate!
tend = 0;
corr_flag = 0;
ax = 'xyz';
fprintf('\n Slab = %d | Averaging in %s over %d timestep(s)\n\n', slab, ax(mean_index),ntavg);

try
    cpb = progressbar();
catch ME
    cpb = [];
end

% initialize correction terms
for i=0:iend-1
    % FROM mod_movie.m - propagate changes back
    % start and count arrays for ncread : corrected to account for stride
    
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt,volu);
    
    % Read extra timesteps to account for averaging.
    if i > 0
        read_start(end) = read_start(end) - ntavg + 1;
        read_count(end) = read_count(end) + ntavg - 1;
    else
        APE = roms_ape(fname);
    end
    
    if isempty(cpb), fprintf('\nReading Data...\n'); end
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt,volu);
     u   = ncread(fname,'u',read_start,read_count,stride); pbar(cpb,i+1,1,iend,5);
     
     [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt,volv);
     v   = ncread(fname,'v',read_start,read_count,stride); pbar(cpb,i+1,2,iend,5);
     
    % [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt,volw);
    % w   = ncread(fname,'w',read_start,read_count,stride); pbar(cpb,i+1,3,iend,5);
     
     [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt,volr);
     rho = ncread(fname,'rho',read_start,read_count,stride); pbar(cpb,i+1,4,iend,5);
    %zeta = ncread(fname,'zeta',[read_start(1:2) read_start(end)],[read_count(1:2) read_count(end)],[stride(1:2) stride(end)]); pbar(cpb,i+1,5,iend,5);
	if isempty(cpb), fprintf('\n Done reading data... \n'); end
    
    % mean fields - average over ntavg timesteps
    um = time_mean2(u,ntavg,mean_index);
    vm = time_mean2(v,ntavg,mean_index);
   % wm = time_mean2(w,ntavg,mean_index);
    rm = time_mean2(rho,ntavg,mean_index);
    
    s = size(u);
    
    % eddy fields
    if mod(ntavg,2) == 0
        ind1 = 2:1:s(4)-2;
        ind2 = 3:1:s(4)-1;
    else if ntavg ~=1
        ind1 = 2:1:s(4)-1;
        ind2 = ind1;
        else
            ind1 = 1:size(um,4);
            ind2 = ind1;
        end
    end
    
    % pull out rho & zeta at timesteps where i'm calculating eddy fields.
    rho  = (rho(2:end-1,2:end-1,:,ind1) + rho(2:end-1,2:end-1,:,ind2))/2;
    %zeta = (zeta(2:end-1,2:end-1,ind1) + zeta(2:end-1,2:end-1,ind2))/2;
    
    up = bsxfun(@minus,(u(:,:,:,ind1) + u(:,:,:,ind2))/2,um);
    vp = bsxfun(@minus,(v(:,:,:,ind1) + v(:,:,:,ind2))/2,vm);
   % wp = bsxfun(@minus,(w(:,:,:,ind1) + w(:,:,:,ind2))/2,wm);
    
    % corrections for initial perturbation field ~= 0
    if i == 0 & tindices(1) == 1
        if abs(max(max(max(up(:,:,:,1))))) >= 0.01, corr_u = up(:,:,:,1); corr_flag = 1; end
        if abs(max(max(max(vp(:,:,:,1))))) >= 0.01, corr_v = vp(:,:,:,1); corr_flag = 1; end        
      %  if abs(max(max(max(wp(:,:,:,1))))) >= 0.01, corr_w = wp(:,:,:,1); corr_flag = 1; end
    end      
    
    if exist('corr_u','var'), up = bsxfun(@minus,up,corr_u); um = bsxfun(@plus,um,corr_u); end
    if exist('corr_v','var'), vp = bsxfun(@minus,vp,corr_v); vm = bsxfun(@plus,vm,corr_v); end
   % if exist('corr_w','var'), wp = bsxfun(@minus,wp,corr_w); wm = bsxfun(@plus,wm,corr_w);end
    
    % correct size of mean fields if not corrected for perturbation
    
    um = correct_size(um,mean_index,size(up));
    vm = correct_size(vm,mean_index,size(vp));
   % wm = correct_size(wm,mean_index,size(wp));
    
    clear u v% w
    
    % average so that everything lands up on interior-rho points
    up = (up(1:end-1,2:end-1,:,:) + up(2:end,2:end-1,:,:))/2;    
    vp = (vp(2:end-1,1:end-1,:,:) + vp(2:end-1,2:end,:,:))/2;
   % wp = (wp(2:end-1,2:end-1,1:end-1,:) + wp(2:end-1,2:end-1,2:end,:))/2;
    
    um = (um(1:end-1,2:end-1,:,:) + um(2:end,2:end-1,:,:))/2;    
    vm = (vm(2:end-1,1:end-1,:,:) + vm(2:end-1,2:end,:,:))/2;
   % wm = (wm(2:end-1,2:end-1,1:end-1,:) + wm(2:end-1,2:end-1,2:end,:))/2;
    
    tstart = tend+1;
    tend = tstart + s(4)-ntavg;
    
    % now calculate energy terms
    % 1/2 R0 u^2 = energy / unit vol in each cell. 
    % Integrate over domain and then divide by total mass in domain = energy / unit mass
    eke = 0.5*R0*(up.^2 + vp.^2); %bsxfun(@rdivide,,cellvol); % Boussinesq + wp.^2
    mke = 0.5*R0*(um.^2 + vm.^2); %bsxfun(@rdivide,,cellvol); % again Boussinesq + wm.^2
    pe  = 9.81*bsxfun(@times,rho,permute(zr,[2 3 1 4])); %bsxfun(@rdivide,,cellvol);  
    
    %oke = rho.*(up.*um + vp.*vm)./area;% -> should average (integrate) to zero theoretically | smaller order term
    
    t_en(tstart:tend,1) = (time(read_start(end)+ind1-1) + time(read_start(end)+ind2-1))/2;
    EKE(tstart:tend) = domain_integrate(eke,xr,yr,zr)./R0./totalvol;
    MKE(tstart:tend) = domain_integrate(mke,xr,yr,zr)./R0./totalvol;
    PE(tstart:tend)  = domain_integrate( pe,xr,yr,zr)./R0./totalvol;
    %PE(tstart:tend)  = domain_integrate2(-R0*9.81*(zeta.^2)/2,grid.x_rho(1,2:end-1),grid.y_rho(2:end-1,1))./area ...  % OTHER PE CODE - OLD      
    %                      + domain_integrate2(R0*9.81*h.^2/2,grid.x_rho(1,2:end-1),grid.y_rho(2:end-1,1))./area ...
    %                           + domain_integrate(pe,grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1)); 
    
    % if calculating growth rate only for u perturbations
    if growthrate_u
        ekeu = bsxfun(@rdivide,0.5*(up.^2),cellvol);
        EKEu(tstart:tend) = domain_integrate(ekeu,xr,yr,zr);
    end
    
    read_start(end) = tstart;
    if write_out
        % Write to netcdf file here
        ncwrite(outname,'eke',eke,read_start);   
        ncwrite(outname,'mke',mke,read_start);  
        ncwrite(outname,'pe' , pe,read_start);  
    end
end

if ~isempty(cpb)
    cpb.stop();
end

if write_out
    try 
        nccreate(outname,tname,'Dimensions',{tname length(t_en)});
        ncwriteatt(outname,tname,'units','s');
    catch ME
    end

    ncwrite(outname,xname,xr);
    ncwrite(outname,yname,yr);
    ncwrite(outname,zname,zr);
    ncwrite(outname,tname,t_en);
end

if exist('corr_u','var'), fprintf('\n Correction(s) applied to u perturbation fields \n\n'); end
if exist('corr_v','var'), fprintf('\n Correction(s) applied to v perturbation fields \n\n'); end
    
%% Calculate growth rate - definitely works!
k=1;
jump = 5; % fit 5 consecutive points

if growthrate_u, fitEKE = EKEu; else fitEKE = EKE; end

for i=1:1:length(fitEKE)-jump
    A(k,:) = polyfit(t_en(i:i+jump),log(fitEKE(i:i+jump)),1);
    time_A(k,:) = (t_en(i+jump)+t_en(i))/2;
    k=k+1;
end

time_A = time_A(2:end);
A = A(2:end,:);

figure
title('Growth Rate Curve')
subplot(211)
plot(time_A/86400,A(:,1)*86400,'b*-')
liney(0);
ylabel('Growth Rate (d^{-1})')
xlabel('Time (days)');
beautify; box on

% Verify
eke2 = exp(A(:,1).*time_A + A(:,2));
subplot(212)
plot(t_en/86400,(fitEKE),'b*');
hold on
plot(time_A/86400,eke2,'r');
ylabel('Energy');
xlabel('Time (days)');
title('Verification');
legend('Original','Fit');
beautify; box on

A = A(:,1);

%% plot

figure;
subplot(211)
title('Energy Curves')
hold on;
plot(t_en/86400,EKE,'r');
plot(t_en/86400,MKE,'k');
plot(t_en/86400,PE-min(PE),'b');
%plot(time,OKE,'m');
ylabel('Energy / unit mass (m^2/s^2)');
xlabel('Time (days)');
legend('EKE','MKE','PE - PE_{min}','Location','Best');
beautify;

subplot(212)
hold on;
plot(t_en/86400,EKE./APE,'r');
plot(t_en/86400,MKE./APE,'k');
plot(t_en/86400,(PE-min(PE))./APE,'b');
%plot(time,OKE,'m');
ylabel('Energy / APE');
xlabel('Time (days)');
title(['APE = ' num2str(APE) ' J/kg']);
beautify; 
% figure;
% plot(t_en/86400,PE);
% ylabel('Energy');
% xlabel('Time (days)');
% legend('PE');

% write to file
fname = ['energy-avg-' ax(mean_index) '.mat']; 
save(fname,'t_en','PE','EKE','MKE','A','time_A','ntavg','volume','commands','APE');

%% local functions

function [datam] = time_mean(data,n,mean_index)
    for ii = 1:n:size(data,4)-n+1
        ind = ceil(ii/n);
        datam(:,:,:,ind) = mean(mean(data(:,:,:,ii:ii+n-1),4),mean_index);
    end
        
function [datam] = time_mean2(data,n,mean_index)
    for ii = 1:size(data,4)-n+1
        datam(:,:,:,ii) = mean(mean(data(:,:,:,ii:ii+n-1),4),mean_index);
    end
    
function [datam] = time_mean3(data,n,mean_index)
    datam = mean(mean(data,4),mean_index);

function [out] = domain_integrate2(in,xax,yax)
    out = squeeze(trapz(xax,trapz(yax,in,2),1));
    
function [out] = correct_size(data,mean_index,s)
    if size(data,mean_index) == 1
        sd = ones([1 length(size(data))]);
        sd(mean_index) = s(mean_index);
        out = repmat(data,sd);
    else
        out = data;
    end
       
function [] = pbar(cpb,i,j,imax,jmax)
    if ~isempty(cpb)
        txt = sprintf(' Progress: i=%d, j=%d',i,j);
        progressbarupdate(cpb,(jmax*(i-1)+j)/(imax*jmax)*100,txt);
    end