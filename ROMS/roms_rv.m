% calculates Ertel PV at interior RHO points (horizontal plane) and midway between rho points in the vertical
%       [rv,xrv,yrv,zrv] = roms_rv(fname,tindices,outname)

function [rv,xrv,yrv,zrv] = roms_rv(fname,tindices,outname)

% parameters
%lam = 'rho';
vinfo = ncinfo(fname,'u');
s     = vinfo.Size;
dim   = length(s); 
slab  = roms_slab(fname,0)-3;

warning off
grid = roms_get_grid(fname,fname,1,1);
warning on

% parse input
if ~exist('tindices','var'), tindices = []; end

[iend,tindices,dt,~,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

rho0  = ncread(fname,'R0');
trv = ncread(fname,'ocean_time');
trv = trv([tindices(1):tindices(2)]);
f   = ncread(fname,'f',[1 1],[Inf Inf]);

xname = 'x_rv'; yname = 'y_rv'; zname = 'z_rv'; tname = 'ocean_time';

grid1.xv = repmat(grid.x_v',[1 1 grid.N]);
grid1.yv = repmat(grid.y_v',[1 1 grid.N]);
grid1.zv = permute(grid.z_v,[3 2 1]);

grid1.xu = repmat(grid.x_u',[1 1 grid.N]);
grid1.yu = repmat(grid.y_u',[1 1 grid.N]);
grid1.zu = permute(grid.z_u,[3 2 1]);

grid1.xr = repmat(grid.x_rho',[1 1 grid.N]);
grid1.yr = repmat(grid.y_rho',[1 1 grid.N]);
grid1.zr = permute(grid.z_r,[3 2 1]);

grid1.zw = grid.z_w;
grid1.s_w = grid.s_w;
grid1.s_rho = grid.s_rho;

%% setup netcdf file

if ~exist('outname','var') || isempty(outname), outname = 'ocean_rv.nc'; end
if exist(outname,'file')
    %in = input('File exists. Do you want to overwrite (1/0)? ');
    in = 1;
    if in == 1, delete(outname); end
end
try
    nccreate(outname,'rv','Dimensions', {xname s(1)-1 yname s(2)-2 zname s(3) tname length(trv)});
    nccreate(outname,xname,'Dimensions',{xname s(1)-1 yname s(2)-2 zname s(3)});
    nccreate(outname,yname,'Dimensions',{xname s(1)-1 yname s(2)-2 zname s(3)});
    nccreate(outname,zname,'Dimensions',{xname s(1)-1 yname s(2)-2 zname s(3)});
    nccreate(outname,tname,'Dimensions',{tname length(trv)});
    
    ncwriteatt(outname,'rv','Description','Relative vorticity (in z) calculated from ROMS output');
    ncwriteatt(outname,'rv','coordinates',[xname ' ' yname ' ' zname ' ' ocean_time]);
    ncwriteatt(outname,'rv','units','N/A');
    ncwriteatt(outname,xname,'units',ncreadatt(fname,'x_u','units'));
    ncwriteatt(outname,yname,'units',ncreadatt(fname,'y_u','units'));
    ncwriteatt(outname,zname,'units','m');
    ncwriteatt(outname,tname,'units','s');
    fprintf('\n Created file : %s\n', outname);
catch ME
    fprintf('\n Appending to existing file.\n');
end

%% calculate pv

misc = roms_load_misc(fname);

for i=0:iend-1
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    tstart = read_start(end);
    tend   = read_start(end) + read_count(end) -1;
    
    u      = ncread(fname,'u',read_start,read_count,stride);
    v      = ncread(fname,'v',read_start,read_count,stride);

    [rv,xrv,yrv,zrv] = rvor_cgrid(grid1,u,v);

    ncwrite(outname,'rv',rv,read_start); 
    
    % write now so that file is still usable in case of crash
    if i == 0
        ncwrite(outname,xname,xrv);
        ncwrite(outname,yname,yrv);
        ncwrite(outname,zname,zrv);
        ncwrite(outname,'ocean_time',trv);
    end
    
    intRV(tstart:tend) = domain_integrate(rv,xrv,yrv,zrv); 
end

save rv.mat rv xrv yrv zrv trv intRV
fprintf('\n Wrote file : %s \n\n',outname);