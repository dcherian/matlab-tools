% calculates Ertel PV
%       [pv] = roms_pv(fname,tindices)

function [pv,xpv,ypv,zpv] = roms_pv(fname,tindices,outname)

% parameters
lam = 'rho';
vinfo = ncinfo(fname,'u');
s     = vinfo.Size;
dim   = length(s); 
slab  = 40;

warning off
grid = roms_get_grid(fname,fname,0,1);
warning on

% parse input
if ~exist('tindices','var'), tindices = []; end

[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

rho0  = ncread(fname,'R0');
tpv = ncread(fname,'ocean_time');
tpv = tpv([tindices(1):tindices(2)]);
f   = ncread(fname,'f',[1 1],[Inf Inf]);
f   = mean(f(:));

xpv = grid.x_rho(1,2:end-1)';
ypv = grid.y_rho(2:end-1,1)';
zpv = avg1(grid.z_r(:,1,1));

xname = 'x_pv'; yname = 'y_pv'; zname = 'z_pv'; tname = 'ocean_time';

%% setup netcdf file
if ~exist('outname','var') || isempty(outname), outname = 'ocean_pv.nc'; end
if exist(outname,'file')
    in = input('File exists. Do you want to overwrite (1/0)? ');
    if in == 1, delete(outname); end
end
try
    nccreate(outname,'pv','Dimensions', {xname s(1)-1 yname s(2)-2 zname s(3)-1 tname length(tpv)});
    nccreate(outname,xname,'Dimensions',{xname s(1)-1});
    nccreate(outname,yname,'Dimensions',{yname s(2)-2});
    nccreate(outname,zname,'Dimensions',{zname s(3)-1});
    nccreate(outname,tname,'Dimensions',{tname length(tpv)});
    
    ncwriteatt(outname,'pv','Description','Ertel PV calculated from ROMS output');
    ncwriteatt(outname,'pv','coordinates','x_pv y_pv z_pv ocean_time');
    ncwriteatt(outname,'pv','units','N/A');
    ncwriteatt(outname,xname,'units',ncreadatt(fname,'x_u','units'));
    ncwriteatt(outname,yname,'units',ncreadatt(fname,'y_u','units'));
    ncwriteatt(outname,zname,'units','m');
    ncwriteatt(outname,tname,'units','s');
    fprintf('\n Created file : %s\n', outname);
catch ME
    fprintf('\n Appending to existing file.\n');
end

ncwrite(outname,xname,xpv);
ncwrite(outname,yname,ypv);
ncwrite(outname,zname,zpv);
ncwrite(outname,'ocean_time',tpv);

%% calculate pv
pv = nan([s(1)-1 s(2)-2 s(3)-1 tindices(2)-tindices(1)+1]);

for i=0:iend-1
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    tstart = read_start(end);
    tend   = read_start(end) + read_count(end) -1;
    
    u      = ncread(fname,'u',read_start,read_count,stride);
    v      = ncread(fname,'v',read_start,read_count,stride);
    lambda = ncread(fname,lam,read_start,read_count,stride); % theta

    % calculate gradients
    vx    = bsxfun(@rdivide,diff(v,1,1),diff(grid.x_v',1,1)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
    vy    = bsxfun(@rdivide,diff(v,1,2),diff(grid.y_v',1,2)); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    vz    = bsxfun(@rdivide,diff(v,1,3),permute(diff(grid.z_v,1,1),[3 2 1])); %diff(v,1,3)./repmat(permute(diff(grid.z_v,1,1),[3 2 1]),[1 1 1 s(4)]);

    ux    = bsxfun(@rdivide,diff(u,1,1),diff(grid.x_u',1,1)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
    uy    = bsxfun(@rdivide,diff(u,1,2),diff(grid.y_u',1,2)); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    uz    = bsxfun(@rdivide,diff(u,1,3),permute(diff(grid.z_u,1,1),[3 2 1])); %diff(v,1,3)./repmat(permute(diff(grid.z_v,1,1),[3 2 1]),[1 1 1 s(4)]);

    tx    = bsxfun(@rdivide,diff(lambda,1,1),diff(grid.x_rho',1,1)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
    ty    = bsxfun(@rdivide,diff(lambda,1,2),diff(grid.y_rho',1,2)); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    tz    = bsxfun(@rdivide,diff(lambda,1,3),permute(diff(grid.z_r,1,1),[3 2 1])); %diff(v,1,3)./repmat(permute(diff(grid.z_v,1,1),[3 2 1]),[1 1 1 s(4)]);
    
    % PV calculated at interior rho points
                                % f + vx - uy                      (rho)_z
    pv(:,:,:,tstart:tend) = -1* double((avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  tz(2:end-1,2:end-1,:,:) ...
                   - avgy(vz(2:end-1,:,:,:)).*avgz(avgx(tx(:,2:end-1,:,:))) ... % vz * (rho)_x
                   + avgx(uz(:,2:end-1,:,:)).*avgz(avgy(ty(2:end-1,:,:,:))))./rho0);%avgz(lambda(2:end-1,2:end-1,:,:))); % uz*(rho)_y

    ncwrite(outname,'pv',pv(:,:,:,tstart:tend),read_start); 
    
    debug = 1;
    
    if debug
        pv1 = -avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  tz(2:end-1,2:end-1,:,:);
        pv2 = avgy(vz(2:end-1,:,:,:)).*avgz(avgx(tx(:,2:end-1,:,:)));
        pv3 = avgx(uz(:,2:end-1,:,:)).*avgz(avgy(ty(2:end-1,:,:,:)));
        
        tind = 1;
        yind = 3;
        
        figure;
        contourf(xpv,zpv,squeeze(pv1(:,yind,:,tind))',20);colorbar
        title('(f + v_x -u_y)\rho_z');
        figure;
        contourf(xpv,zpv,squeeze(pv2(:,yind,:,tind))',20);colorbar
        title('v_z \rho_x');
        figure;
        contourf(xpv,zpv,squeeze(pv3(:,yind,:,tind))',20);colorbar
        title('u_x \rho_y');
        figure;
        contourf(xpv,zpv,squeeze(pv(:,yind,:,tind))',20);colorbar
        title('Full PV');
        colormap(hsv);
        pause;
    end
end
intPV = domain_integrate(pv,xpv,ypv,zpv);
save pv.mat pv xpv ypv zpv tpv intPV
fprintf('\n Wrote file : %s \n\n',outname);

function [um] = avgy(um)
    um = (um(:,1:end-1,:,:)+um(:,2:end,:,:))/2;

function [um] = avgx(um)
    um = (um(1:end-1,:,:,:)+um(2:end,:,:,:))/2;

function [um] = avgz(um)
    um = (um(:,:,1:end-1,:)+um(:,:,2:end,:))/2;

    %% old code
    
%     pv1    = avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  (tz(2:end-1,2:end-1,:,:));
%     pv2    = (-1)*;
%     pv3    = uz.*avgz(tx);
    %pv = double((pv1 + avgy(pv2(2:end-1,:,:,:)) + avgx(pv3(:,2:end-1,:,:)))./av