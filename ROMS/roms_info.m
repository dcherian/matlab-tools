% Prints out run info and parameters
%       [] = roms_info(fname,plot)
%            fname - filename
%            plot  - plot initial condition? (0/1)

function [] = roms_info(fname,plot)

    %fname = find_file(fname);
    if isdir(fname)
        fnames = roms_find_file(fname,'his');
        fname = [fname '/' char(fnames(:,1))];
    end
    if ~exist('plot','var')
        plot = 0;
    end
 
    if ~exist('fname','var')
        fname = 'ocean_his.nc';
    end
    
    grid = roms_get_grid(fname,fname,0,1);

    title_file  = ncreadatt(fname,'/','title');
    cpp    = ncreadatt(fname,'/','CPP_options');
    lbc    = ncreadatt(fname,'/','NLM_LBC');
    hist   = ncreadatt(fname,'/','history');
    ntimes = ncread(fname,'ntimes');
    dtfast = ncread(fname,'dtfast');
    dt     = ncread(fname,'dt');
    nhis   = ncread(fname,'nHIS');
    navg   = 0;%ncread(fname,'nAVG');
    ndia   = 0;%ncread(fname,'nDIA');
    nrst   = ncread(fname,'nRST');
    X      = ncread(fname,'xl');
    Y      = ncread(fname,'el');
    Z      = max(abs(grid.z_u(:)));
    f      = ncread(fname,'f');
    R0     = ncread(fname,'R0');
    TCOEF  = ncread(fname,'Tcoef');
    SCOEF  = ncread(fname,'Scoef');
    
    % viscosities
    rdrg   = ncread(fname,'rdrg');
    rdrg2  = ncread(fname,'rdrg2');
    ii = 0;
    if strfind(cpp,'UV_VIS2')
        tnu   = ncread(fname,'nl_tnu2');
        visc  = ncread(fname,'nl_visc2');
        ii = 2;
    end
    if strfind(cpp,'UV_VIS4')
        tnu   = ncread(fname,'nl_tnu4');
        visc  = ncread(fname,'nl_visc4');
        ii = 4;
    end
    if strfind(cpp,'UV_QDRAG')
        rdrg = 0;
    else
        rdrg2 = 0;
    end
    
    Vtransform  = ncread(fname,'Vtransform');
    Vstretching = ncread(fname,'Vstretching');
    theta_s     = ncread(fname,'theta_s');
    theta_b     = ncread(fname,'theta_b');
    Tcline      = ncread(fname,'Tcline');
    spherical   = ncread(fname,'spherical');
    
    % grids
    ocean_time = ncread(fname,'ocean_time');
    ocean_time = ocean_time/(24*3600);  
    
    Lm = size(grid.x_rho,2)-2;
    Mm = size(grid.x_rho,1)-2;
    N  = size(grid.z_r,1);
    
    xdiff = diff(grid.x_rho,1,2);
    ydiff = diff(grid.y_rho,1,1);
    zdiff = diff(grid.z_r,1,1);
    fdiff = diff(f,1,2);
        
    xdiff = [min(xdiff(:)) max(xdiff(:))];
    ydiff = [min(ydiff(:)) max(ydiff(:))];
    zdiff = [min(zdiff(:)) max(zdiff(:))];
    f0    = mean(f(:));
    beta  = mean(fdiff./diff(grid.y_rho',1,2));
    beta  = mean(beta(:));
    
    xunit = 'm'; yunit = 'm';
    if X > 1000, X = X/1000; xunit = 'km'; end
    if Y > 1000, Y = Y/1000; yunit = 'km'; end
    
    %% fix lbc    
    loc = strfind(lbc,char(10));
    lbcf = [''];
    for i=1:length(loc)-1
        lbcf = sprintf('%s \n\t\t\t   %s',lbcf, lbc(loc(i)+1:loc(i+1)-1));
    end
    lbcf = sprintf('%s \n\t\t\t   %s',lbcf, lbc(loc(end)+1:end));
    %%
    % pretty print stuff
    fprintf(['\n\t\t\t %40s \n', ...
             '\n\t\t\t %40s \n', ...
             '\n\t\t %40s \n', ...
             '\n\t\t %100s \n', ...
             '\n\t\t %10s: %6.2f %2s \t %3d \t (%6.2f : %6.2f) m     ', ...
             '\n\t\t %10s: %6.2f %2s \t %3d \t (%6.2f : %6.2f) m     ', ...
             '\n\t\t %10s: %6.2f %2s \t %3d \t (%6.2f : %6.2f) m     ', ...
             '\n\t\t %10s: (n = %d) %.2f d : %.2f s : %.2f d |  dtfast = %.4f s ', ...
             '\n', ...
             '\n\t\t %10s: [%d %d %d %d %.2f]  | %s: %d      ', ...
             '\n', ...
             '\n\t\t %10s: visc%d = %.2e  | tnu%d = [%.2e'], ...
             title_file,fname,hist, lbcf, ...
             'Domain: X',X,xunit,Lm,xdiff(1),xdiff(2), 'Y', Y,yunit,Mm, ydiff(1),ydiff(2), 'Z', Z,'m',N,zdiff(1),zdiff(2), ...
             'T',length(ocean_time),ocean_time(1),dt,ocean_time(end), dtfast,    ...
             'Stretching', Vtransform, Vstretching,theta_s,theta_b,Tcline, 'Spherical', spherical, ...
             'Viscosity', ii, visc(1),ii,tnu(1));
    for i=2:length(tnu)
       fprintf(' %.2e',tnu(i)); 
    end
    fprintf(']');
    fprintf(['\n\t\t %10s: rdrg  = %.2e  | rdrg2 = %.2e      ', ...
             '\n', ...
             '\n\t\t %10s:   R0 = %.2f   | TCOEF = %1.3e | SCOEF = %1.3e ', ...
             '\n\t\t %10s:    f = %1.2e  | beta  = %1.3e               ', ...
             '\n', ...
             '\n\t\t %10s: %4d = %5.2f hrs | %3s: %4d = %5.2f hrs ', ...
             '\n\t\t %10s: %4d = %5.2f hrs | %3s: %4d = %5.2f hrs ', ...
             '\n\n'], ...
             'Bott Fric', rdrg,rdrg2,...
             'EOS',R0, TCOEF, SCOEF,'Coriolis',f0,beta, ...
             'Output HIS',nhis,nhis*dt/3600,'DIA',ndia,ndia*dt/3600, 'AVG',navg,navg*dt/3600,'RST',nrst,double(nrst)*dt/3600);
    
    %% plot initial condition
    
    if plot        
        % read data
        u    = ncread(fname,'u',[1 1 1 1],[Inf Inf Inf 1]);
        v    = ncread(fname,'v',[1 1 1 1],[Inf Inf Inf 1]);
        temp = ncread(fname,'temp',[1 1 1 1],[Inf Inf Inf 1]);
        zeta = ncread(fname,'zeta',[1 1 1],[Inf Inf 1]);
        
        ymid = Mm/2;
        zmid = N/2;
        
        % plot
        figure;
        subplot(241)
        contourf(squeeze(grid.x_u(1,:)),squeeze(grid.z_u(:,1,1)),squeeze(u(:,ymid,:))');
        colorbar;
        title('u');
        xlabel('x'); ylabel('z');

        subplot(242)
        contourf(grid.x_v(1,:),squeeze(grid.z_v(:,1,1)),squeeze(v(:,ymid,:))',20);
        colorbar;
        title('v');
        xlabel('x'); ylabel('z');

        subplot(243)
        contourf(grid.x_rho(1,:),squeeze(grid.z_r(:,1,1)),squeeze(temp(:,ymid,:))',20);
        colorbar;
        title('temp');
        xlabel('x'); ylabel('z');

        subplot(245)
        contourf(grid.x_u(1,:),grid.y_u(:,1),squeeze(u(:,:,zmid))');
        colorbar;
        title('u');
        xlabel('x'); ylabel('y');

        subplot(246)
        contourf(grid.x_v(1,:),grid.y_v(:,1),squeeze(v(:,:,zmid))');
        colorbar;
        title('v');
        xlabel('x'); ylabel('y');

        subplot(247)
        contourf(grid.x_rho(1,:),grid.y_rho(:,1),squeeze(temp(:,:,zmid))');
        colorbar;
        title('temp');
        xlabel('x'); ylabel('y');

        subplot(248)
        contourf(grid.x_rho(1,:),squeeze(grid.y_rho(:,1)),squeeze(zeta(:,:,1))');
        colorbar;
        title('zeta');
        xlabel('x'); ylabel('y');
    end