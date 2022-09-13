function [] = gcm_info(fname)

    fname = find_file(fname);

    vars = rdmnc(fname);
    
    title  = ncreadatt(fname,'/','title');
    lbc    = ncreadatt(fname,'/','NLM_LBC');
    hist   = ncreadatt(fname,'/','history');
    ntimes = ncread(fname,'ntimes');
    dt     = ncread(fname,'dt');
    nhis   = ncread(fname,'nHIS');
    navg   = ncread(fname,'nAVG');
    ndia   = ncread(fname,'nDIA');
    nrst   = ncread(fname,'nRST');
    X      = ncread(fname,'xl');
    Y      = ncread(fname,'el');
    Z      = max(abs(grid.z_u(:)));
    f      = ncread(fname,'f');
    R0     = ncread(fname,'R0');
    TCOEF  = ncread(fname,'Tcoef');
    SCOEF  = ncread(fname,'Scoef');
    
    Vtransform  = ncread(fname,'Vtransform');
    Vstretching = ncread(fname,'Vstretching');
    theta_s     = ncread(fname,'theta_s');
    theta_b     = ncread(fname,'theta_b');
    Tcline      = ncread(fname,'Tcline');
    spherical   = ncread(fname,'spherical');
    
    % grids
    ocean_time = ncread(fname,'T');
    ocean_time = ocean_time/(24*3600);  
    
    Lm = size(grid.x_rho,2)-2;
    Mm = size(grid.x_rho,1)-2;
    N  = size(grid.z_r,1);
    
    xdiff = diff(vars.X,1,2);
    ydiff = diff(vars.Y,1,1);
    zdiff = diff(vars.Z,1,1);
    %fdiff = diff(f,1,2);
        
    xdiff = [min(xdiff(:)) max(xdiff(:))];
    ydiff = [min(ydiff(:)) max(ydiff(:))];
    zdiff = [min(zdiff(:)) max(zdiff(:))];
    f0    = mean(f(:));
    beta  = mean(fdiff(:))./ydiff(1);
    
    xunit = 'm'; yunit = 'm';
    if X > 1000, X = X/1000; xunit = 'km'; end
    if Y > 1000, Y = Y/1000; yunit = 'km'; end
    
    fprintf(['\n\t\t\t %40s \n', ...
            '\n\t\t %40s \n', ...
            '\n\t\t %100s \n', ...
            '\n\t\t %10s: %6.2f %2s \t %3d \t (%6.2f : %6.2f) m     ', ...
            '\n\t\t %10s: %6.2f %2s \t %3d \t (%6.2f : %6.2f) m     ', ...
            '\n\t\t %10s: %6.2f %2s \t %3d \t (%6.2f : %6.2f) m     ', ...
            '\n\t\t %10s: %.2f d : %.2f s : %.2f d |  dtfast = %.4f s ', ...
            '\n', ...
            '\n\t\t %10s:   R0 = %.2f   | TCOEF = %1.3e | SCOEF = %1.3e ', ...
            '\n\t\t %10s:    f = %1.2e | beta  = %1.3e               ', ...
            '\n', ...
            '\n\t\t %10s: %4d = %5.2f hrs', ...
            '\n\n'], ...
            title,hist, lbcf, ...
            'Domain: X',X,xunit,Lm,xdiff, 'Y', Y,yunit,Mm, ydiff, 'Z', Z,'m',N,zdiff(1),zdiff(2), ...
            'T',ocean_time(1),dt,ocean_time(end), dtfast,    ...
            'EOS',R0, TCOEF, SCOEF,'Coriolis',f0,beta, ...
            'Output HIS',nhis,nhis*dt/3600);