% Warning: this doesn't do variations in y very well.
%       [APE] = roms_ape(fname)

function [APE] = roms_ape(fname)
    warning off
    roms_grid = roms_get_grid(fname,fname,0,1);
    warning on
    
    rho = ncread(fname,'rho',[1 1 1 1],[Inf Inf Inf 1]);
    
    APE = 0;
    ny = size(roms_grid.x_rho,1);
    for i=1:ny
        APE = APE + calc_ape(roms_grid.x_psi(1,:),roms_grid.z_w(:,1,1),rho(:,i,:));
    end
    
    APE = APE./ny;
    