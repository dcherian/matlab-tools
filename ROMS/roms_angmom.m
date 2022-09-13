% returns angular momentum calculated about (x0,y0,z0) at RHO points
%   [angmom,TAM] = roms_angmom(fname,x0,y0,z0,tindices,vol)
% TAM is domain integrated value


function [angmom,TAM] = roms_angmom(fname,x0,y0,z0,tindices,vol)

    grid1 = dc_roms_get_grid(fname);
    
    center.x = x0;
    center.y = y0;
    center.z = z0;
    
    u = double(ncread(fname,'u'));
    v = double(ncread(fname,'v'));
    rho = double(ncread(fname,'rho'));
    
    [angmom,TAM] = angmom_cgrid(grid1,u,v,rho,center);