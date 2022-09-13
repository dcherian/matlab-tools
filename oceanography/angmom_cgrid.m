% returns angular momentum calculated at RHO points
%   [angmom] = angmom_cgrid(grid,u,v,rho,center)

function [angmom,TAM] = angmom_cgrid(grid,u,v,rho,center)

    x = grid.xrmat - center.x;
    y = grid.yrmat - center.y;
    z = grid.zrmat - center.z;
    
    
    % ignoring w
    angmom_x = bsxfun(@times, avg1(- v(2:end-1,:,:,:),2), z(2:end-1,2:end-1,:));
    angmom_y = bsxfun(@times, avg1(u(:,2:end-1,:,:,:),1), z(2:end-1,2:end-1,:));
    angmom_z = bsxfun(@times, avg1(v(2:end-1,:,:,:),2), x(2:end-1,2:end-1,:)) ...
               - bsxfun(@times,avg1(u(:,2:end-1,:,:,:),1), y(2:end-1,2:end-1,:));
    
    angmom = rho(2:end-1,2:end-1,:,:) .* sqrt(angmom_x.^2 + angmom_y.^2 + angmom_z.^2);
    
    % domain integrate
    xax = squeeze(grid.xrmat(2:end-1,1,1));
    yax = squeeze(grid.yrmat(1,2:end-1,1));
    zax = squeeze(grid.zrmat(1,1,:));

    TAM = domain_integrate(angmom,xax,yax,zax); 